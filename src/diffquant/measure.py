"""Per-cell mean intensities and background subtraction.

Background (per site and channel) = median of cell-free pixels, i.e. pixels farther than ``far_px`` from any mask.
A site uses its own value only if cell-free pixels cover >= ``min_cellfree_frac`` of the field; dense sites use the plate
value = median over qualifying sites (if fewer than ``min_sites`` qualify: the ``min_sites`` sites with the largest
cell-free fraction). Background-subtracted cell mean = max(raw mean - background, 0).
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure as skm

from . import io
from ._parallel import pmap


def _site(args):
    well, site, mpath, files, far_px = args
    m = io.read_mask(mpath)
    H, W = m.shape
    far = ndi.distance_transform_edt(m == 0) > far_px
    props = skm.regionprops(m)
    lab = np.array([p.label for p in props], int)
    df = pd.DataFrame(dict(well=well, site=site, label=lab, area_px=[p.area for p in props],
                           centroid_y=[round(p.centroid[0], 1) for p in props], centroid_x=[round(p.centroid[1], 1) for p in props],
                           solidity=[round(p.solidity, 4) for p in props],
                           touches_border=[p.bbox[0] == 0 or p.bbox[1] == 0 or p.bbox[2] == H or p.bbox[3] == W for p in props]))
    bg = dict(well=well, site=site, frac_cellfree=round(float(far.mean()), 4))
    if len(lab):
        idx = m.ravel()
        n = np.bincount(idx)[lab]
    for marker, path in files.items():
        im = io.read_image(path)
        bg[f"{marker}_bg_cellfree"] = float(np.median(im[far])) if far.sum() > 1000 else np.nan
        if len(lab):
            df[f"{marker}_mean_raw"] = np.round(np.bincount(idx, weights=im.ravel())[lab] / n, 3)
    return df, bg


def site_backgrounds(bgs: pd.DataFrame, markers, b: dict):
    """Apply the background rule to the per-site table; returns (table with <marker>_background, plate values, rule)."""
    bgs = bgs.copy()
    ok = bgs.frac_cellfree >= b["min_cellfree_frac"]
    rule = f">= {b['min_cellfree_frac']:.0%} cell-free"
    if ok.sum() < b["min_sites"]:
        ok = bgs.frac_cellfree >= bgs.frac_cellfree.nlargest(b["min_sites"]).min()
        rule = f"top {b['min_sites']} sites by cell-free fraction"
    plate_bg = {}
    for mk in markers:
        col = f"{mk}_bg_cellfree"
        plate_bg[mk] = float(bgs.loc[ok & bgs[col].notna(), col].median())
        bgs[f"{mk}_background"] = np.where(ok & bgs[col].notna(), bgs[col], plate_bg[mk])
    bgs["bg_source"] = np.where(ok, "site", "plate_median")
    return bgs, plate_bg, rule


def run(cfg: dict, plate: dict) -> pd.DataFrame:
    sites = io.find_sites(plate, cfg)
    masks = io.find_masks(cfg, plate)
    b = cfg["background"]
    args = [(w, s, masks[(w, s)], files, b["far_px"]) for (w, s), files in sites.items() if (w, s) in masks]
    if not args:
        raise SystemExit(f"[measure] {plate['name']}: no site has both images and a mask (run 'diffquant segment' first)")
    res = pmap(_site, args, int(cfg.get("jobs", 4)), 4)
    cells = pd.concat([r[0] for r in res], ignore_index=True)
    bgs = pd.DataFrame([r[1] for r in res])
    markers = list(cfg["channels"])
    od = os.path.join(cfg["output"], plate["name"])
    if b.get("enabled", True):
        bgs, plate_bg, rule = site_backgrounds(bgs, markers, b)
    else:
        plate_bg, rule = {m: 0.0 for m in markers}, "disabled"
        for mk in markers:
            bgs[f"{mk}_background"] = 0.0
        bgs["bg_source"] = "none"
    bgs.to_csv(os.path.join(od, "site_background.csv"), index=False)
    keep = ["well", "site", "frac_cellfree", "bg_source"] + [f"{mk}_background" for mk in markers]
    cells = cells.merge(bgs[keep], on=["well", "site"], how="left")
    for mk in markers:
        cells[f"{mk}_background"] = cells[f"{mk}_background"].round(3)
        cells[f"{mk}_mean"] = np.maximum(cells[f"{mk}_mean_raw"] - cells[f"{mk}_background"], 0).round(3)
    cells.to_csv(os.path.join(od, "cells.csv.gz"), index=False, compression="gzip")
    with open(os.path.join(od, "measure_info.json"), "w") as fh:
        json.dump(dict(plate=plate["name"], n_sites=len(args), n_sites_without_mask=len(sites) - len(args), n_cells=len(cells),
                       background=dict(b, rule=rule, plate_values={k: round(v, 1) for k, v in plate_bg.items()},
                                       n_sites_plate_value=int((bgs.bg_source != "site").sum())),
                       mean="max(raw mean - background, 0) per cell"), fh, indent=1)
    print(f"[measure] {plate['name']}: {len(args)} sites, {len(cells)} cells, background {rule}")
    return cells
