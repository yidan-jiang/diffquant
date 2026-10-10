"""Blinded single-cell scoring sets for calibration (scored in Fiji with fiji/score_cells.py).

For every plate and marker a stack of 160 x 160 px crops (background-subtracted marker channel, one cell per slice, outline as
an ImageJ overlay, fixed display range) is written in random order, with no intensity values shown to the scorer.
  readout set: gated cells, stratified by readout-intensity decile (``per_decile`` each) + ``n_random`` unstratified cells
               (at most one per site); defaults 6 x 10 + 40 = 100 cells per plate
  gate set:    all cells, stratified by distance below the gate peak (5 log10 bands, ``per_band`` each; default 60 cells)
A few cells are shown twice (``repeats``) to measure scorer consistency.
Output: <output>/scoring/<plate>_<marker>/crops.tif, index.csv (hidden key), display.txt.
"""
from __future__ import annotations

import os
import zlib

import numpy as np
import pandas as pd
import tifffile
from skimage import measure as skm

from . import io
from .config import gate_marker
from .threshold import log10p

B = 160
GATE_BANDS = [-9, -1.2, -0.75, -0.55, -0.4, -0.1]      # log10 below the gate peak


def _seed(name, k=0):
    return (zlib.crc32(str(name).encode()) + k) % (2 ** 31)


def run(cfg: dict, plate: dict, per_decile: int = 6, n_random: int = 40, per_band: int = 12,
        repeats_readout: int = 10, repeats_gate: int = 8) -> None:
    import roifile
    name = plate["name"]
    od = os.path.join(cfg["output"], name)
    cells = pd.read_csv(os.path.join(od, "cells_called.csv.gz"))
    thr = pd.read_csv(os.path.join(od, "thresholds.csv")).iloc[0]
    sites = io.find_sites(plate, cfg)
    masks = io.find_masks(cfg, plate)
    cells = cells[cells.area_px >= 0.25 * cells.area_px.median()]
    cells = cells[~cells.touches_border.astype(bool)]
    plans = {}
    gated = cells[cells.gated.astype(bool)]
    for r in cfg["readouts"]:
        mk = r["marker"]
        c = gated.assign(stratum=pd.qcut(log10p(gated[f"{mk}_mean"]), 10, labels=False, duplicates="drop"))
        plans[mk] = (c, per_decile, n_random, repeats_readout)
    gm = gate_marker(cfg)
    if gm:
        rel = log10p(cells[f"{gm}_mean"]) - np.log10(thr.gate_peak + 1)
        c = cells.assign(stratum=pd.cut(rel, GATE_BANDS, labels=False)).dropna(subset=["stratum"])
        plans[gm] = (c, per_band, 0, repeats_gate)
    cache = {}

    def get(well, site, mk):
        k = (well, site, mk)
        if k not in cache:
            if len(cache) > 40:
                cache.clear()
            cache[k] = (io.read_mask(masks[(well, site)]), io.read_image(sites[(well, site)][mk]))
        return cache[k]

    def fits(r, mk):
        if (r.well, int(r.site)) not in masks or (r.well, int(r.site)) not in sites:
            return False
        m, _ = get(r.well, int(r.site), mk)
        rr, cc = np.nonzero(m == int(r.label))
        return len(rr) > 0 and rr.max() - rr.min() <= B - 8 and cc.max() - cc.min() <= B - 8

    for mk, (cand, per, nrand, nrep) in plans.items():
        acc = []
        for st, grp in cand.groupby("stratum"):
            got = 0
            for _, r in grp.sample(frac=1, random_state=_seed(name)).head(per * 6).iterrows():
                if fits(r, mk):
                    acc.append(r)
                    got += 1
                    if got == per:
                        break
        if nrand:
            taken = {(r.well, int(r.site), int(r.label)) for r in acc}
            seen, got = set(), 0
            for _, r in cand.sample(frac=1, random_state=_seed(name, 7)).iterrows():
                k = (r.well, int(r.site))
                if k in seen or (k[0], k[1], int(r.label)) in taken or not fits(r, mk):
                    continue
                acc.append(r.copy().rename(None))
                acc[-1]["stratum"] = -1
                seen.add(k)
                got += 1
                if got == nrand:
                    break
        df = pd.DataFrame(acc).reset_index(drop=True)
        df["repeat"] = 0
        df = pd.concat([df, df.sample(min(nrep, len(df)), random_state=_seed(name, 1)).assign(repeat=1)])
        df = df.sample(frac=1, random_state=_seed(name, 2)).reset_index(drop=True)
        stack, rois, pix = [], [], []
        for i, r in df.iterrows():
            m, im = get(r.well, int(r.site), mk)
            cell = m == int(r.label)
            rr, cc = np.nonzero(cell)
            cy, cx = (rr.min() + rr.max()) // 2, (cc.min() + cc.max()) // 2
            y0 = int(np.clip(cy - B // 2, 0, m.shape[0] - B))
            x0 = int(np.clip(cx - B // 2, 0, m.shape[1] - B))
            crop = np.clip(im[y0:y0 + B, x0:x0 + B] - r[f"{mk}_background"], 0, None)
            sub = cell[y0:y0 + B, x0:x0 + B]
            stack.append(crop.astype(np.float32))
            pix.append(crop[sub])
            cont = max(skm.find_contours(np.pad(sub, 1).astype(float), 0.5), key=len) - 1
            roi = roifile.ImagejRoi.frompoints(cont[:, ::-1] + 0.5)
            roi.roitype = roifile.ROI_TYPE.POLYGON
            roi.position = i + 1
            roi.stroke_color = b"\xff\xff\xff\x00"
            roi.stroke_width = 1
            rois.append(roi.tobytes())
        sd = os.path.join(cfg["output"], "scoring", f"{name}_{mk}")
        os.makedirs(sd, exist_ok=True)
        tifffile.imwrite(os.path.join(sd, "crops.tif"), np.stack(stack), imagej=True, metadata={"axes": "ZYX", "Overlays": rois})
        df.insert(0, "slice", np.arange(1, len(df) + 1))
        keep = ["slice", "well", "site", "label", "repeat", "stratum"] + [f"{m}_mean" for m in cfg["channels"] if f"{m}_mean" in df]
        df[keep].to_csv(os.path.join(sd, "index.csv"), index=False)
        vmax = float(np.percentile(np.concatenate(pix), 90))
        with open(os.path.join(sd, "display.txt"), "w") as fh:
            fh.write(f"# {name} {mk}: background-subtracted, linear; black = 0, white = 90th percentile of the shown cells\n"
                     f"marker {mk}\nrole {'gate' if mk == gm else 'readout'}\nmin 0\nmax {vmax:.1f}\n")
        print(f"[score-sets] {name} {mk}: {len(df)} crops ({int(df.repeat.sum())} repeats) -> {sd}")
