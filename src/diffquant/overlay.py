"""Visual check of the calls: ImageJ TIFFs (all channels, background-subtracted) with cell outlines coloured by call.

Outline colours: green = readout+ (gated), magenta = readout- (gated), cyan = excluded (below the gate or removed by filters).
Open in Fiji (Image > Color > Channels Tool to switch channels; the overlay shows the outlines).
"""
from __future__ import annotations

import os
import random

import numpy as np
import pandas as pd
import tifffile
from skimage import measure as skm

from . import io

COL = {"pos": b"\xff\x00\xff\x00", "neg": b"\xff\xff\x00\xff", "low": b"\xff\x00\xff\xff"}   # ARGB


def run(cfg: dict, plate: dict, n_sites: int = 3, sites=None, marker: str | None = None, seed: int = 1, crop: int | None = 1080) -> list:
    import roifile
    mk = marker or cfg["readouts"][0]["marker"]
    od = os.path.join(cfg["output"], plate["name"])
    cells = pd.read_csv(os.path.join(od, "cells_called.csv.gz"))
    files = io.find_sites(plate, cfg)
    masks = io.find_masks(cfg, plate)
    keys = [(w, int(s)) for w, s in sites] if sites else random.Random(seed).sample(
        [k for k in files if k in masks], min(n_sites, len([k for k in files if k in masks])))
    outd = os.path.join(od, "overlays")
    os.makedirs(outd, exist_ok=True)
    written = []
    for well, site in keys:
        c = cells[(cells.well == well) & (cells.site == site)]
        m = io.read_mask(masks[(well, site)])
        H = min(crop or m.shape[0], m.shape[0])
        y0, x0 = (m.shape[0] - H) // 2, (m.shape[1] - H) // 2
        mc = m[y0:y0 + H, x0:x0 + H]
        order = [mk] + [x for x in cfg["channels"] if x != mk]
        chans = [np.clip(io.read_image(files[(well, site)][x])[y0:y0 + H, x0:x0 + H] - c[f"{x}_background"].iloc[0], 0, None)
                 for x in order]
        cls = {int(r.label): ("low" if (r.gate_low or r.get("filtered", False)) else ("pos" if r[f"{mk}_pos"] else ("neg" if r.gated else None)))
               for _, r in c.iterrows()}
        rois = []
        for p in skm.regionprops(mc):
            k = cls.get(p.label)
            if k is None:
                continue
            cs = skm.find_contours(np.pad(p.image, 1).astype(float), 0.5)
            if not cs:
                continue
            cm = max(cs, key=len) - 1 + np.array(p.bbox[:2])
            roi = roifile.ImagejRoi.frompoints(cm[:, ::-1] + 0.5)
            roi.roitype = roifile.ROI_TYPE.POLYGON
            roi.stroke_color = COL[k]
            roi.stroke_width = 2
            rois.append(roi.tobytes())
        f = os.path.join(outd, f"overlay_{well}_s{site}.tif")
        tifffile.imwrite(f, np.stack(chans).astype(np.float32), imagej=True,
                         metadata={"axes": "CYX", "mode": "color", "Overlays": rois, "Labels": order})
        written.append(f)
        print(f"[overlay] {plate['name']} {well} s{site}: {sum(v == 'pos' for v in cls.values())} {mk}+, "
              f"{sum(v == 'neg' for v in cls.values())} {mk}-, {sum(v == 'low' for v in cls.values())} below gate -> {f}")
    return written
