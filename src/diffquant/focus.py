"""Focus QC per site: normalised variance of the Laplacian.

focus = var(Laplacian(Gaussian(image, sigma = 1 px))) / median(image)^2 on a 2x down-sampled image (scale-free; lower = blurrier).
Sites below 0.5 x the plate median are flagged (``focus_flag``); whole plates with many flagged sites are worth checking by eye.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from . import io
from ._parallel import pmap


def focus_value(path: str) -> float:
    im = io.read_image(path)[::2, ::2]
    lap = ndi.laplace(ndi.gaussian_filter(im, 1.0))
    return float(lap.var() / max(np.median(im), 1.0) ** 2)


def _job(a):
    well, site, paths = a
    return dict(well=well, site=site, **{f"focus_{k}": focus_value(p) for k, p in paths.items()})


def run(cfg: dict, plate: dict) -> pd.DataFrame:
    chans = cfg["focus"].get("channels") or list(cfg["channels"])
    sites = io.find_sites(plate, cfg)
    args = [(w, s, {m: f[m] for m in chans}) for (w, s), f in sites.items()]
    df = pd.DataFrame(pmap(_job, args, int(cfg.get("jobs", 4)), 8))
    flag = np.zeros(len(df), bool)
    for m in chans:
        rel = df[f"focus_{m}"] / df[f"focus_{m}"].median()
        df[f"focus_{m}_rel"] = rel.round(3)
        flag |= (rel < 0.5).values
    df["focus_flag"] = flag
    df.insert(0, "plate", plate["name"])
    df.to_csv(os.path.join(cfg["output"], plate["name"], "focus.csv"), index=False)
    print(f"[focus] {plate['name']}: {len(df)} sites, {int(flag.sum())} flagged (< 0.5 x plate median)")
    return df
