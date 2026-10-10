"""Finding images and masks on disk."""
from __future__ import annotations

import os
import re

import numpy as np
import tifffile


def channels(cfg: dict, plate: dict) -> dict:
    """marker -> file-name channel token; a plate entry may override the global mapping (``channels:`` under the plate)."""
    return dict(cfg["channels"], **(plate.get("channels") or {}))


def find_sites(plate: dict, cfg: dict) -> dict:
    """{(well, site): {marker: file}} for every site that has all configured channels."""
    im = cfg["images"]
    plate_path = plate["path"]
    pat = re.compile(im["pattern"])
    chmap = channels(cfg, plate)
    token2marker = {v: k for k, v in chmap.items()}
    found = {}
    for f in sorted(os.listdir(plate_path)):
        if not f.lower().endswith(tuple(im["extensions"])) or any(s in f for s in im["ignore"]):
            continue
        m = pat.search(f)
        if not m or m["channel"] not in token2marker:
            continue
        found.setdefault((m["well"], int(m["site"])), {})[token2marker[m["channel"]]] = os.path.join(plate_path, f)
    need = set(chmap)
    return {k: v for k, v in sorted(found.items()) if need <= set(v)}


def read_image(path: str) -> np.ndarray:
    return tifffile.imread(path).astype(np.float32)


def mask_path(cfg: dict, plate: str, well: str, site: int) -> str:
    return os.path.join(cfg["output"], plate, "masks", f"{well}_s{site}_masks.tif")


def find_masks(cfg: dict, plate: dict) -> dict:
    """{(well, site): mask file}. Uses segmentation.masks (existing masks, e.g. PNG from an earlier run) if set,
    else the masks written by ``diffquant segment``."""
    seg = cfg["segmentation"]
    out = {}
    if seg.get("masks"):
        d = seg["masks"]
        d = os.path.join(d, plate["name"]) if os.path.isdir(os.path.join(d, plate["name"])) else d
        pat = re.compile(seg["mask_pattern"])
        for f in sorted(os.listdir(d)):
            m = pat.search(f)
            if m and f.lower().endswith((".tif", ".tiff", ".png")):
                out[(m["well"], int(m["site"]))] = os.path.join(d, f)
        return out
    d = os.path.join(cfg["output"], plate["name"], "masks")
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            m = re.match(r"([A-P]\d{2})_s(\d+)_masks\.tif$", f)
            if m:
                out[(m[1], int(m[2]))] = os.path.join(d, f)
    return out


def read_mask(path: str) -> np.ndarray:
    if path.lower().endswith(".png"):
        try:
            import imageio.v3 as iio
            return np.asarray(iio.imread(path)).astype(np.int32)
        except ImportError:                       # older imageio
            from skimage import io as skio
            return np.asarray(skio.imread(path)).astype(np.int32)
    return tifffile.imread(path).astype(np.int32)


def write_mask(path: str, labels: np.ndarray) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    dt = np.uint16 if labels.max() < 65535 else np.uint32
    tmp = path + ".tmp.tif"
    tifffile.imwrite(tmp, labels.astype(dt), compression="zlib")
    os.replace(tmp, path)


def inspect(plate: dict, cfg: dict, n_lines: int = 6) -> list:
    """One line per channel token found in the plate folder: example file + wavelength / filter lines from the TIFF
    description (MetaMorph writes these), to decide which token holds which marker."""
    im = cfg["images"]
    pat = re.compile(im["pattern"])
    first = {}
    for f in sorted(os.listdir(plate["path"])):
        if f.lower().endswith(tuple(im["extensions"])) and not any(s in f for s in im["ignore"]):
            m = pat.search(f)
            if m and m["channel"] not in first:
                first[m["channel"]] = f
    out = []
    for tok, f in sorted(first.items()):
        with tifffile.TiffFile(os.path.join(plate["path"], f)) as t:
            desc = (t.pages[0].description or "").replace("\r", "\n")
        keep = [l.strip() for l in desc.split("\n") if re.search(r"(?i)wave|filter|illum|exposure|emission|excitation|_name", l)]
        out.append(f"{tok}: {f} | " + (" | ".join(keep[:n_lines]) or "no metadata"))
    return out
