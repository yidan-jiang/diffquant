"""Configuration: one YAML file describes the plates, channels, and every rule of the pipeline.

Relative paths in the YAML are resolved against the folder that contains the YAML file.
Run ``diffquant init`` to write a commented template.
"""
from __future__ import annotations

import copy
import os

import yaml

DEFAULTS = {
    "output": "results",
    "plates": [],                      # [{name, path[, channels]}] - one folder of images per plate
    "images": {
        # ImageXpress / MetaMorph naming: <prefix>_<well>_s<site>_w<n><GUID>.tif
        "pattern": r"_(?P<well>[A-P]\d{2})_s(?P<site>\d+)_(?P<channel>w\d)",
        "ignore": ["_thumb"],
        "extensions": [".tif", ".tiff"],
        "pixel_size_um": None,
    },
    "channels": {},                    # marker name -> channel token in the file name, e.g. {DAPI: w1, KRT14: w3, IVL: w2}
    "segmentation": {
        "model": "cpsam",              # any Cellpose model name or a path to a custom-trained model
        "channels": [],                # markers given to Cellpose, e.g. [KRT14, DAPI] (cytoplasm, nucleus)
        "gpu": "auto",
        "diameter": None,
        "flow_threshold": 0.4,
        "cellprob_threshold": 0.0,
        "masks": None,                 # folder with existing label masks (skips segmentation)
        "mask_pattern": r"_(?P<well>[A-P]\d{2})_s(?P<site>\d+)_",
    },
    "background": {
        "enabled": True,
        "far_px": 15,                  # cell-free pixel = farther than this from any mask
        "min_cellfree_frac": 0.20,     # a site uses its own background only if >= this fraction is cell-free
        "min_sites": 10,               # otherwise: plate median over qualifying sites (top-N sites if too few qualify)
    },
    "filters": {                       # object filter before gating (default off = all Cellpose objects are kept)
        "min_area_px": None,           # drop objects smaller than this (debris, fragments)
        "max_area_px": None,           # drop objects larger than this (merged clumps)
        "exclude_border": False,       # drop cells cut by the image edge
    },
    "gate": {                          # cell-type gate; set "gate: null" to score all segmented cells
        "marker": None,
        "rule": "peak_minus",          # peak_minus | otsu | fixed
        "delta": 0.8,                  # peak_minus: log10 units below the main (cell-type) peak
        "value": None,                 # fixed: intensity cut (background-subtracted)
        "upper_mad": 4.0,              # drop very bright outliers > median + k * MAD of gated cells (null = off)
    },
    "readouts": [],                    # [{marker, rule: otsu_times | fixed, factor, value}]
    "min_gated_cells": 200,            # wells with fewer gated cells are flagged invalid
    "focus": {"enabled": True, "channels": []},
    "plate_map": None,                 # optional CSV: plate, well, condition[, vehicle]
    "jobs": 4,
}

HIST = {"lo": -0.1, "hi": 4.6, "bins": 470, "sigma": 3}   # smoothed log10 histogram used to find the gate peak


def _merge(base: dict, new: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (new or {}).items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load(path: str) -> dict:
    with open(path) as fh:
        user = yaml.safe_load(fh) or {}
    cfg = _merge(DEFAULTS, user)
    if user.get("gate", "missing") is None:
        cfg["gate"] = None
    root = os.path.dirname(os.path.abspath(path))
    rel = lambda p: p if p is None or os.path.isabs(p) else os.path.normpath(os.path.join(root, p))
    cfg["output"] = rel(cfg["output"])
    cfg["plate_map"] = rel(cfg["plate_map"])
    cfg["segmentation"]["masks"] = rel(cfg["segmentation"]["masks"])
    for p in cfg["plates"]:
        p["name"] = str(p["name"])
        p["path"] = rel(p["path"])
    validate(cfg)
    cfg["_config_path"] = os.path.abspath(path)
    return cfg


def validate(cfg: dict) -> None:
    if not cfg["plates"]:
        raise ValueError("config: 'plates' is empty")
    if not cfg["channels"]:
        raise ValueError("config: 'channels' must map marker names to file-name channel tokens")
    markers = set(cfg["channels"])
    for p in cfg["plates"]:
        if set(p.get("channels") or {}) - markers:
            raise ValueError(f"config: plate {p['name']} overrides unknown markers {sorted(set(p['channels']) - markers)}")
    used = [r["marker"] for r in cfg["readouts"]] + list(cfg["segmentation"]["channels"])
    if cfg["gate"] is not None and cfg["gate"].get("marker"):
        used.append(cfg["gate"]["marker"])
    missing = [m for m in used if m not in markers]
    if missing:
        raise ValueError(f"config: markers {missing} are used but not listed under 'channels'")
    if not cfg["readouts"]:
        raise ValueError("config: at least one readout marker is needed")
    for r in cfg["readouts"]:
        if r.get("rule", "otsu_times") not in ("otsu_times", "fixed"):
            raise ValueError(f"config: unknown readout rule {r.get('rule')}")
    if cfg["gate"] is not None and cfg["gate"].get("marker") and cfg["gate"]["rule"] not in ("peak_minus", "otsu", "fixed"):
        raise ValueError(f"config: unknown gate rule {cfg['gate']['rule']}")


def gate_marker(cfg: dict):
    g = cfg.get("gate")
    return g["marker"] if g and g.get("marker") else None


def plate_out(cfg: dict, plate: str) -> str:
    d = os.path.join(cfg["output"], str(plate))
    os.makedirs(d, exist_ok=True)
    return d


TEMPLATE = """\
# diffquant configuration (example: involucrin+ % among KRT14+ keratinocytes)
# Relative paths are relative to this file.

output: results

plates:                      # one folder of images per plate
  - name: plate01
    path: /path/to/plate01
  - name: plate02
    path: /path/to/plate02
    channels: {IVL: w3, KRT14: w2}   # optional per-plate override if the wavelength order differs ('diffquant inspect')

images:
  # ImageXpress / MetaMorph naming <prefix>_<well>_s<site>_w<n>....tif ; adapt the regex for other microscopes
  pattern: '_(?P<well>[A-P]\\d{2})_s(?P<site>\\d+)_(?P<channel>w\\d)'
  ignore: ['_thumb']
  pixel_size_um: 0.648

channels:                    # marker name: channel token in the file name
  DAPI: w1
  IVL: w2
  KRT14: w3

segmentation:
  model: cpsam               # any Cellpose model that segments the cells reliably (cpsam, cyto3, or a custom model path)
  channels: [KRT14, DAPI]    # cytoplasm/cell marker, nucleus
  gpu: auto
  masks: null                # or a folder with existing label masks (skips segmentation)

background:
  far_px: 15
  min_cellfree_frac: 0.20
  min_sites: 10

filters:                     # optional object filter before gating (debris, clumps, edge cells); off by default
  min_area_px: null          # e.g. 200 (check the area histogram in cells.csv.gz first)
  max_area_px: null
  exclude_border: false

gate:                        # cell-type gate (set 'gate: null' to score all cells)
  marker: KRT14
  rule: peak_minus           # peak_minus | otsu | fixed
  delta: 0.8                 # default from the example set-up - recalibrate (docs/recalibration.md)
  upper_mad: 4

readouts:                    # one or more markers scored among gated cells
  - marker: IVL
    rule: otsu_times         # otsu_times | fixed
    factor: 1.3              # default from the example set-up - recalibrate (docs/recalibration.md)

min_gated_cells: 200

focus:
  channels: [DAPI, IVL]

plate_map: null              # optional CSV with columns plate, well, condition[, vehicle]
jobs: 4
"""
