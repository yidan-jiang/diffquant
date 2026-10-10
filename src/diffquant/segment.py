"""Cell segmentation with Cellpose (optional dependency: pip install "diffquant[seg]").

Any Cellpose model can be used; choose the one that segments the cells most reliably (inspect a few sites visually).
The example uses cpsam (Cellpose 4) on the cell-type marker (KRT14) + DAPI.
"""
from __future__ import annotations

import json
import os

import numpy as np

from . import io


def _model(seg: dict):
    try:
        from cellpose import models
    except ImportError as e:                      # pragma: no cover
        raise SystemExit("Cellpose is not installed: pip install 'diffquant[seg]'  (or provide segmentation.masks)") from e
    gpu = seg.get("gpu", "auto")
    if gpu == "auto":
        try:
            from cellpose import core
            gpu = bool(core.use_gpu())
        except Exception:
            gpu = False
    name = seg.get("model", "cpsam")
    try:                                          # Cellpose >= 4 (cpsam) and custom model paths
        return models.CellposeModel(gpu=bool(gpu), pretrained_model=name)
    except Exception:                             # Cellpose 3 built-in names (cyto3, nuclei, ...)
        return models.CellposeModel(gpu=bool(gpu), model_type=name)


def run(cfg: dict, plate: dict, overwrite: bool = False) -> int:
    seg = cfg["segmentation"]
    if seg.get("masks"):
        print(f"[segment] {plate['name']}: using existing masks in {seg['masks']}")
        return 0
    chans = seg["channels"]
    if not chans:
        raise ValueError("segmentation.channels is empty")
    sites = io.find_sites(plate, cfg)
    todo = [(k, v) for k, v in sites.items() if overwrite or not os.path.exists(io.mask_path(cfg, plate["name"], *k))]
    if not todo:
        return 0
    model = _model(seg)
    from importlib.metadata import version
    for (well, site), files in todo:
        ims = [io.read_image(files[m]) for m in chans]
        while len(ims) < 3:
            ims.append(np.zeros_like(ims[0]))
        masks = model.eval(np.stack(ims[:3], -1), channel_axis=-1, diameter=seg.get("diameter"),
                           flow_threshold=seg.get("flow_threshold", 0.4), cellprob_threshold=seg.get("cellprob_threshold", 0.0))[0]
        io.write_mask(io.mask_path(cfg, plate["name"], well, site), masks)
    info = dict(plate=plate["name"], model=seg.get("model"), cellpose=version("cellpose"), input_channels=chans,
                diameter=seg.get("diameter"), flow_threshold=seg.get("flow_threshold"), cellprob_threshold=seg.get("cellprob_threshold"),
                n_sites=len(sites))
    with open(os.path.join(cfg["output"], plate["name"], "masks", "segmentation_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    print(f"[segment] {plate['name']}: {len(todo)} sites segmented ({len(sites)} total)")
    return len(todo)
