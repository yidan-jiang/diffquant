"""End-to-end run on the bundled example field (precomputed cpsam masks, no Cellpose needed)."""
import os

import numpy as np
import pandas as pd

from conftest import ROOT
from diffquant import calibrate, focus, measure, overlay, scoring, summarize, threshold

EXP = os.path.join(ROOT, "tests", "expected")


def test_example_end_to_end(example_cfg):
    cfg = example_cfg
    pl = cfg["plates"][0]
    os.makedirs(os.path.join(cfg["output"], pl["name"]), exist_ok=True)
    cells = measure.run(cfg, pl)
    assert len(cells) == 521
    w = threshold.run(cfg, pl)
    exp = pd.read_csv(os.path.join(EXP, "example_wells.csv"))
    pd.testing.assert_frame_equal(w.reset_index(drop=True), exp, check_dtype=False)
    thr = pd.read_csv(os.path.join(cfg["output"], pl["name"], "thresholds.csv"))
    exp_t = pd.read_csv(os.path.join(EXP, "example_thresholds.csv"))
    pd.testing.assert_frame_equal(thr, exp_t, check_dtype=False)
    called = pd.read_csv(os.path.join(cfg["output"], pl["name"], "cells_called.csv.gz"))
    exp_c = pd.read_csv(os.path.join(EXP, "example_cells_called.csv.gz"))
    pd.testing.assert_frame_equal(called, exp_c, check_dtype=False)
    focus.run(cfg, pl)
    summarize.run(cfg)
    assert overlay.run(cfg, pl, sites=[("B08", 1)])


def test_scoring_and_calibration_recover_rule(example_cfg):
    """Fake 'manual' calls made with a known rule (gate delta 0.6, readout 1.5 x Otsu) are recovered by calibrate."""
    cfg = example_cfg
    pl = cfg["plates"][0]
    os.makedirs(os.path.join(cfg["output"], pl["name"]), exist_ok=True)
    measure.run(cfg, pl)
    threshold.run(cfg, pl)
    scoring.run(cfg, pl, per_decile=6, n_random=10, per_band=12)
    thr = pd.read_csv(os.path.join(cfg["output"], pl["name"], "thresholds.csv")).iloc[0]
    for mk in ("IVL", "KRT14"):
        sd = os.path.join(cfg["output"], "scoring", f"example_{mk}")
        ix = pd.read_csv(os.path.join(sd, "index.csv"))
        assert os.path.exists(os.path.join(sd, "crops.tif")) and len(ix) > (40 if mk == "IVL" else 10)   # one field: few cells near the gate
        if mk == "IVL":
            pos = ix.IVL_mean > 1.5 * thr.IVL_otsu
        else:
            pos = np.log10(ix.KRT14_mean + 1) - np.log10(thr.gate_peak + 1) > -0.6
        pd.DataFrame(dict(slice=ix.slice, call=np.where(pos, "positive", "negative"), scorer="TEST", seconds=1.0,
                          timestamp="x")).to_csv(os.path.join(sd, "calls_TEST.csv"), index=False)
    s = calibrate.run(cfg, scorer="TEST")
    assert abs(s["IVL"]["recommended_factor"] - 1.5) <= 0.1
    assert s["IVL"]["MAE_at_recommended"] < 3
    assert s["gate"]["plateau"][0] <= 0.6 <= s["gate"]["plateau"][1] + 0.05
    assert s["consistency"]["IVL"]["identical"] == 1.0
