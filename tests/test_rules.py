import numpy as np
import pandas as pd

from diffquant import measure, summarize, threshold as T


def test_gate_peak_minus_finds_main_population():
    rng = np.random.default_rng(0)
    lg = np.concatenate([rng.normal(2.6, 0.15, 5000), rng.normal(1.2, 0.2, 500)])   # cell type + K14-low other cells
    gated, low, info = T.gate_cells(lg, dict(rule="peak_minus", delta=0.8, upper_mad=4))
    assert abs(info["peak"] - 2.6) < 0.05
    assert abs(info["lo"] - 1.8) < 0.05
    assert 0.88 < gated.mean() < 0.92 and low[5000:].mean() > 0.99


def test_gate_fixed_and_no_upper():
    lg = T.log10p(np.array([5.0, 50.0, 500.0]))
    gated, low, _ = T.gate_cells(lg, dict(rule="fixed", value=40, upper_mad=None))
    assert gated.tolist() == [False, True, True] and low.tolist() == [True, False, False]


def test_readout_factor_scales_linear_cut():
    rng = np.random.default_rng(1)
    li = np.concatenate([rng.normal(2.0, 0.15, 4000), rng.normal(2.9, 0.15, 1000)])
    t1, to = T.readout_cut(li, dict(rule="otsu_times", factor=1.0))
    t13, _ = T.readout_cut(li, dict(rule="otsu_times", factor=1.3))
    assert np.isclose(t1, to)
    assert np.isclose(10 ** t13 - 1, 1.3 * (10 ** to - 1))
    tf, _ = T.readout_cut(li, dict(rule="fixed", value=99.0))
    assert np.isclose(tf, 2.0)


def test_background_rule_site_vs_plate():
    b = dict(min_cellfree_frac=0.2, min_sites=2)
    bgs = pd.DataFrame(dict(well="A01", site=[1, 2, 3], frac_cellfree=[0.5, 0.3, 0.05], X_bg_cellfree=[100.0, 110.0, 300.0]))
    out, plate, rule = measure.site_backgrounds(bgs, ["X"], b)
    assert out.X_background.tolist() == [100.0, 110.0, 105.0]          # dense site 3 gets the plate median
    assert out.bg_source.tolist() == ["site", "site", "plate_median"] and plate["X"] == 105.0


def test_bh():
    q = summarize.bh([0.01, 0.04, 0.03, np.nan])
    assert np.allclose(q[:3], [0.03, 0.04, 0.04]) and np.isnan(q[3])


def test_object_filter_excludes_small_and_border(example_cfg):
    import os
    from diffquant import measure, threshold
    cfg = example_cfg
    cfg["filters"] = dict(min_area_px=800, max_area_px=None, exclude_border=True)
    pl = cfg["plates"][0]
    os.makedirs(os.path.join(cfg["output"], pl["name"]), exist_ok=True)
    measure.run(cfg, pl)
    w = threshold.run(cfg, pl)
    c = pd.read_csv(os.path.join(cfg["output"], pl["name"], "cells_called.csv.gz"))
    bad = (c.area_px < 800) | c.touches_border.astype(bool)
    assert c.filtered.equals(bad) and not (c.gated & bad).any() and not (c.IVL_pos & bad).any()
    assert int(w.n_filtered.iloc[0]) == int(bad.sum()) > 0
