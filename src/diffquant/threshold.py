"""Per-plate, data-driven cell-type gate and readout cut-offs, then % positive per well.

All rules work on log10(mean + 1) of the background-subtracted cell means and are computed for every plate separately
(so plate-to-plate differences in staining or exposure are absorbed by the rule, not by a fixed number).

Object filter (optional, before gating): min / max area, cells touching the image border. Filtered objects are never gated
and do not enter the per-plate rules.
Gate (cell type, e.g. KRT14+ keratinocytes):
    peak_minus  cut = mode of the smoothed log10 histogram - delta   (the cell type of interest forms the main peak)
    otsu        cut = Otsu threshold of log10 values
    fixed       cut = log10(value + 1)
    upper_mad   cells brighter than median + k * MAD (of cells above the cut) are dropped as outliers
Readout (e.g. involucrin), among gated cells:
    otsu_times  cut = factor x Otsu threshold (Otsu on log10 values, factor applied on the linear intensity scale)
    fixed       cut = value
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from skimage.filters import threshold_otsu

from .config import HIST, gate_marker


def log10p(x):
    return np.log10(np.clip(np.asarray(x, float), 0, None) + 1.0)


def smooth_hist(v, lo=HIST["lo"], hi=HIST["hi"], nb=HIST["bins"], sig=HIST["sigma"]):
    h, e = np.histogram(v, bins=nb, range=(lo, hi))
    return gaussian_filter1d(h.astype(float), sig), (e[:-1] + e[1:]) / 2


def peak(v) -> float:
    h, c = smooth_hist(v)
    return float(c[int(np.argmax(h))])


def gate_cells(lg, g: dict):
    """Returns (gated, low, info) for log10 gate-marker values."""
    pk = peak(lg)
    if g["rule"] == "peak_minus":
        lo = pk - float(g["delta"])
    elif g["rule"] == "otsu":
        lo = float(threshold_otsu(lg, nbins=256))
    else:
        lo = float(log10p(g["value"]))
    pos0 = lg > lo
    hi = np.inf
    if g.get("upper_mad") is not None and pos0.any():
        med = np.median(lg[pos0])
        mad = 1.4826 * np.median(np.abs(lg[pos0] - med))
        hi = float(med + float(g["upper_mad"]) * mad)
    return pos0 & (lg <= hi), lg <= lo, dict(peak=pk, lo=lo, hi=hi)


def readout_cut(li_gated, r: dict):
    """Returns (cut_log10, otsu_log10 or nan)."""
    if r.get("rule", "otsu_times") == "fixed":
        return float(log10p(r["value"])), np.nan
    t = float(threshold_otsu(li_gated, nbins=256))
    return float(np.log10(float(r.get("factor", 1.0)) * (10 ** t - 1) + 1)), t


def _expand(sub, keep):
    out = np.zeros(len(keep), bool)
    out[keep] = sub
    return out


def lin(t):
    return round(float(10 ** t - 1), 1) if np.isfinite(t) else np.nan


def run(cfg: dict, plate: dict) -> pd.DataFrame:
    od = os.path.join(cfg["output"], plate["name"])
    cells = pd.read_csv(os.path.join(od, "cells.csv.gz"))
    gm = gate_marker(cfg)
    f = cfg.get("filters") or {}
    keep = np.ones(len(cells), bool)
    if f.get("min_area_px") is not None:
        keep &= cells.area_px.values >= float(f["min_area_px"])
    if f.get("max_area_px") is not None:
        keep &= cells.area_px.values <= float(f["max_area_px"])
    if f.get("exclude_border"):
        keep &= ~cells.touches_border.astype(bool).values
    cells["filtered"] = ~keep
    row = dict(plate=plate["name"], n_cells=len(cells), n_filtered=int((~keep).sum()))
    hist = {}
    if gm:
        lg = log10p(cells[f"{gm}_mean"].values)
        gated, low, gi = gate_cells(lg[keep], cfg["gate"])
        gated = _expand(gated, keep)
        low = _expand(low, keep)
        row.update(gate_marker=gm, gate_rule=cfg["gate"]["rule"], gate_delta=cfg["gate"].get("delta"),
                   gate_peak=lin(gi["peak"]), gate_cut=lin(gi["lo"]), gate_upper=lin(gi["hi"]),
                   frac_gated=round(float(gated.mean()), 4), frac_gate_low=round(float(low.mean()), 4))
        hist["gate"] = lg[keep].astype(np.float32)
        hist["gate_cuts"] = np.array([gi["peak"], gi["lo"], gi["hi"]])
    else:
        gated = keep.copy()
        low = np.zeros(len(cells), bool)
        row.update(gate_marker="none", frac_gated=1.0)
    cells["gated"] = gated
    cells["gate_low"] = low
    for r in cfg["readouts"]:
        mk = r["marker"]
        li = log10p(cells[f"{mk}_mean"].values)
        t, to = readout_cut(li[gated], r)
        pos = gated & (li > t)
        cells[f"{mk}_pos"] = pos
        fpos = float(pos.sum() / max(gated.sum(), 1))
        check = r.get("rule", "otsu_times") == "otsu_times" and not (0.02 <= fpos <= 0.95)
        row.update({f"{mk}_rule": r.get("rule", "otsu_times"), f"{mk}_factor": r.get("factor"), f"{mk}_otsu": lin(to),
                    f"{mk}_cut": lin(t), f"{mk}_frac_pos": round(fpos, 4), f"{mk}_check_histogram": bool(check)})
        if check:
            print(f"[threshold] {plate['name']}: {mk}+ {100 * fpos:.1f} % - check thresholds.png / overlays: with (almost) one "
                  f"population the Otsu threshold may split it; if unimodal, use a fixed cut from bimodal plates")
        hist[mk] = li[gated].astype(np.float32)
        hist[f"{mk}_cuts"] = np.array([t, to])
    g = cells.groupby("well")
    w = pd.DataFrame(dict(n_cells=g.size(), n_sites=g.site.nunique(), n_filtered=g.filtered.sum(), n_gated=g.gated.sum(),
                          n_gate_low=g.gate_low.sum())).reset_index()
    w["frac_gate_low"] = (w.n_gate_low / w.n_cells.clip(lower=1)).round(4)
    for r in cfg["readouts"]:
        mk = r["marker"]
        w[f"n_{mk}_pos"] = g[f"{mk}_pos"].sum().values
        w[f"pct_{mk}_pos"] = (100 * w[f"n_{mk}_pos"] / w.n_gated.clip(lower=1)).round(2)
    w["valid"] = w.n_gated >= int(cfg["min_gated_cells"])
    w.insert(0, "plate", plate["name"])
    cells.to_csv(os.path.join(od, "cells_called.csv.gz"), index=False, compression="gzip")
    w.to_csv(os.path.join(od, "wells.csv"), index=False)
    pd.DataFrame([row]).to_csv(os.path.join(od, "thresholds.csv"), index=False)
    np.savez_compressed(os.path.join(od, "histograms.npz"), **hist)
    _plot(cfg, plate["name"], hist, os.path.join(od, "thresholds.png"))
    msg = ", ".join(f"{r['marker']}+ {100 * row[r['marker'] + '_frac_pos']:.1f} %" for r in cfg["readouts"])
    print(f"[threshold] {plate['name']}: {row['frac_gated'] * 100:.1f} % gated ({gm or 'no gate'}), {msg}; "
          f"{int(w.valid.sum())}/{len(w)} wells valid")
    return w


def _plot(cfg, name, hist, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    gm = gate_marker(cfg)
    panels = ([("gate", gm, "#c51b7d")] if gm else []) + [(r["marker"], r["marker"], "#1a9850") for r in cfg["readouts"]]
    fig, axs = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 2.4), squeeze=False)
    for ax, (key, label, col) in zip(axs[0], panels):
        ax.hist(hist[key], bins=150, color="0.6", lw=0)
        for i, t in enumerate(hist[f"{key}_cuts"]):
            if np.isfinite(t):
                ax.axvline(t, color=col, lw=0.8, ls=":" if (key == "gate" and i == 0) or (key != "gate" and i == 1) else "-")
        ax.set_xlabel(f"log10 {label} (cell mean)" + ("" if key == "gate" else f", {gm}+ cells" if gm else ""))
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
    axs[0][0].set_title(f"{name}  (solid = cut, dotted = {'peak / ' if gm else ''}Otsu)", loc="left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
