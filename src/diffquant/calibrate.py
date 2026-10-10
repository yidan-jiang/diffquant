"""Calibrate the gate delta and the readout Otsu factor against blinded manual calls.

Input: scoring sets (``diffquant score-sets``) scored in Fiji (``fiji/score_cells.py`` -> calls_<scorer>.csv), plus the current
per-cell values and plate thresholds. Only first presentations (repeat = 0) are used; repeats give scorer consistency.

Gate:    rel = log10(value + 1) - log10(plate peak + 1); a cell is gate+ if rel > -delta. Agreement with the calls
         (unsure dropped) for delta = 0.40 ... 1.20; the recommended delta is the middle of the plateau of best agreement.
Readout: per plate, manual % positive = mean over intensity-decile strata of the mean call score (positive 1, unsure 0.5,
         negative 0) - this undoes the stratified sampling; the rule % is computed on the same cells, weighted the same way,
         for cut = factor x Otsu, factor = 0.80 ... 2.00. Bias and mean absolute error (MAE) across plates;
         the recommended factor minimises the MAE (bootstrap over plates gives its stability).
"""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

SCORE = {"positive": 1.0, "negative": 0.0, "unsure": 0.5}
DELTAS = np.round(np.arange(0.40, 1.201, 0.05), 2)
FACTORS = np.round(np.arange(0.80, 2.001, 0.05), 2)


def collect(cfg: dict, scorer: str | None = None) -> pd.DataFrame:
    """Tidy table: one row per scored crop with the current value of the scored marker and the plate thresholds."""
    rows = []
    for sd in sorted(glob.glob(os.path.join(cfg["output"], "scoring", "*"))):
        meta = {}
        if os.path.exists(os.path.join(sd, "display.txt")):
            for line in open(os.path.join(sd, "display.txt")):
                if line.strip() and not line.startswith("#"):
                    k, *v = line.split()
                    meta[k] = v[0] if v else ""
        calls = sorted(glob.glob(os.path.join(sd, "calls_*.csv")))
        if scorer:
            calls = [c for c in calls if os.path.basename(c) == f"calls_{scorer}.csv"]
        if not calls or "marker" not in meta:
            continue
        plate = os.path.basename(sd)[: -(len(meta["marker"]) + 1)]
        ix = pd.read_csv(os.path.join(sd, "index.csv"))
        cells = pd.read_csv(os.path.join(cfg["output"], plate, "cells_called.csv.gz"))
        thr = pd.read_csv(os.path.join(cfg["output"], plate, "thresholds.csv")).iloc[0]
        mk = meta["marker"]
        for cf in calls:
            cl = pd.read_csv(cf)
            d = ix[["slice", "well", "site", "label", "repeat", "stratum"]].merge(cl[["slice", "call", "scorer"]], on="slice")
            d = d.merge(cells[["well", "site", "label", f"{mk}_mean"]].rename(columns={f"{mk}_mean": "value"}), on=["well", "site", "label"])
            d["plate"], d["marker"], d["role"] = plate, mk, meta.get("role", "readout")
            d["peak"] = thr.get("gate_peak", np.nan)
            d["otsu"] = thr.get(f"{mk}_otsu", np.nan)
            rows.append(d)
    if not rows:
        raise SystemExit("[calibrate] no scored sets found (output/scoring/*/calls_<scorer>.csv)")
    return pd.concat(rows, ignore_index=True)


def youden(x, y):
    best = (-1.0, np.nan)
    for t in np.unique(x):
        p = x > t
        j = (p & y).sum() / max(y.sum(), 1) + (~p & ~y).sum() / max((~y).sum(), 1) - 1
        if j > best[0]:
            best = (float(j), float(t))
    return best


def gate_curve(calls: pd.DataFrame, deltas=DELTAS):
    k = calls[(calls.role == "gate") & (calls.repeat == 0) & (calls.call != "unsure")].copy()
    if k.empty:
        return None, None, {}
    k["rel"] = np.log10(k.value + 1) - np.log10(k.peak + 1)
    y = (k.call == "positive").values
    curve = pd.DataFrame({"delta": deltas, "agreement": [((k.rel > -d) == y).mean() for d in deltas]})
    curve["n_cells"], curve["n_plates"] = len(k), k.plate.nunique()
    best = curve.agreement.max()
    on = curve.agreement >= best - 0.006                 # plateau: within 0.6 % of the best agreement
    ib = int(curve.agreement.values.argmax())
    lo = ib
    while lo > 0 and on.iloc[lo - 1]:
        lo -= 1
    hi = ib
    while hi < len(curve) - 1 and on.iloc[hi + 1]:
        hi += 1
    rec = round(round((curve.delta.iloc[lo] + curve.delta.iloc[hi]) / 2 / 0.05) * 0.05, 2)
    per = []
    for pid, g in k.groupby("plate"):
        j, t = youden(g.rel.values, (g.call == "positive").values)
        per.append(dict(plate=pid, n=len(g), n_pos=int((g.call == "positive").sum()), youden_delta=round(-t, 2),
                        agreement_at_recommended=round(float(((g.rel > -rec) == (g.call == "positive")).mean()), 3)))
    j, t = youden(k.rel.values, y)
    info = dict(n_cells=int(len(k)), n_plates=int(k.plate.nunique()), best_agreement=round(float(best), 3),
                plateau=[float(curve.delta.iloc[lo]), float(curve.delta.iloc[hi])], recommended_delta=float(rec),
                agreement_at_recommended=round(float(curve.set_index("delta").agreement.get(rec, np.nan)), 3),
                youden_delta=round(-t, 3), youden_J=round(j, 3))
    return curve.round(4), pd.DataFrame(per), info


def readout_curve(calls: pd.DataFrame, factors=FACTORS, n_boot: int = 2000):
    out = {}
    for mk, iv in calls[(calls.role == "readout") & (calls.repeat == 0)].groupby("marker"):
        iv = iv.assign(s=iv.call.map(SCORE))
        st = iv[iv.stratum >= 0]
        man = st.groupby(["plate", "stratum"]).s.mean().groupby("plate").mean() * 100
        R = pd.DataFrame({f: st.assign(p=st.value > f * st.otsu).groupby(["plate", "stratum"]).p.mean().groupby("plate").mean() * 100
                          for f in factors})
        err = R.sub(man, axis=0)
        d = iv[iv.call != "unsure"]
        y = (d.call == "positive").values
        S = pd.DataFrame({"factor": factors, "bias": err.mean().values, "MAE": err.abs().mean().values,
                          "max_abs": err.abs().max().values,
                          "cell_agreement": [((d.value > f * d.otsu) == y).mean() for f in factors]}).round(3)
        rng = np.random.default_rng(1)
        bm = [err.iloc[rng.integers(0, len(err), len(err))].abs().mean().idxmin() for _ in range(n_boot)]
        rec = float(err.abs().mean().idxmin())
        plates = pd.DataFrame({"manual_pct": man.round(1), "auto_pct_at_recommended": R[rec].round(1), "auto_pct_at_1.0": R[1.0].round(1)
                               if 1.0 in R else np.nan, "n_scored": iv.groupby("plate").size()}).reset_index()
        info = dict(n_plates=int(len(man)), n_cells=int(len(d)), recommended_factor=rec,
                    zero_bias_factor=float(err.mean().abs().idxmin()),
                    bootstrap_median=float(np.median(bm)), bootstrap_95=[float(x) for x in np.percentile(bm, [2.5, 97.5])],
                    bias_at_recommended=round(float(err[rec].mean()), 1), MAE_at_recommended=round(float(err[rec].abs().mean()), 1),
                    cell_agreement_at_recommended=round(float(((d.value > rec * d.otsu) == y).mean()), 3))
        out[mk] = (S, plates, info)
    return out


def consistency(calls: pd.DataFrame) -> dict:
    """Agreement between first and repeated presentation of the same cell."""
    rep = {}
    for mk, g in calls.groupby("marker"):
        a = g[g.repeat == 0].set_index(["plate", "well", "site", "label", "scorer"]).call
        b = g[g.repeat == 1].set_index(["plate", "well", "site", "label", "scorer"]).call
        j = a.to_frame("first").join(b.to_frame("second"), how="inner")
        rep[mk] = dict(n_repeats=int(len(j)), identical=round(float((j["first"] == j["second"]).mean()), 3) if len(j) else None)
    return rep


def run(cfg: dict, scorer: str | None = None, exclude=(), calls: pd.DataFrame | None = None, out: str | None = None) -> dict:
    calls = collect(cfg, scorer) if calls is None else calls
    calls = calls[~calls.plate.astype(str).isin([str(e) for e in exclude])]
    out = out or os.path.join(cfg["output"], "calibration")
    os.makedirs(out, exist_ok=True)
    calls.to_csv(os.path.join(out, "calls.csv"), index=False)
    summary = {"consistency": consistency(calls)}
    curve, per, gi = gate_curve(calls)
    if curve is not None:
        curve.to_csv(os.path.join(out, "gate_agreement_curve.csv"), index=False)
        per.to_csv(os.path.join(out, "gate_per_plate.csv"), index=False)
        summary["gate"] = gi
    ro = readout_curve(calls)
    for mk, (S, plates, info) in ro.items():
        S.to_csv(os.path.join(out, f"{mk}_factor_curve.csv"), index=False)
        plates.to_csv(os.path.join(out, f"{mk}_per_plate.csv"), index=False)
        summary[mk] = info
    with open(os.path.join(out, "calibration_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    plot(curve, ro, gi, os.path.join(out, "calibration.png"))
    print(f"[calibrate] -> {out}")
    if gi:
        print(f"  gate: recommended delta = {gi['recommended_delta']} (plateau {gi['plateau'][0]}-{gi['plateau'][1]}, "
              f"agreement {gi['agreement_at_recommended']}, {gi['n_cells']} cells / {gi['n_plates']} plates)")
    for mk, (_, _, info) in ro.items():
        print(f"  {mk}: recommended factor = {info['recommended_factor']} x Otsu (bias {info['bias_at_recommended']:+.1f}, "
              f"MAE {info['MAE_at_recommended']:.1f} % points, {info['n_plates']} plates; bootstrap 95 % {info['bootstrap_95']})")
    return summary


def plot(curve, ro, gi, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = (curve is not None) + 2 * len(ro)
    fig, axs = plt.subplots(1, max(n, 1), figsize=(2.6 * max(n, 1), 2.4), squeeze=False)
    axs = list(axs[0])
    for mk, (S, plates, info) in ro.items():
        a = axs.pop(0)
        a.plot(S.factor, S.bias, color="0.3", lw=1, label="bias")
        a.plot(S.factor, S.MAE, color="#1a9850", lw=1.2, label="MAE")
        a.axhline(0, color="0.6", lw=0.5, ls=":")
        a.axvline(info["recommended_factor"], color="#1a9850", lw=0.8, ls="--")
        a.set_xlabel(f"{mk} cut (x Otsu)")
        a.set_ylabel("auto - manual (% points)")
        a.legend(frameon=False, fontsize=7)
        a = axs.pop(0)
        a.plot([0, 100], [0, 100], color="0.6", lw=0.5, ls=":")
        a.scatter(plates.manual_pct, plates.auto_pct_at_recommended, s=10, color="0.3")
        m = max(plates.manual_pct.max(), plates.auto_pct_at_recommended.max()) * 1.1
        a.set_xlim(0, m)
        a.set_ylim(0, m)
        a.set_xlabel(f"manual {mk}+ % (blinded)")
        a.set_ylabel(f"automatic {mk}+ % ({info['recommended_factor']} x Otsu)")
    if curve is not None:
        a = axs.pop(0)
        a.plot(curve.delta, curve.agreement, color="#c51b7d", lw=1.2, marker="o", ms=2)
        a.axvspan(*gi["plateau"], color="#c51b7d", alpha=0.1, lw=0)
        a.axvline(gi["recommended_delta"], color="#c51b7d", lw=0.8, ls="--")
        a.set_xlabel("gate delta (log10 below peak)")
        a.set_ylabel("agreement with manual calls")
    for a in fig.axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
