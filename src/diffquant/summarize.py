"""Collect wells of all plates; with a plate map, summarise conditions and the change vs the matched vehicle.

Plate map (CSV): plate, well, condition, [vehicle]
    vehicle = name of the condition on the same plate that this condition is compared with (blank for the vehicle itself).
Per plate and condition: mean % positive of valid wells (>= min_gated_cells). Delta = condition - its vehicle, per plate.
Across plates (plate = replicate): mean Delta, SD, 95 % CI, two-sided one-sample t-test, Benjamini-Hochberg q.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy import stats


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        q_ok = np.empty(len(pv))
        q_ok[o] = np.minimum(r, 1)
        q[ok] = q_ok
    return q


def run(cfg: dict) -> None:
    out = cfg["output"]
    wells, thr = [], []
    for p in cfg["plates"]:
        d = os.path.join(out, p["name"])
        if os.path.exists(os.path.join(d, "wells.csv")):
            wells.append(pd.read_csv(os.path.join(d, "wells.csv"), dtype={"plate": str}))
            thr.append(pd.read_csv(os.path.join(d, "thresholds.csv"), dtype={"plate": str}))
    if not wells:
        raise SystemExit("[summarize] no wells.csv found - run 'diffquant threshold' first")
    W = pd.concat(wells, ignore_index=True)
    pd.concat(thr, ignore_index=True).to_csv(os.path.join(out, "thresholds_all.csv"), index=False)
    fq = [os.path.join(out, p["name"], "focus.csv") for p in cfg["plates"]]
    fq = [pd.read_csv(f, dtype={"plate": str}) for f in fq if os.path.exists(f)]
    if fq:
        F = pd.concat(fq).groupby(["plate", "well"]).focus_flag.sum().rename("n_sites_focus_flag").reset_index()
        W = W.merge(F, on=["plate", "well"], how="left")
    if not cfg.get("plate_map"):
        W.to_csv(os.path.join(out, "wells_all.csv"), index=False)
        print(f"[summarize] {len(W)} wells -> wells_all.csv (no plate_map: no condition summary)")
        return
    pm = pd.read_csv(cfg["plate_map"], dtype={"plate": str})
    W = W.merge(pm, on=["plate", "well"], how="left")
    W.to_csv(os.path.join(out, "wells_all.csv"), index=False)
    V = W[W.valid & W.condition.notna()]
    rows, deltas = [], []
    for r in cfg["readouts"]:
        col = f"pct_{r['marker']}_pos"
        C = V.groupby(["plate", "condition"])[col].agg(["mean", "size"]).rename(columns={"mean": "pct", "size": "n_wells"}).reset_index()
        C["marker"] = r["marker"]
        rows.append(C)
        if "vehicle" in pm:
            veh = pm.dropna(subset=["vehicle"]).drop_duplicates(["plate", "condition"])[["plate", "condition", "vehicle"]]
            D = veh.merge(C, on=["plate", "condition"]).merge(
                C[["plate", "condition", "pct"]].rename(columns={"condition": "vehicle", "pct": "pct_vehicle"}), on=["plate", "vehicle"])
            D["delta"] = D.pct - D.pct_vehicle
            deltas.append(D)
    pd.concat(rows).to_csv(os.path.join(out, "conditions_per_plate.csv"), index=False)
    if deltas:
        D = pd.concat(deltas, ignore_index=True)
        D.to_csv(os.path.join(out, "delta_per_plate.csv"), index=False)
        S = []
        for (mk, cond, veh), g in D.groupby(["marker", "condition", "vehicle"]):
            v = g.delta.values
            n = len(v)
            se = v.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
            hw = stats.t.ppf(0.975, n - 1) * se if n > 1 else np.nan
            p = stats.ttest_1samp(v, 0).pvalue if n > 1 else np.nan
            S.append(dict(marker=mk, condition=cond, vehicle=veh, n_plates=n, mean_delta=v.mean(), sd=v.std(ddof=1) if n > 1 else np.nan,
                          ci95_lo=v.mean() - hw, ci95_hi=v.mean() + hw, p=p))
        S = pd.DataFrame(S)
        S["q_BH"] = np.nan
        for mk, g in S.groupby("marker"):
            S.loc[g.index, "q_BH"] = bh(g.p.values)
        S.round(4).to_csv(os.path.join(out, "delta_summary.csv"), index=False)
        print(f"[summarize] {len(W)} wells, {len(S)} condition comparisons -> delta_summary.csv")
    else:
        print(f"[summarize] {len(W)} wells -> conditions_per_plate.csv (plate map has no 'vehicle' column: no deltas)")
