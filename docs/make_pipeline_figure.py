"""Pipeline overview figure (docs/img/pipeline.png) from the bundled example.
Run from the repository root after:  diffquant run examples/config.yaml
    python docs/make_pipeline_figure.py
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy import ndimage as ndi

from diffquant import config as C, io

cfg = C.load("examples/config.yaml")
pl = cfg["plates"][0]
od = os.path.join(cfg["output"], pl["name"])
files = io.find_sites(pl, cfg)[("B08", 1)]
mask_full = io.read_mask(io.find_masks(cfg, pl)[("B08", 1)])
cells = pd.read_csv(os.path.join(od, "cells_called.csv.gz"))
thr = pd.read_csv(os.path.join(od, "thresholds.csv")).iloc[0]
H = np.load(os.path.join(od, "histograms.npz"))

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 6,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "xtick.major.size": 2, "pdf.fonttype": 42,
                     "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
GREEN, MAG, BLUE, CYAN = np.array([0.1, 0.85, 0.1]), np.array([0.85, 0.1, 0.85]), np.array([0.2, 0.35, 1.0]), "#00b4d8"
ZY, ZX, ZS, PX = 620, 1060, 500, cfg["images"]["pixel_size_um"]
sl = (slice(ZY, ZY + ZS), slice(ZX, ZX + ZS))
raw = {k: io.read_image(files[k])[sl] for k in ("DAPI", "KRT14", "IVL")}
bg = {k: float(cells[f"{k}_background"].iloc[0]) for k in raw}
m = mask_full[sl]
far = (ndi.distance_transform_edt(mask_full == 0) > cfg["background"]["far_px"])[sl]


def nrm(a, lo, hi):
    return np.clip((a - lo) / (hi - lo), 0, 1)


def bnd(lab):
    b = np.zeros(lab.shape, bool)
    b[:-1] |= lab[:-1] != lab[1:]; b[1:] |= lab[1:] != lab[:-1]
    b[:, :-1] |= lab[:, :-1] != lab[:, 1:]; b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


rn = {k: nrm(raw[k], np.percentile(raw[k], 0.5), np.percentile(raw[k], 99.5)) for k in raw}
pr = {k: np.clip(raw[k] - bg[k], 0, None) for k in raw}
pn = {k: nrm(pr[k], 0, np.percentile(pr[k], 99.5)) for k in raw}

W, Hh = 183, 44
fig = plt.figure(figsize=(W / 25.4, Hh / 25.4))
T, X0, GAP, TOP = 22.5, 5.0, 8.2, 6.5
xs = [X0 + i * (T + GAP) for i in range(6)]
axs = [fig.add_axes([x / W, 1 - (TOP + T) / Hh, T / W, T / Hh]) for x in xs]
steps = [("1  Image", "nucleus · cell type · marker"),
         ("2  Segment", "Cellpose + optional\nsize / edge filter"),
         ("3  Background", "cell-free median"),
         ("4  Cell-type gate", f"drops other cells + debris\n(default: peak − {cfg['gate']['delta']})"),
         ("5  Marker cut", f"default: {cfg['readouts'][0]['factor']} × Otsu threshold"),
         ("6  % positive", "IVL+ of KRT14+, per well")]
axs[0].imshow(np.clip(rn["IVL"][..., None] * GREEN + rn["KRT14"][..., None] * MAG + rn["DAPI"][..., None] * BLUE, 0, 1))
L = 100 / PX
axs[0].add_patch(Rectangle((ZS - 25 - L, ZS - 40), L, 12, color="white", lw=0))
seg = np.clip(pn["KRT14"][..., None] * 0.75 + pn["DAPI"][..., None] * BLUE, 0, 1)
seg[bnd(m)] = [1.0, 0.85, 0.0]
axs[1].imshow(seg)
g3 = np.repeat(rn["IVL"][..., None], 3, -1) * 0.85
g3[far] = 0.55 * g3[far] + 0.45 * np.array(mpl.colors.to_rgb(CYAN))
axs[2].imshow(g3)
NL = int(max(mask_full.max(), cells.label.max())) + 1
call = np.zeros(NL, np.int8)
call[cells.label.values] = np.where(cells.gate_low.values | ~cells.gated.values, 3, np.where(cells.IVL_pos.values, 1, 2))
c6 = np.clip(pn["IVL"][..., None] * GREEN + pn["KRT14"][..., None] * MAG + pn["DAPI"][..., None] * BLUE, 0, 1)
bb, cl = bnd(m), call[m]
c6[bb & (cl == 1)] = [1, 1, 0.2]; c6[bb & (cl == 2)] = [0.6, 0.6, 0.6]; c6[bb & (cl == 3)] = [1, 0.3, 0.3]
axs[5].imshow(c6)
for a in (axs[0], axs[1], axs[2], axs[5]):
    a.set_xticks([]); a.set_yticks([])
    for s_ in a.spines.values():
        s_.set_visible(False)


def hist(a, v, cuts, col, xlim, xlab, span):
    a.hist(v, bins=np.linspace(*xlim, 60), color="0.65", lw=0)
    for t, ls in cuts:
        a.axvline(t, color=col, lw=0.8, ls=ls)
    a.axvspan(*span, color=col, alpha=0.08, lw=0)
    a.set_xlim(*xlim); a.set_yticks([]); a.spines["left"].set_visible(False); a.set_xticks([1, 2, 3])
    a.set_xlabel(xlab, labelpad=1)


pk, lo, hi = H["gate_cuts"]
hist(axs[3], H["gate"], [(lo, "-"), (hi, "-"), (pk, ":")], "#c51b7d", (0.5, 4.0), "log10 KRT14", (lo, hi))
t, to = H["IVL_cuts"]
hist(axs[4], H["IVL"], [(to, ":"), (t, "-")], "#1a9850", (1.0, 3.6), "log10 IVL", (t, 3.6))
w = pd.read_csv(os.path.join(od, "wells.csv")).iloc[0]
axs[5].text(0.04, 0.96, f"IVL+ {w.pct_IVL_pos:.1f} %\nof KRT14+", transform=axs[5].transAxes, ha="left", va="top", fontsize=5.5,
            color="white", bbox=dict(boxstyle="square,pad=0.15", fc="black", ec="none", alpha=0.6))
for i, (a, (t_, sub)) in enumerate(zip(axs, steps)):
    a.set_title(t_, loc="left", fontsize=6.5, fontweight="bold", pad=3)
    yoff = -0.36 if i in (3, 4) else -0.08
    a.text(0, yoff, sub, transform=a.transAxes, ha="left", va="top", fontsize=5.5, linespacing=1.2)
    if i < 5:
        y = 1 - (TOP + T / 2) / Hh
        fig.add_artist(FancyArrowPatch(((xs[i] + T + 1.0) / W, y), ((xs[i + 1] - 1.0) / W, y), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=6, lw=0.7, color="0.35"))
os.makedirs("docs/img", exist_ok=True)
fig.savefig("docs/img/pipeline.png", dpi=300)
fig.savefig("docs/img/pipeline.svg")
print("docs/img/pipeline.png")
