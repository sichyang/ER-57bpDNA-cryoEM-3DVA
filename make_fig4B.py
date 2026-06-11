#!/usr/bin/env python3
"""
Fig. 4B - flanking-DNA conformational landscape from 3DVA.

Reproduces the latent-landscape panel: black scatter of all consensus
particles, per-class translucent density fields for the C/S/I shapes, and
italic letters at the class centroids.

INPUT (deposited on Zenodo, DOI: <INSERT>):
  3DVA_J2538_8A_latent_with_class.csv.gz
  columns: uid, PC1, PC2, PC3, class   (class is 'C', 'S', 'I', or empty)

USAGE:
  python make_fig4B.py                       # expects the CSV in the same folder
  python make_fig4B.py path/to/latent.csv.gz # or pass the path

OUTPUT: Fig4B_latent_landscape.pdf / .png / .svg

Dependencies: numpy, pandas, matplotlib  (see requirements.txt)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe

# ----------------------------------------------------------------- INPUT
CSV = sys.argv[1] if len(sys.argv) > 1 else "3DVA_J2538_8A_latent_with_class.csv.gz"
OUT = "Fig4B_latent_landscape"

# ----------------------------------------------------------------- STYLE (edit)
HUE = {"C": "#e8998d", "S": "#8fce8f", "I": "#5b9bd5"}  # C/S/I colours (match Fig. 1)
AXLIM = 35              # axis range: -AXLIM .. +AXLIM on both axes
SCATTER_N = 150_000    # number of background particles drawn (subsample)
SCATTER_ALPHA = 0.06   # transparency of the black scatter
DENS_BINS = 28         # 2D histogram bins for the per-class density fields
DENS_CUT = 0.16        # drop class density below this fraction of its max
SEED = 0

# ----------------------------------------------------------------- LOAD
df = pd.read_csv(CSV)                    # pandas auto-detects .gz
df["class"] = df["class"].fillna("").astype(str)
zx, zy = df["PC1"].to_numpy(), df["PC2"].to_numpy()
lab = df["class"].to_numpy()
cen = {n: (zx[lab == n].mean(), zy[lab == n].mean()) for n in ["C", "S", "I"]}
for n in ["C", "S", "I"]:
    print(f"{n}: {(lab == n).sum()} particles")

# ----------------------------------------------------------------- PLOT
plt.rcParams.update({"font.size": 13, "font.family": "sans-serif",
                     "axes.linewidth": 0.9, "svg.fonttype": "none"})
xr = (-AXLIM, AXLIM); yr = (-AXLIM, AXLIM)
fig, ax = plt.subplots(figsize=(6.2, 6.0))

idx = np.random.default_rng(SEED).choice(len(zx), min(SCATTER_N, len(zx)), replace=False)
ax.scatter(zx[idx], zy[idx], s=1.3, c="k", alpha=SCATTER_ALPHA,
           rasterized=True, linewidths=0, zorder=1)

for name in ["I", "S", "C"]:
    m = lab == name
    H, _, _ = np.histogram2d(zx[m], zy[m], bins=DENS_BINS, range=[xr, yr], density=True)
    H = H.T
    a = H / H.max(); a[a < DENS_CUT] = 0
    cmap = LinearSegmentedColormap.from_list(name, [(1, 1, 1), HUE[name]])
    ax.imshow(H, origin="lower", extent=[xr[0], xr[1], yr[0], yr[1]], aspect="auto",
              cmap=cmap, interpolation="gaussian", alpha=a * 0.92, vmin=0, zorder=2)

for name in ["C", "S", "I"]:
    ax.text(*cen[name], name, fontsize=22, fontweight="bold", style="italic",
            ha="center", va="center", color="k", zorder=7,
            path_effects=[pe.withStroke(linewidth=3.4, foreground="white")])

ax.set_xlim(xr); ax.set_ylim(yr)
ax.set_xlabel("variability component 1 (PC1)")
ax.set_ylabel("variability component 2 (PC2)")
fig.tight_layout()
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"{OUT}.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {OUT}.pdf/.png/.svg")
