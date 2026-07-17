#!/usr/bin/env python3
"""
Fig. 4B - flanking-DNA conformational landscape from 3DVA.

Reproduces the latent-landscape panel: a grey scatter of all consensus
particles (the continuous distribution), plus a per-class density field for
each of the C-, S-, and I-shape classes with coloured centroid markers and
labels. Opacity is weighted by local class dominance, so each state's
distinctive region is emphasised while the heavily overlapping core is shown
faintly - representative geometries sampled from a single continuous landscape.

The output SVG is LAYERED: each element is written as a separately named SVG
group (layer_allparticles, layer_<C/S/I>_density, layer_<C/S/I>_centroid,
layer_<C/S/I>_label) so the colours can be toggled, recoloured, or reordered
independently in Illustrator / Inkscape.

INPUT (deposited on Zenodo, DOI: 10.5281/zenodo.20649656):
  3DVA_J2619_8A_latent_with_class.csv
  columns: uid, PC1, PC2, PC3, class   (class is 'C', 'S', 'I', or empty)
  3DVA run: flanking-DNA-ends mask, 8 A filter resolution, 3 components, C1.

The figure plots variability components PC2 (x) and PC3 (y) - the pair that
best resolves the (heavily overlapping) C/S/I classes in this run. Edit
XCOL/YCOL below to change.

USAGE:
  python make_fig4B.py                    # expects the CSV in the same folder
  python make_fig4B.py path/to/latent.csv # or pass the path (.csv or .csv.gz)

OUTPUT: Fig4B_latent_landscape.pdf / .png / .svg  (.svg is the layered version)

Dependencies: numpy, pandas, scipy, matplotlib  (see requirements.txt)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
import matplotlib.patheffects as pe

# ----------------------------------------------------------------- INPUT
CSV = sys.argv[1] if len(sys.argv) > 1 else "3DVA_J2619_8A_latent_with_class.csv"
OUT = "Fig4B_latent_landscape_bright"

# ----------------------------------------------------------------- STYLE (edit)
HUE = {"C": "#e8241c", "S": "#12b312", "I": "#1f7ff5"}   # brighter, more saturated C/S/I
XCOL, YCOL = "PC2", "PC3"        # components plotted (best C/S/I separation)
AXLIM = 25                       # axis range: -AXLIM .. +AXLIM on both axes
SCATTER_N = 120_000              # number of background (grey) particles drawn
SCATTER_ALPHA = 0.06             # transparency of the grey scatter (lowered so colours pop)
BINS = 52                        # 2D histogram bins per class (higher = sharper)
SMOOTH = 0.6                     # gaussian smoothing of each density (lower = sharper)
DENS_CUT = 0.08                  # floor: faint base only where density exceeds this
BASE_FLOOR = 0.28                # faint opacity in the shared/overlap core (raised)
RAMP_GAIN = 1.0                  # peak opacity where a class locally dominates (maxed)
RAMP_EXP = 0.62                  # dominance ramp exponent (<1 = contrast rises early)
DOT_SIZE = 140
LABEL_OFF = {"I": (-24, 4), "C": (2, 24), "S": (24, -6)}  # label offsets (points)
SEED = 0

# ----------------------------------------------------------------- LOAD
df = pd.read_csv(CSV)                    # pandas auto-detects .gz
df["class"] = df["class"].fillna("").astype(str)
zx, zy = df[XCOL].to_numpy(), df[YCOL].to_numpy()
lab = df["class"].to_numpy()
cen = {n: (zx[lab == n].mean(), zy[lab == n].mean()) for n in ["C", "S", "I"]}
for n in ["C", "S", "I"]:
    print(f"{n}: {(lab == n).sum()} particles")

# per-class normalized density on a common grid
xr = (-AXLIM, AXLIM); yr = (-AXLIM, AXLIM)
H = {}
for n in ["C", "S", "I"]:
    m = lab == n
    h, _, _ = np.histogram2d(zx[m], zy[m], bins=BINS, range=[xr, yr], density=True)
    h = ndimage.gaussian_filter(h.T, SMOOTH)
    H[n] = h / h.max()
tot = H["C"] + H["S"] + H["I"] + 1e-9     # local total (for dominance)

# ----------------------------------------------------------------- PLOT
# composite_image:False keeps each imshow a separate <image> so SVG layers survive
plt.rcParams.update({"font.size": 13, "font.family": "sans-serif",
                     "axes.linewidth": 0.9, "svg.fonttype": "none",
                     "image.composite_image": False})
fig, ax = plt.subplots(figsize=(6.2, 6.0))

# layer: all particles (grey continuous distribution)
idx = np.random.default_rng(SEED).choice(len(zx), min(SCATTER_N, len(zx)), replace=False)
sc = ax.scatter(zx[idx], zy[idx], s=1.2, c="0.82", alpha=SCATTER_ALPHA,
                rasterized=True, linewidths=0, zorder=1)
sc.set_gid("layer_allparticles")

# layer: per-class density, dominance-weighted opacity
for name in ["I", "C", "S"]:
    h = H[name]
    frac = h / tot                                    # local class share (1/3 = shared)
    g = np.clip((frac - 1/3) / (1 - 1/3), 0, 1)        # 0 in shared core -> 1 where exclusive
    ramp = g ** RAMP_EXP
    base = np.where(h > DENS_CUT, BASE_FLOOR, 0.0)
    alpha = np.clip(h * base + ramp * RAMP_GAIN, 0, 1)
    cmap = LinearSegmentedColormap.from_list(name, [(1, 1, 1), HUE[name]])
    im = ax.imshow(h, origin="lower", extent=[xr[0], xr[1], yr[0], yr[1]], aspect="auto",
                   cmap=cmap, interpolation="bilinear", alpha=alpha, vmin=0, zorder=2)
    im.set_gid(f"layer_{name}_density")

# layer: labels (centroid dots removed)
for name in ["I", "C", "S"]:
    tx = ax.annotate(name, cen[name], fontsize=16, fontweight="bold", style="italic",
                     color=HUE[name], xytext=LABEL_OFF[name], textcoords="offset points",
                     ha="center", va="center", zorder=8,
                     path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
    tx.set_gid(f"layer_{name}_label")

ax.set_xlim(xr); ax.set_ylim(yr)
ax.set_xticks([-20, -10, 0, 10, 20]); ax.set_yticks([-20, -10, 0, 10, 20])
ax.set_xlabel("variability component 2 (PC2)")
ax.set_ylabel("variability component 3 (PC3)")
fig.tight_layout()
fig.savefig(f"{OUT}.svg", bbox_inches="tight")          # layered
fig.savefig(f"{OUT}.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}.png", dpi=300, bbox_inches="tight")
print(f"wrote {OUT}.svg (layered) / .pdf / .png")
