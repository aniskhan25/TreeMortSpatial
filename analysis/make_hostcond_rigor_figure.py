"""Figure: tree-level vs centroid host-conditioned L_inhom(r)-r per block, out to 3 km.
Reads the spatstat envelope curves dumped to analysis/out/linhom_{centroids,trees}.csv.
Each panel: observed L_inhom(r)-r against the host-conditioned simulation envelope (shaded);
the 0 line is the host-conditioned expectation. Curves above the band = residual aggregation
beyond what forest host explains. Both units, both >1 km, confirm the same finding.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, "analysis")
from figstyle import apply_style, PALETTE, despine, panel_label, save

apply_style()
cen = pd.read_csv("analysis/out/linhom_centroids.csv")
tre = pd.read_csv("analysis/out/linhom_trees.csv")
blocks = sorted(cen.block.unique())

fig, axes = plt.subplots(2, 4, figsize=(12.5, 6.2), sharex=True)
axes = axes.ravel()
for ax, b in zip(axes, blocks):
    c = cen[cen.block == b]
    t = tre[tre.block == b]
    # host-conditioned envelopes (shaded) + observed curves, in metres -> km
    ax.fill_between(c.r / 1000, c.lo, c.hi, color=PALETTE["muted"], alpha=0.35, lw=0, zorder=1)
    ax.fill_between(t.r / 1000, t.lo, t.hi, color=PALETTE["secondary"], alpha=0.18, lw=0, zorder=1)
    ax.plot(c.r / 1000, c.obs, color=PALETTE["primary"], lw=1.6, zorder=3, label="clusters")
    ax.plot(t.r / 1000, t.obs, color=PALETTE["accent"], lw=1.6, zorder=3, label="trees")
    ax.axhline(0, color=PALETTE["neutral"], lw=0.8, ls="--", zorder=2)
    despine(ax)
    ax.set_title(f"Block {b}", fontsize=10)
    ax.set_xlim(0, 3)
    ax.set_ylim(bottom=min(-50, c.lo.min() * 1.1))

# legend / key in the unused 8th panel
key = axes[len(blocks)]
key.axis("off")
key.plot([], [], color=PALETTE["primary"], lw=1.6, label="Observed (clusters, $\\leq$3 km)")
key.plot([], [], color=PALETTE["accent"], lw=1.6, label="Observed (trees, $\\leq$2 km)")
key.fill_between([], [], [], color=PALETTE["muted"], alpha=0.35, label="Host-conditioned envelope (clusters)")
key.fill_between([], [], [], color=PALETTE["secondary"], alpha=0.18, label="Host-conditioned envelope (trees)")
key.axhline(0, color=PALETTE["neutral"], lw=0.8, ls="--", label="Host expectation (0)")
key.legend(loc="center", fontsize=9, frameon=False, handlelength=1.8)

for ax in list(axes[4:]) + [axes[3]]:
    if ax is not key:
        ax.set_xlabel("distance $r$ (km)")
key.set_xlabel("")
for r in (0, 4):
    axes[r].set_ylabel(r"$L_{\mathrm{inhom}}(r)-r$ (m)")

fig.suptitle("Residual aggregation above the host-conditioned null persists at the tree level and beyond 1 km",
             fontsize=11.5, y=1.00)
for letter, ax in zip("abcdefg", axes[:len(blocks)]):
    panel_label(ax, letter)
fig.tight_layout(rect=[0, 0, 1, 0.97])
save(fig, "figures/fig_hostcond_rigor")
print("wrote figures/fig_hostcond_rigor.png/.pdf")
