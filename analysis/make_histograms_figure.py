"""
Publication-grade regeneration of fig2_histograms.png.

Design choices:
- Semantic colour scheme: panels are coloured by metric family (size / shape /
  structure) rather than arbitrarily, so colour carries meaning.
- Heavy right-skewed size metrics (a-c: trees, area, perimeter) on a log x-axis with
  log-spaced bins; the rest linear. This reveals the ~lognormal shape of area/perimeter.
- Clean editorial styling: despined axes, light horizontal gridlines behind the bars,
  a dashed median reference line with a value label on every panel, 'Count' only on the
  left column, embedded TrueType fonts for print.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation, MaxNLocator

# ---- global style -------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#444444",
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "axes.labelcolor": "#222222",
    "text.color": "#222222",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

# ---- palette: colour by metric family ----------------------------------------
SIZE      = "#2F6690"   # steel blue   - extent / size
SHAPE     = "#C76C3F"   # terracotta   - shape / form
STRUCTURE = "#5E8C6A"   # sage green   - internal structure / complexity
MEDLINE   = "#23303A"   # near-black for the median reference

df = pd.read_csv("data/clusters_with_metrics.csv")

# (col, label, log_x, colour, family)
panels = [
    ("n_points",          "Number of trees",   True,  SIZE,      "size"),
    ("area_m2",           "Area (m$^2$)",       True,  SIZE,      "size"),
    ("perimeter_m",       "Perimeter (m)",      True,  SIZE,      "size"),
    ("compactness",       "Form factor",        False, SHAPE,     "shape"),
    ("aspect_ratio",      "Elongation ratio",   False, SHAPE,     "shape"),
    ("fractal_dim",       "Fractal dimension",  False, STRUCTURE, "structure"),
    ("clark_evans",       "Clark–Evans index",  False, STRUCTURE, "structure"),
    ("alpha_compactness", "Alpha compactness",  False, STRUCTURE, "structure"),
]


def fmt_med(v):
    return f"{v:,.0f}" if v >= 10 else f"{v:.2f}"


fig, axes = plt.subplots(2, 4, figsize=(15, 6.6), dpi=150)

for i, (ax, (col, label, logx, color, fam)) in enumerate(zip(axes.flat, panels)):
    data = df[col].dropna()
    if logx:
        data = data[data > 0]
        bins = np.logspace(np.log10(data.min()), np.log10(data.max()), 38)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(LogFormatterSciNotation())
    else:
        bins = 38

    # gridlines behind bars
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.7, zorder=0)

    ax.hist(data, bins=bins, color=color, edgecolor="white",
            linewidth=0.35, alpha=0.92, zorder=3)

    # median reference line + label
    med = data.median()
    ax.axvline(med, color=MEDLINE, linestyle=(0, (4, 2)), linewidth=1.1, zorder=4)
    ax.annotate(f"median {fmt_med(med)}", xy=(med, 0.97), xycoords=("data", "axes fraction"),
                xytext=(4, -2), textcoords="offset points", ha="left", va="top",
                fontsize=8, color=MEDLINE,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.75))

    # despine
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3, width=0.8)

    if not logx:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))

    # panel letter
    ax.text(-0.02, 1.07, f"({chr(97 + i)})", transform=ax.transAxes,
            fontsize=11.5, fontweight="bold", va="top", ha="right", color="#222222")
    ax.set_xlabel(label + ("  (log)" if logx else ""), labelpad=3)
    if i % 4 == 0:
        ax.set_ylabel("Number of clusters")
    else:
        ax.set_ylabel("")

# family legend (color key)
from matplotlib.patches import Patch
legend_handles = [
    Patch(facecolor=SIZE, label="Size / extent"),
    Patch(facecolor=SHAPE, label="Shape / form"),
    Patch(facecolor=STRUCTURE, label="Internal structure"),
]
fig.legend(handles=legend_handles, loc="lower center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.02), fontsize=9, handlelength=1.2, columnspacing=2.2)

fig.tight_layout(rect=[0, 0.03, 1, 1], w_pad=1.6, h_pad=2.4)
out = "figures/fig2_histograms.png"
fig.savefig(out, dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig("figures/fig2_histograms.pdf", bbox_inches="tight", facecolor="white")
print(f"n = {len(df):,} clusters; saved {out} (400 dpi) + .pdf (vector)")
