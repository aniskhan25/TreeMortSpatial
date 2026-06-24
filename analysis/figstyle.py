"""Shared house style for the custom analysis figures (matches fig2_histograms).

Usage:
    from figstyle import apply_style, PALETTE, despine, panel_label, save
    apply_style()
    ...
    save(fig, "figures/fig_x")   # writes .png (400 dpi) + .pdf (vector)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Semantic / role palette (same hues as fig2_histograms)
PALETTE = {
    "primary":   "#2F6690",   # steel blue   - main distribution / size
    "secondary": "#C76C3F",   # terracotta   - secondary series / shape
    "tertiary":  "#5E8C6A",   # sage green   - structure
    "accent":    "#C2452D",   # alert red    - observed value / emphasis
    "neutral":   "#23303A",   # near-black   - reference lines
    "muted":     "#9AB0BE",   # light blue   - envelopes / null mass
    "grid":      "#E8E8E8",
}


def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 10,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#444444",
        "axes.titlesize": 10.5,
        "axes.labelsize": 10,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.color": "#444444",
        "ytick.color": "#444444",
        "axes.labelcolor": "#222222",
        "text.color": "#222222",
        "legend.fontsize": 8.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def despine(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3, width=0.8)
    ax.set_axisbelow(True)
    if grid_axis:
        ax.grid(axis=grid_axis, color=PALETTE["grid"], linewidth=0.7, zorder=0)


def panel_label(ax, letter):
    ax.text(-0.02, 1.08, f"({letter})", transform=ax.transAxes,
            fontsize=11.5, fontweight="bold", va="top", ha="right", color="#222222")


def save(fig, stem):
    fig.savefig(f"{stem}.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
