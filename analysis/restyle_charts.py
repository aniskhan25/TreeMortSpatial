"""
Restyle the notebook-generated CHART figures to the house style (figstyle.py),
and update stale type labels/colours to the gradient framing:
  Type0 -> terracotta, "elongated/dispersed";  Type1 -> blue, "compact/aggregated".

Figures: fig3_correlation, fig4_scatters, fig_orientation_rose, fig_shape_space,
         fig_size_dist, fig_moran_heatmap.  Reads data/clusters_typed.csv.
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import seaborn as sns
from figstyle import apply_style, PALETTE, despine, panel_label, save

apply_style()
P = PALETTE
TYPE_COLOR  = {"Type0": P["secondary"], "Type1": P["primary"]}
TYPE_MARKER = {"Type0": "^",            "Type1": "o"}
TYPE_LABEL  = {"Type0": "Type 0 (elongated/dispersed)",
               "Type1": "Type 1 (compact/aggregated)"}

df = pd.read_csv("data/clusters_typed.csv")


# ---- fig3: correlation matrix -------------------------------------------------
def fig_correlation():
    cols = ["n_points", "area_m2", "perimeter_m", "compactness",
            "aspect_ratio", "fractal_dim", "clark_evans", "alpha_compactness"]
    labels = ["n trees", "area", "perimeter", "form factor",
              "elongation", "fractal dim", "Clark–Evans", "alpha comp."]
    corr = df[cols].corr()
    mask = np.triu(np.ones_like(corr, dtype=bool))
    fig, ax = plt.subplots(figsize=(7.5, 6.2), dpi=150)
    sns.heatmap(corr, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
                vmin=-1, vmax=1, ax=ax, square=True,
                linewidths=0.6, linecolor="white",
                annot_kws={"size": 8.5},
                cbar_kws={"label": "Pearson r", "shrink": 0.75})
    ax.set_xticklabels(labels, rotation=40, ha="right")
    ax.set_yticklabels(labels, rotation=0)
    ax.tick_params(length=0)
    fig.tight_layout()
    save(fig, "figures/fig3_correlation")
    plt.close(fig)


# ---- fig4: pairwise scatters --------------------------------------------------
def fig_scatters():
    pairs = [
        ("clark_evans", "fractal_dim", "Clark–Evans index", "Fractal dimension"),
        ("clark_evans", "aspect_ratio", "Clark–Evans index", "Elongation ratio"),
        ("fractal_dim", "alpha_compactness", "Fractal dimension", "Alpha compactness"),
        ("fractal_dim", "compactness", "Fractal dimension", "Form factor"),
        ("fractal_dim", "perimeter_m", "Fractal dimension", "Perimeter (m)"),
        ("alpha_compactness", "perimeter_m", "Alpha compactness", "Perimeter (m)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4), dpi=150)
    for i, (ax, (xc, yc, xl, yl)) in enumerate(zip(axes.flat, pairs)):
        despine(ax, grid_axis="both")
        for t in ("Type1", "Type0"):   # plot Type1 first so Type0 sits on top
            sub = df[df.morph_type == t][[xc, yc]].dropna()
            ax.scatter(sub[xc], sub[yc], s=4, alpha=0.18,
                       color=TYPE_COLOR[t], edgecolors="none", rasterized=True)
        r = df[[xc, yc]].dropna().corr().iloc[0, 1]
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.text(0.04, 0.95, f"r = {r:.2f}", transform=ax.transAxes,
                va="top", ha="left", fontsize=9, color="#333",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
        panel_label(ax, chr(97 + i))
        if yc == "perimeter_m":
            ax.set_yscale("log")
    handles = [plt.Line2D([], [], marker="o", ls="", color=TYPE_COLOR[t],
                          label=TYPE_LABEL[t], markersize=7) for t in ("Type0", "Type1")]
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=[0, 0.03, 1, 1], w_pad=2, h_pad=2.2)
    save(fig, "figures/fig4_scatters")
    plt.close(fig)


# ---- fig_orientation_rose -----------------------------------------------------
def fig_orientation():
    import pycircstat
    d = df.dropna(subset=["morph_type", "orientation"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2),
                             subplot_kw={"projection": "polar"}, dpi=150)
    for ax, t in zip(axes, ["Type0", "Type1"]):
        sub = d[d.morph_type == t]
        doubled = np.deg2rad((2 * sub.orientation.values) % 360)
        p_val, _ = pycircstat.tests.rayleigh(doubled)
        r_bar = pycircstat.resultant_vector_length(doubled)
        mean_dir = (np.rad2deg(pycircstat.mean(doubled)) % 360) / 2 % 180
        mirror = np.deg2rad(np.concatenate([sub.orientation.values,
                                            sub.orientation.values + 180]) % 360)
        bins = np.linspace(0, 2 * np.pi, 37)
        counts, _ = np.histogram(mirror, bins=bins)
        centers = 0.5 * (bins[:-1] + bins[1:])
        ax.bar(centers, counts, width=2 * np.pi / 36, color=TYPE_COLOR[t],
               alpha=0.85, edgecolor="white", lw=0.3, zorder=3)
        ax.axvline(np.deg2rad(mean_dir), color=P["neutral"], lw=1.8,
                   label=f"mean axis {mean_dir:.0f}°")
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
        ax.set_title(f"{TYPE_LABEL[t]}\n$r$={r_bar:.3f},  $p$={p_val:.2g}",
                     fontsize=10, pad=16)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=8, loc="lower right", bbox_to_anchor=(1.18, -0.05))
    fig.tight_layout(w_pad=3)
    save(fig, "figures/fig_orientation_rose")
    plt.close(fig)


# ---- fig_shape_space (zones now illustrative / author-defined) ----------------
def fig_shape_space():
    d = df.dropna(subset=["morph_type", "aspect_ratio", "compactness"])
    zones = [
        {"label": "elongated (wind-like)",   "x": (3.5, 10.0), "y": (0.05, 0.45), "c": P["secondary"]},
        {"label": "intermediate (beetle-like)", "x": (1.0, 5.0), "y": (0.25, 0.70), "c": P["neutral"]},
        {"label": "compact (fire-like)",     "x": (1.0, 3.5),  "y": (0.40, 0.85), "c": P["tertiary"]},
    ]
    fig, ax = plt.subplots(figsize=(9, 5.6), dpi=150)
    despine(ax, grid_axis="both")
    for z in zones:
        ax.add_patch(Rectangle((z["x"][0], z["y"][0]), z["x"][1]-z["x"][0], z["y"][1]-z["y"][0],
                     lw=1.2, edgecolor=z["c"], facecolor=z["c"], alpha=0.08, zorder=1, ls="--"))
        ax.text(np.mean(z["x"]), z["y"][1]+0.012, z["label"], ha="center", va="bottom",
                color=z["c"], fontsize=8, fontstyle="italic")
    for t in ("Type1", "Type0"):
        sub = d[d.morph_type == t]
        ax.scatter(sub.aspect_ratio, sub.compactness, color=TYPE_COLOR[t],
                   s=6, alpha=0.22, edgecolors="none", label=TYPE_LABEL[t],
                   zorder=2, rasterized=True)
    ax.set_xlim(0.5, 11); ax.set_ylim(0, 1)
    ax.set_xlabel("Aspect ratio (elongation)"); ax.set_ylabel("Compactness (form factor)")
    leg = ax.legend(loc="upper right", markerscale=2)
    for lh in leg.legend_handles:
        lh.set_alpha(1)
    ax.text(0.015, 0.03, "Reference regions are illustrative and author-defined (see caption).",
            transform=ax.transAxes, fontsize=7.5, color="#777", style="italic")
    fig.tight_layout()
    save(fig, "figures/fig_shape_space")
    plt.close(fig)


# ---- fig_size_dist (CCDF; lognormal preferred) --------------------------------
def fig_size_dist():
    import powerlaw as pl
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), dpi=150)
    for ax, t in zip(axes, ["Type0", "Type1"]):
        despine(ax, grid_axis="both")
        data = df[df.morph_type == t]["n_points"].dropna().values.astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = pl.Fit(data, xmin=5, discrete=True, verbose=False)
            R, p_val = fit.distribution_compare("power_law", "lognormal")
            fit.plot_ccdf(ax=ax, color=TYPE_COLOR[t], lw=2.2, label="empirical CCDF")
            fit.lognormal.plot_ccdf(ax=ax, color=P["neutral"], lw=1.5, ls="-",
                                    label="lognormal (preferred)")
            fit.power_law.plot_ccdf(ax=ax, color="#999", lw=1.4, ls="--",
                                    label=f"power law (α={fit.power_law.alpha:.2f})")
        ax.set_xlabel("Cluster size ($n$ trees)"); ax.set_ylabel(r"$P(X \geq x)$")
        ax.set_title(f"{TYPE_LABEL[t]}  ($n$={len(data):,})\n"
                     f"lognormal favoured (LR $R$={R:.2f}, $p$={p_val:.3f})", fontsize=9.5)
        ax.legend(loc="lower left")
    panel_label(axes[0], "a"); panel_label(axes[1], "b")
    fig.tight_layout(w_pad=2.4)
    save(fig, "figures/fig_size_dist")
    plt.close(fig)


# ---- fig_moran_heatmap --------------------------------------------------------
def fig_moran():
    from libpysal.weights import DistanceBand
    from esda.moran import Moran
    THR = [200, 500, 1000, 2000, 5000]
    MET = ["fractal_dim", "clark_evans", "aspect_ratio", "compactness"]
    MLAB = ["Fractal dim.", "Clark–Evans", "Aspect ratio", "Compactness"]
    N = 2000
    t0, t1 = df[df.morph_type == "Type0"], df[df.morph_type == "Type1"]
    rng = np.random.default_rng(42)
    n0 = int(N * len(t0) / len(df)); n1 = N - n0
    dm = pd.concat([t0.iloc[rng.choice(len(t0), min(n0, len(t0)), replace=False)],
                    t1.iloc[rng.choice(len(t1), min(n1, len(t1)), replace=False)]]).reset_index(drop=True)
    coords = dm[["cx", "cy"]].values
    alpha_b = 0.05 / (len(MET) * len(THR))
    I = np.full((len(MET), len(THR)), np.nan); pv = np.full_like(I, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for j, th in enumerate(THR):
            w = DistanceBand.from_array(coords, threshold=th, silence_warnings=True)
            w.transform = "r"
            for i, col in enumerate(MET):
                # analytical p (normal approx): a permutation p_sim floor (>=1e-3)
                # cannot resolve the Bonferroni alpha* = 0.0025, so use p_norm.
                m = Moran(dm[col].fillna(dm[col].median()).values, w, permutations=0)
                I[i, j], pv[i, j] = m.I, m.p_norm
    fig, ax = plt.subplots(figsize=(8.5, 3.8), dpi=150)
    im = ax.imshow(I, aspect="auto", cmap="RdBu_r", vmin=-0.3, vmax=0.3)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02); cb.set_label("Moran's $I$")
    cb.outline.set_visible(False)
    for i in range(len(MET)):
        for j in range(len(THR)):
            p = pv[i, j]; sig = "*" if p < alpha_b else ("†" if p < 0.05 else "")
            ax.text(j, i, f"{I[i,j]:.2f}{sig}", ha="center", va="center", fontsize=9.5,
                    color="white" if abs(I[i, j]) > 0.15 else "#222")
    ax.set_xticks(range(len(THR))); ax.set_xticklabels([f"{t} m" for t in THR])
    ax.set_yticks(range(len(MET))); ax.set_yticklabels(MLAB)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"Global Moran's $I$  (*$p<${alpha_b:.4f} Bonferroni; † $p<$0.05)",
                 fontsize=10, pad=8)
    fig.tight_layout()
    save(fig, "figures/fig_moran_heatmap")
    plt.close(fig)


if __name__ == "__main__":
    fig_correlation(); print("fig3_correlation")
    fig_scatters();    print("fig4_scatters")
    fig_orientation(); print("fig_orientation_rose")
    fig_shape_space(); print("fig_shape_space")
    fig_size_dist();   print("fig_size_dist")
    fig_moran();       print("fig_moran_heatmap")
    print("charts done")
