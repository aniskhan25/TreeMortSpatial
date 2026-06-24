"""
fig_shap_heldout.png — restyled to the house style.

The shap beeswarm hardcodes its own red/blue colormap and layout and cannot be brought
into the house palette; its scientific content here is the *ranking* of mean |SHAP|
(compactness dominates; reach_200m, the only shape-orthogonal feature, ranks last). We
therefore render a clean horizontal mean-|SHAP| bar chart, with reach_200m emphasised.

Values are the E10 held-out-feature Random Forest results reported in the manuscript
(spatial block CV, AUC = 0.976 +/- 0.002).
"""
import matplotlib.pyplot as plt
from figstyle import apply_style, PALETTE, despine, save

apply_style()
P = PALETTE

# mean |SHAP| from NB06 E10 (held-out features), as reported in the manuscript
shap_vals = [
    ("compactness",        0.226),
    ("core_to_edge_ratio", 0.068),
    ("nnd_mean",           0.037),
    ("orientation",        0.034),
    ("alpha_compactness",  0.019),
    ("nnd_std",            0.010),
    ("reach_200m",         0.005),
]
shap_vals.sort(key=lambda x: x[1])          # ascending for barh (largest on top)
names = [n for n, _ in shap_vals]
vals = [v for _, v in shap_vals]
colors = [P["accent"] if n == "reach_200m" else P["primary"] for n in names]

fig, ax = plt.subplots(figsize=(7.6, 4.2), dpi=150)
despine(ax, grid_axis="x")
bars = ax.barh(names, vals, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
for n, v in zip(names, vals):
    ax.text(v + 0.004, n, f"{v:.3f}", va="center", ha="left", fontsize=8.5,
            color=P["accent"] if n == "reach_200m" else "#333")
ax.set_xlim(0, 0.25)
ax.set_xlabel("mean |SHAP value|  (impact on predicted type)")
ax.set_title("Held-out Random Forest feature importance (E10)",
             fontsize=11, pad=8, loc="left", fontweight="bold")
ax.text(0.99, 0.06,
        "reach_200m — the only shape-orthogonal feature — ranks last:\n"
        "the high AUC is geometric redundancy, not independent signal.",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.8,
        color="#555", style="italic")
fig.tight_layout()
save(fig, "figures/fig_shap_heldout")
print("fig_shap_heldout done")
