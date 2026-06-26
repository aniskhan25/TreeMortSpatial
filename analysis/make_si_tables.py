"""Supplementary Information tables/figure: shape-feature correlation matrix and univariate AUCs.
Descriptive diagnostics that support the morphological-gradient framing (correlated, low-AUC features
=> a continuum the discrete labels summarise). Writes:
  analysis/out/si_feature_corr.csv      Spearman correlation matrix
  analysis/out/si_univariate_auc.csv    per-feature AUC for separating the two morphological endpoints
  figures/fig_si_corr.{png,pdf}         correlation heatmap
and prints booktabs LaTeX for the AUC table (pasted into report/supplementary.tex).
"""
import os, sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
sys.path.insert(0, "analysis")
from figstyle import apply_style, PALETTE, despine, save

OUT = "analysis/out"; os.makedirs(OUT, exist_ok=True)
KMEANS = ["fractal_dim", "clark_evans", "aspect_ratio"]          # the features that DEFINE the partition
FEATURES = ["aspect_ratio", "elongation", "compactness", "alpha_compactness", "core_to_edge_ratio",
            "fractal_dim", "clark_evans", "area_m2", "perimeter_m", "n_points", "nnd_mean"]
LABEL = {"aspect_ratio": "aspect ratio", "elongation": "elongation", "compactness": "compactness",
         "alpha_compactness": r"$\alpha$-compactness", "core_to_edge_ratio": "core:edge",
         "fractal_dim": "fractal dim", "clark_evans": "Clark--Evans", "area_m2": "area",
         "perimeter_m": "perimeter", "n_points": "$n$ points", "nnd_mean": "mean NND"}

d = pd.read_csv("data/clusters_typed.csv")
y = (d.morph_type == "Type1").astype(int).values                # 1 = compact endpoint (Type1, 76.5%)
X = d[FEATURES].copy()

# --- Spearman correlation matrix ---
C = X.corr(method="spearman")
C.to_csv(f"{OUT}/si_feature_corr.csv")

# --- univariate AUC (orient so AUC >= 0.5; record sign) ---
rows = []
for f in FEATURES:
    m = np.isfinite(X[f].values)
    a = roc_auc_score(y[m], X[f].values[m]); sign = "+" if a >= 0.5 else "-"
    rows.append((f, max(a, 1 - a), sign, int((~m).sum()), f in KMEANS))
auc = pd.DataFrame(rows, columns=["feature", "auc", "direction_with_compact", "n_missing", "kmeans_input"])
auc = auc.sort_values("auc", ascending=False).reset_index(drop=True)
auc.to_csv(f"{OUT}/si_univariate_auc.csv", index=False)

# --- heatmap ---
apply_style()
fig, ax = plt.subplots(figsize=(7.2, 6.2))
order = FEATURES
M = C.loc[order, order].values
im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
ax.set_xticklabels([LABEL[f] for f in order], rotation=45, ha="right")
ax.set_yticklabels([LABEL[f] for f in order])
for i in range(len(order)):
    for j in range(len(order)):
        v = M[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color="white" if abs(v) > 0.55 else "#222222", fontsize=7)
# mark the three K-means-input features
for f in KMEANS:
    k = order.index(f)
    ax.add_patch(plt.Rectangle((k - .5, -.5), 1, len(order), fill=False, ec=PALETTE["accent"], lw=1.6))
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("Spearman $\\rho$")
ax.set_title("Shape-feature correlation (red box = K-means input)")
fig.tight_layout()
save(fig, "figures/fig_si_corr")

# --- LaTeX for the AUC table ---
print("clusters:", len(d), "| Type1(compact):", int(y.sum()), "Type0(elongated):", int((1-y).sum()))
print("\n% --- univariate AUC table (paste into supplementary.tex) ---")
for _, r in auc.iterrows():
    star = r"$^\dagger$" if r.kmeans_input else ""
    print(f"{LABEL[r.feature]}{star} & {r.auc:.3f} & {r.direction_with_compact} \\\\")
print("\nsaved si_feature_corr.csv, si_univariate_auc.csv, figures/fig_si_corr.*")
