"""
Sampling-design follow-ups (Vasquez + Anand):
(1) Re-run the held-out RF coherence check with the 7 real super-site BLOCKS as CV groups
    (leave-one-block-out), replacing the arbitrary 100 km grid — the correct leakage guard.
(2) Per-block coverage table (n_tiles proxy via clusters, latitude, Type-0%, within-block
    compactness–AR correlation).
(3) Block-level sign test for the connectivity-reach direction (pseudoreplication-honest:
    7 blocks, not 14,582 clusters).
"""
import warnings
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
from scipy.stats import binomtest
warnings.filterwarnings("ignore")
SEED = 42

d = pd.read_csv("data/clusters_typed.csv")
d["block"] = DBSCAN(eps=15000, min_samples=1).fit_predict(d[["cx", "cy"]].values)
def lat(y): return 60 + (y - 6660000) / 111320.0
y = (d.morph_type == "Type0").astype(int).to_numpy()

# ---- (1) block-grouped (leave-one-block-out) RF CV ---------------------------
feats = ["compactness", "alpha_compactness", "core_to_edge_ratio",
         "nnd_mean", "nnd_std", "orientation", "reach_200m"]
X = d[feats].fillna(d[feats].median()).to_numpy()
Xs = StandardScaler().fit_transform(X)
groups = d.block.values
logo = LeaveOneGroupOut()
aucs = []
print("=== (1) Held-out RF, leave-one-BLOCK-out CV (7 folds) ===")
for tr, te in logo.split(Xs, y, groups):
    if len(np.unique(y[te])) < 2:
        continue
    rf = RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=-1).fit(Xs[tr], y[tr])
    a = roc_auc_score(y[te], rf.predict_proba(Xs[te])[:, 1])
    b = d.block.values[te][0]
    aucs.append(a)
    print(f"  held-out block {b} (lat {lat(d[d.block==b].cy.mean()):.1f}N, n={ (groups==b).sum() }): AUC={a:.3f}")
print(f"  mean AUC = {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}   (100km-grid CV was 0.976)")

# ---- (2) per-block coverage table -------------------------------------------
print("\n=== (2) per-block coverage ===")
print(f"{'blk':>3} {'n_clus':>7} {'lat':>5} {'Type0%':>7} {'r(AR,comp)':>11}")
for b, g in d.groupby("block"):
    r = g[["aspect_ratio", "compactness"]].corr().iloc[0, 1]
    print(f"{b:>3} {len(g):>7} {lat(g.cy.mean()):>5.1f} {100*(g.morph_type=='Type0').mean():>6.1f}% {r:>11.3f}")

# ---- (3) reach direction: block-level sign test -----------------------------
print("\n=== (3) connectivity reach: block-level sign test (Type1>Type0?) ===")
dirs = []
for b, g in d.groupby("block"):
    m0 = g[g.morph_type == "Type0"].reach_200m.median()
    m1 = g[g.morph_type == "Type1"].reach_200m.median()
    dirs.append(m1 > m0)
    print(f"  block {b}: Type1 med reach={m1:.0f} vs Type0={m0:.0f}  -> {'Type1>Type0' if m1>m0 else 'no'}")
k = sum(dirs); n = len(dirs)
p = binomtest(k, n, 0.5).pvalue
print(f"  {k}/{n} blocks show Type1>Type0; sign test p={p:.3f}")
print("  (honest test on 7 PSUs, vs the pseudoreplicated Mann-Whitney p=1.7e-10 on 14,582 clusters)")
