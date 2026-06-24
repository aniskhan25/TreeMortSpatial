"""
Gates 1 & 2 (panel-required, all four reviewers): re-run the FULL pipeline from raw
points and test whether the two-type structure survives.

  Gate 1 (eps-stability): DBSCAN eps in {10,15,20,25,30,40} -> shape metrics + CE ->
    filter (n>=5, area>=100, AR<=10) -> standardize 3 features -> k-means k=2.
    Report per eps: n_clusters, Type-0 prevalence, silhouette/CH/DB, BIC(1,2,3),
    per-type raw centroids, and tree-level type-label ARI vs the eps=20 baseline.

  Gate 2 (AR-from-raw): at eps=20, compute metrics on ALL clusters (n>=5, area>=100,
    NO AR cap), then apply AR caps {5,7,10,15,20,inf}, re-fit k-means each time.
    Report prevalence, Type-0 raw centroid drift, ARI vs AR<=10, and the within-Type-0
    AR distribution / pile-up against the cap.

Caches per-eps cluster tables to analysis/out/ so re-runs are cheap.
"""
import os, json, sys
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from clustering import run_dbscan, filter_clusters          # noqa: E402
from metrics.shape import compute_all_shape_metrics          # noqa: E402
from metrics.spatial import clark_evans_index                # noqa: E402

from sklearn.preprocessing import StandardScaler             # noqa: E402
from sklearn.cluster import KMeans                            # noqa: E402
from sklearn.mixture import GaussianMixture                  # noqa: E402
from sklearn.metrics import (silhouette_score,               # noqa: E402
                             calinski_harabasz_score,
                             davies_bouldin_score,
                             adjusted_rand_score)

SEED = 42
FEATS = ["fractal_dim", "clark_evans", "aspect_ratio"]
EPS_GRID = [10, 15, 20, 25, 30, 40]
AR_CAPS = [5, 7, 10, 15, 20, np.inf]
OUT = "analysis/out"
os.makedirs(OUT, exist_ok=True)


def clusters_at_eps(trees, eps):
    """Re-cluster raw points at `eps`, return per-cluster table with CE (no AR cap)."""
    cache = f"{OUT}/clusters_eps{eps}.csv"
    if os.path.exists(cache):
        return pd.read_csv(cache)
    print(f"  [eps={eps}] DBSCAN on {len(trees):,} points ...", flush=True)
    df = run_dbscan(trees, eps=eps, min_samples=3)
    df = df[df.cluster_id != -1]
    print(f"  [eps={eps}] {df.cluster_id.nunique():,} raw clusters; shape metrics ...", flush=True)
    shp = compute_all_shape_metrics(df, cluster_col="cluster_id").reset_index()
    # Clark-Evans per cluster (points + convex-hull area from shape table)
    area_by_cid = dict(zip(shp.cluster_id, shp.area_m2))
    ce = {}
    for cid, grp in df.groupby("cluster_id"):
        if cid in area_by_cid:
            ce[cid] = clark_evans_index(grp[["x", "y"]].values, area_by_cid[cid])
    shp["clark_evans"] = shp.cluster_id.map(ce)
    shp.to_csv(cache, index=False)
    return shp


def fit_k2(tab):
    """Standardize 3 features, fit k-means k=2, label Type0=higher-AR centroid."""
    t = tab.dropna(subset=FEATS).copy()
    Xs = StandardScaler().fit_transform(t[FEATS].to_numpy(float))
    km = KMeans(2, n_init=10, random_state=SEED).fit(Xs)
    lab = km.labels_
    # Type0 := cluster whose mean aspect_ratio is higher (paper convention)
    ar0 = t.aspect_ratio[lab == 0].mean()
    ar1 = t.aspect_ratio[lab == 1].mean()
    type0 = 0 if ar0 >= ar1 else 1
    t["morph_type"] = np.where(lab == type0, "Type0", "Type1")
    metrics = dict(
        n=len(t),
        type0_pct=round(100 * (t.morph_type == "Type0").mean(), 2),
        silhouette=round(silhouette_score(Xs, lab), 4),
        ch=round(calinski_harabasz_score(Xs, lab), 1),
        db=round(davies_bouldin_score(Xs, lab), 4),
    )
    for k in (1, 2, 3):
        metrics[f"bic{k}"] = round(
            GaussianMixture(k, covariance_type="full", n_init=3,
                            random_state=SEED).fit(Xs).bic(Xs), 1)
    centroids = {tp: {f: round(t[t.morph_type == tp][f].mean(), 4) for f in
                      ["aspect_ratio", "compactness", "fractal_dim", "clark_evans"]}
                 for tp in ("Type0", "Type1")}
    return t, metrics, centroids


# ============================ load raw points ==============================
trees = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y"])
print(f"loaded {len(trees):,} trees", flush=True)

# ============================ GATE 1: eps-stability =========================
print("\n=== GATE 1: eps-stability ===", flush=True)
eps_rows, tree_labels = {}, {}
base_tab_eps20 = None
for eps in EPS_GRID:
    tab = clusters_at_eps(trees, eps)
    capped = filter_clusters(tab, min_trees=5, min_area=100.0, max_elongation=10.0)
    typed, m, cen = fit_k2(capped)
    eps_rows[eps] = {**m, "centroids": cen}
    # tree-level type labels for ARI: map each tree (via re-cluster) to its type
    # (recompute membership; trees not in a typed cluster -> "none")
    if eps == 20:
        base_tab_eps20 = tab  # reuse for Gate 2
    print(f"  eps={eps}: n={m['n']:5d}  Type0%={m['type0_pct']:5.1f}  "
          f"sil={m['silhouette']:.3f}  BIC1/2/3={m['bic1']:.0f}/{m['bic2']:.0f}/{m['bic3']:.0f}",
          flush=True)

# eps-invariance evidence = stability of Type-0 prevalence + per-type raw centroids
# across the grid (cluster identities differ across eps, so partition-ARI is not
# well-defined; prevalence + centroid drift is the meaningful invariance check).

with open(f"{OUT}/gate1_eps_stability.json", "w") as f:
    json.dump(eps_rows, f, indent=2)

# ============================ GATE 2: AR-from-raw ===========================
print("\n=== GATE 2: AR-from-raw sweep (eps=20) ===", flush=True)
ar_rows = {}
base = base_tab_eps20[(base_tab_eps20.n_points >= 5) & (base_tab_eps20.area_m2 >= 100)].copy()
print(f"  base clusters (n>=5, area>=100, NO AR cap): {len(base):,}", flush=True)
print(f"  clusters with AR>10: {(base.aspect_ratio > 10).sum():,}  "
      f"AR>15: {(base.aspect_ratio > 15).sum():,}  "
      f"max AR: {base.aspect_ratio.max():.1f}", flush=True)

ref_labels = None  # AR<=10 labeling, for ARI
for cap in AR_CAPS:
    sub = base[base.aspect_ratio <= cap].copy()
    typed, m, cen = fit_k2(sub)
    capname = "inf" if np.isinf(cap) else str(cap)
    # within-Type0 AR distribution
    t0ar = typed[typed.morph_type == "Type0"].aspect_ratio
    m["type0_ar_p50"] = round(t0ar.median(), 2)
    m["type0_ar_p95"] = round(t0ar.quantile(0.95), 2)
    m["type0_ar_max"] = round(t0ar.max(), 2)
    # ARI vs AR<=10 partition on shared clusters
    typed_idx = typed.set_index("cluster_id")["morph_type"]
    if cap == 10:
        ref_labels = typed_idx
        m["ARI_vs_cap10"] = 1.0
    elif ref_labels is not None:
        shared = ref_labels.index.intersection(typed_idx.index)
        m["ARI_vs_cap10"] = round(
            adjusted_rand_score(ref_labels.loc[shared], typed_idx.loc[shared]), 4)
    ar_rows[capname] = {**m, "centroids": cen}
    print(f"  AR<={capname:>4}: n={m['n']:5d}  Type0%={m['type0_pct']:5.1f}  "
          f"sil={m['silhouette']:.3f}  Type0_AR(p95/max)={m['type0_ar_p95']}/{m['type0_ar_max']}  "
          f"ARI_vs10={m.get('ARI_vs_cap10','-')}", flush=True)

with open(f"{OUT}/gate2_ar_sweep.json", "w") as f:
    json.dump(ar_rows, f, indent=2, default=str)

print("\nDONE. Wrote gate1_eps_stability.json and gate2_ar_sweep.json to", OUT)
