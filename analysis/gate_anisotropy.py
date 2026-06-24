"""
Gate 3 (Vasquez, required): does the wind-vs-beetle distinction have any ANISOTROPIC
second-order content, or is "elongated Type 0" just hull geometry?

The first-order pipeline (hull aspect ratio) cannot see whether points are arranged
ALONG a line (wind swath) vs filling an elongated blob (a patch that merely happens to
be non-round). We test second-order linearity on the within-cluster points:

  - per-cluster PCA eigenvalues l1>=l2 of the point scatter:
      scatter_AR = sqrt(l1/l2)         (point-level elongation)
      linearity  = l1 / (l1 + l2)      (0.5 = isotropic, ->1 = collinear/line)
  - pooled pairwise-displacement axial resultant after aligning each cluster to its
    own major axis: high resultant = displacements concentrate along the axis = line-like.

Critically, we compare Type0 vs Type1 BOTH overall AND within matched hull-AR bins, so
the comparison is not just restating that Type0 has higher hull AR by construction.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from scipy.stats import mannwhitneyu                          # noqa: E402

SEED = 42

typed = pd.read_csv("data/clusters_typed.csv")[
    ["cluster_id", "morph_type", "aspect_ratio"]]
trees = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id"])
trees = trees[trees.cluster_id.isin(set(typed.cluster_id))]

type_of = dict(zip(typed.cluster_id, typed.morph_type))
hullar_of = dict(zip(typed.cluster_id, typed.aspect_ratio))

rows = []
for cid, grp in trees.groupby("cluster_id"):
    P = grp[["x", "y"]].to_numpy(float)
    if len(P) < 3:
        continue
    Pc = P - P.mean(0)
    cov = np.cov(Pc, rowvar=False)
    w, V = np.linalg.eigh(cov)          # ascending
    l2, l1 = max(w[0], 1e-9), max(w[1], 1e-9)
    major = V[:, 1]
    # pairwise displacement angles relative to the major axis (axial, doubled)
    # resultant of 2*theta: ->1 means displacements lie along one axis (linear)
    n = len(Pc)
    if n >= 3:
        idx = np.triu_indices(n, 1)
        d = Pc[idx[0]] - Pc[idx[1]]
        ang = np.arctan2(d @ np.array([-major[1], major[0]]), d @ major)  # rel. major
        R_axial = np.abs(np.mean(np.exp(2j * ang)))
    else:
        R_axial = np.nan
    rows.append(dict(cluster_id=cid, morph_type=type_of[cid], n=n,
                     hull_ar=hullar_of[cid],
                     scatter_ar=np.sqrt(l1 / l2),
                     linearity=l1 / (l1 + l2),
                     pair_axial_R=R_axial))

an = pd.DataFrame(rows)
t0, t1 = an[an.morph_type == "Type0"], an[an.morph_type == "Type1"]

def cmp(col):
    a, b = t0[col].dropna(), t1[col].dropna()
    U, p = mannwhitneyu(a, b)
    # rank-biserial effect size
    rbc = 1 - 2 * U / (len(a) * len(b))
    return dict(metric=col,
                Type0_median=round(a.median(), 3), Type1_median=round(b.median(), 3),
                Type0_mean=round(a.mean(), 3), Type1_mean=round(b.mean(), 3),
                rank_biserial=round(rbc, 3), p=f"{p:.2e}")

print(f"clusters analysed: {len(an):,}  (Type0={len(t0):,}, Type1={len(t1):,})\n")
print("=== second-order linearity by type (overall) ===")
for c in ["scatter_ar", "linearity", "pair_axial_R"]:
    r = cmp(c)
    print(f"  {r['metric']:13s} Type0={r['Type0_median']:.3f} Type1={r['Type1_median']:.3f} "
          f"  rank-biserial={r['rank_biserial']:+.3f}  p={r['p']}")

# AR-matched comparison: within hull-AR bins, is Type0 still more linear?
print("\n=== AR-matched: linearity Type0 vs Type1 within hull-AR bins ===")
bins = [1, 1.5, 2, 2.5, 3, 4, 10]
an["arbin"] = pd.cut(an.hull_ar, bins)
matched = []
for b, g in an.groupby("arbin", observed=True):
    g0, g1 = g[g.morph_type == "Type0"].linearity, g[g.morph_type == "Type1"].linearity
    if len(g0) >= 20 and len(g1) >= 20:
        U, p = mannwhitneyu(g0, g1)
        matched.append((str(b), len(g0), len(g1), round(g0.median(), 3),
                        round(g1.median(), 3), f"{p:.1e}"))
        print(f"  AR {str(b):12s} n0={len(g0):4d} n1={len(g1):5d}  "
              f"linearity T0={g0.median():.3f} T1={g1.median():.3f}  p={p:.1e}")
if not matched:
    print("  (no AR bins with both types well-populated -> Type0/Type1 barely overlap in AR,"
          "\n   which itself means 'type' is essentially a hull-AR threshold, not a 2nd-order class)")

an.to_csv("analysis/out/gate3_anisotropy.csv", index=False)
print("\nwrote analysis/out/gate3_anisotropy.csv")
print("\nINTERPRETATION: if Type0 is NOT meaningfully more linear than Type1 (small effect"
      "\nsize, or no difference within AR-matched bins), the elongated 'type' carries no"
      "\ndirectional/2nd-order signal beyond hull shape -> the wind label is unsupported.")
