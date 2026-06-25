"""
Tanaka's gate: propagate detection error into the Moran's I autocorrelation statistics
(not just the type labels). Does the multi-scale spatial-autocorrelation pattern survive
random false-negative detection perturbation?

Design:
- Same stratified 2,000-cluster subsample as the paper (seed 42).
- Baseline = gate metrics recomputed from the FULL points of those clusters.
- N realizations: drop 15% of each cluster's detection points (FN model; keep >=5),
  recompute the 4 Moran metrics + centroid, recompute Moran's I at 5 thresholds.
- Analytical Moran p-value (p_norm) is used throughout: the paper's permutations=199 has a
  p_sim floor of 0.005 and literally cannot reach the Bonferroni alpha*=0.0025, so the
  permutation p is the wrong instrument for a Bonferroni claim. p_norm is continuous.
- PASS (per metric x threshold): sign of I preserved AND p_norm < 0.0025 in >=95% of reps.
  Focus on the 200-500 m scales the paper interprets.

A lean fractal estimate (coarser grid) is used for tractability; baseline and perturbed use
the SAME estimator, so sign/significance comparisons are internally consistent.
"""
import sys, warnings
import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint, box as sbox
from scipy.stats import linregress
from libpysal.weights import DistanceBand
from esda.moran import Moran

sys.path.insert(0, "src")
from metrics.spatial import clark_evans_index  # noqa: E402

warnings.filterwarnings("ignore")
SEED = 42
THR = [200, 500, 1000, 2000, 5000]
MET = ["fractal_dim", "clark_evans", "aspect_ratio", "compactness"]
ALPHA_B = 0.05 / (len(MET) * len(THR))     # 0.0025
N_REP = 100
FRAC = 0.15

typed = pd.read_csv("data/clusters_typed.csv")
# stratified 2000 subsample, same scheme as the Moran heatmap
t0 = typed[typed.morph_type == "Type0"]; t1 = typed[typed.morph_type == "Type1"]
rng = np.random.default_rng(SEED)
n0 = int(2000 * len(t0) / len(typed)); n1 = 2000 - n0
sub_ids = pd.concat([t0.iloc[rng.choice(len(t0), n0, replace=False)],
                     t1.iloc[rng.choice(len(t1), n1, replace=False)]]).cluster_id.values
sub_set = set(sub_ids)
print(f"subsample: {len(sub_ids)} clusters (n0={n0}, n1={n1})", flush=True)

# points for the subsample clusters
trees = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id"])
trees = trees[trees.cluster_id.isin(sub_set)]
pts_by = {cid: g[["x", "y"]].to_numpy(float) for cid, g in trees.groupby("cluster_id")}


def fractal_lean(hull):
    x0, y0, x1, y1 = hull.bounds
    diag = np.hypot(x1 - x0, y1 - y0)
    if diag == 0:
        return np.nan
    scales = np.logspace(np.log10(0.1), np.log10(0.5), 6) * diag
    counts = []
    for s in scales:
        xb = np.arange(x0, x1 + s, s); yb = np.arange(y0, y1 + s, s)
        c = sum(1 for i in range(len(xb)-1) for j in range(len(yb)-1)
                if hull.intersects(sbox(xb[i], yb[j], xb[i+1], yb[j+1])))
        counts.append(c)
    counts = np.array(counts); valid = counts > 0
    if valid.sum() < 2:
        return np.nan
    slope = linregress(np.log(scales[valid]), np.log(counts[valid]))[0]
    return -slope


def metrics(pts):
    if len(pts) < 3:
        return None
    hull = MultiPoint(pts).convex_hull
    if hull.geom_type != "Polygon":
        return None
    area, per = hull.area, hull.length
    if area <= 0 or per <= 0:
        return None
    comp = 4 * np.pi * area / per**2
    rect = hull.minimum_rotated_rectangle.exterior.coords
    sides = [np.hypot(rect[i+1][0]-rect[i][0], rect[i+1][1]-rect[i][1]) for i in range(4)]
    ar = max(sides) / max(min(sides), 1.0)
    ce = clark_evans_index(pts, area)
    fd = fractal_lean(hull)
    return dict(fractal_dim=fd, clark_evans=ce, aspect_ratio=ar, compactness=comp,
                cx=pts[:, 0].mean(), cy=pts[:, 1].mean())


def build_frame(drop_frac, seed):
    r = np.random.default_rng(seed)
    rows = []
    for cid in sub_ids:
        p = pts_by.get(cid)
        if p is None:
            continue
        if drop_frac > 0 and len(p) > 5:
            keep = r.random(len(p)) >= drop_frac
            if keep.sum() >= 5:
                p = p[keep]
        m = metrics(p)
        if m:
            rows.append(m)
    return pd.DataFrame(rows)


def moran_grid(df):
    coords = df[["cx", "cy"]].values
    I = np.full((len(MET), len(THR)), np.nan); pn = np.full_like(I, np.nan)
    for j, th in enumerate(THR):
        w = DistanceBand.from_array(coords, threshold=th, silence_warnings=True)
        w.transform = "r"
        for i, col in enumerate(MET):
            v = df[col].fillna(df[col].median()).values
            m = Moran(v, w, permutations=0)
            I[i, j], pn[i, j] = m.I, m.p_norm
    return I, pn


# baseline
print("baseline (unperturbed)...", flush=True)
base_I, base_p = moran_grid(build_frame(0.0, 0))

# realizations
sign_ok = np.zeros((len(MET), len(THR)))
sig_ok = np.zeros((len(MET), len(THR)))
I_stack = []
for rep in range(1, N_REP + 1):
    I, pn = moran_grid(build_frame(FRAC, 1000 + rep))
    sign_ok += (np.sign(I) == np.sign(base_I)).astype(float)
    sig_ok += ((pn < ALPHA_B) & (np.sign(I) == np.sign(base_I))).astype(float)
    I_stack.append(I)
    if rep % 20 == 0:
        print(f"  {rep}/{N_REP}", flush=True)

I_stack = np.array(I_stack)
sign_pct = 100 * sign_ok / N_REP
sig_pct = 100 * sig_ok / N_REP

print(f"\n=== Moran's I detection-propagation gate (FN={int(FRAC*100)}%, N={N_REP}) ===")
print(f"Bonferroni alpha* = {ALPHA_B}")
print(f"{'metric':14s} {'thr':>6s} {'I_base':>8s} {'I_mean':>8s} {'I_sd':>7s} "
      f"{'sign%':>6s} {'sig%':>6s}")
out = {}
for i, met in enumerate(MET):
    for j, th in enumerate(THR):
        print(f"{met:14s} {th:6d} {base_I[i,j]:8.3f} {I_stack[:,i,j].mean():8.3f} "
              f"{I_stack[:,i,j].std():7.3f} {sign_pct[i,j]:6.0f} {sig_pct[i,j]:6.0f}")
        out[f"{met}@{th}"] = dict(I_base=round(float(base_I[i,j]),4),
                                  I_mean=round(float(I_stack[:,i,j].mean()),4),
                                  sign_pct=float(sign_pct[i,j]), sig_pct=float(sig_pct[i,j]))

print("\n=== VERDICT at 200-500 m (the interpreted scales) ===")
verdict_rows = [(m, t) for m in MET for t in (200, 500)]
allpass = True
for met, th in verdict_rows:
    j = THR.index(th); i = MET.index(met)
    ok = sign_pct[i, j] >= 95 and sig_pct[i, j] >= 95
    allpass &= ok
    print(f"  {met}@{th}m: sign {sign_pct[i,j]:.0f}%  sig {sig_pct[i,j]:.0f}%  -> {'PASS' if ok else 'FAIL'}")
print(f"\nGATE {'PASSED' if allpass else 'MIXED/FAILED'}: "
      f"{'autocorrelation survives detection perturbation' if allpass else 'see per-cell results'}")

import json
json.dump(out, open("analysis/out/gate_moran_propagation.json", "w"), indent=2)
print("saved analysis/out/gate_moran_propagation.json")
