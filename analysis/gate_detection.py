"""
Gate 4 (Tanaka, required): is the typology partly a footprint of the unvalidated
detector, rather than ecology? Two cheap-but-pointed tests using the max_value
detection-confidence field that the paper never uses.

  R4 - confound discriminator: do detection covariates (mean detection confidence,
       cluster size, local density) PREDICT the type / gradient position? If yes,
       the "type" co-varies with detector performance.
  R3-lite - false-negative perturbation: drop a fraction of detections per cluster
       (the FN model), recompute the 3 k-means features + CE, re-fit, and measure how
       many clusters change type. Indicative only; a full Monte-Carlo needs the
       stratified precision/recall (R1/R2) which requires external ground truth.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from metrics.shape import compute_all_shape_metrics          # noqa: E402
from metrics.spatial import clark_evans_index                # noqa: E402
from sklearn.preprocessing import StandardScaler             # noqa: E402
from sklearn.cluster import KMeans                            # noqa: E402
from sklearn.linear_model import LogisticRegression          # noqa: E402
from sklearn.model_selection import cross_val_score          # noqa: E402
from sklearn.metrics import adjusted_rand_score              # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr              # noqa: E402

SEED = 42
FEATS = ["fractal_dim", "clark_evans", "aspect_ratio"]

typed = pd.read_csv("data/clusters_typed.csv")
trees = pd.read_csv("data/trees_clustered.csv",
                    usecols=["x", "y", "cluster_id", "max_value"])
trees = trees[trees.cluster_id.isin(set(typed.cluster_id))]

# per-cluster detection covariates
conf = trees.groupby("cluster_id").max_value.agg(["mean", "std"]).rename(
    columns={"mean": "conf_mean", "std": "conf_std"})
d = typed.merge(conf, on="cluster_id", how="left")
d["density"] = d.n_points / d.area_m2
y = (d.morph_type == "Type0").astype(int).to_numpy()

print("=== R4: do detection covariates predict TYPE? ===")
print(f"  Type0 conf_mean median={d[d.morph_type=='Type0'].conf_mean.median():.4f}  "
      f"Type1={d[d.morph_type=='Type1'].conf_mean.median():.4f}")
U, p = mannwhitneyu(d[d.morph_type == "Type0"].conf_mean.dropna(),
                    d[d.morph_type == "Type1"].conf_mean.dropna())
print(f"  detection-confidence differs by type? Mann-Whitney p={p:.2e}")

det = d[["conf_mean", "conf_std", "n_points", "density"]].fillna(d[
    ["conf_mean", "conf_std", "n_points", "density"]].median())
Xd = StandardScaler().fit_transform(det)
auc_det = cross_val_score(LogisticRegression(max_iter=1000), Xd, y,
                          cv=5, scoring="roc_auc")
print(f"  TYPE predicted by DETECTION covariates only: AUC = {auc_det.mean():.3f} "
      f"+/- {auc_det.std():.3f}")
print("    (AUC ~0.5 => detector does not explain type; AUC high => confound)")

# gradient position (PC1 of standardized features) vs detection confidence
from sklearn.decomposition import PCA
Xs = StandardScaler().fit_transform(d[FEATS].fillna(d[FEATS].median()).to_numpy())
pc1 = PCA(2, random_state=SEED).fit_transform(Xs)[:, 0]
rho, prho = spearmanr(pc1, d.conf_mean.fillna(d.conf_mean.median()))
print(f"  gradient position (PC1) vs detection confidence: Spearman rho={rho:+.3f} (p={prho:.1e})")

# ---------------- R3-lite: false-negative perturbation ----------------------
print("\n=== R3-lite: false-negative perturbation (drop p% detections/cluster) ===")
orig_type = dict(zip(typed.cluster_id, typed.morph_type))

def refit_from_points(tdf):
    shp = compute_all_shape_metrics(tdf, cluster_col="cluster_id").reset_index()
    area = dict(zip(shp.cluster_id, shp.area_m2))
    ce = {cid: clark_evans_index(g[["x", "y"]].values, area[cid])
          for cid, g in tdf.groupby("cluster_id") if cid in area}
    shp["clark_evans"] = shp.cluster_id.map(ce)
    shp = shp[(shp.n_points >= 5) & (shp.area_m2 >= 100) &
              (shp.aspect_ratio <= 10)].dropna(subset=FEATS)
    Z = StandardScaler().fit_transform(shp[FEATS].to_numpy(float))
    km = KMeans(2, n_init=10, random_state=SEED).fit(Z)
    lab = km.labels_
    t0 = 0 if shp.aspect_ratio[lab == 0].mean() >= shp.aspect_ratio[lab == 1].mean() else 1
    return dict(zip(shp.cluster_id, np.where(lab == t0, "Type0", "Type1")))

for frac in (0.10, 0.20):
    rng = np.random.default_rng(SEED)
    keep = trees.groupby("cluster_id", group_keys=False).apply(
        lambda g: g.sample(frac=1 - frac, random_state=SEED) if len(g) > 5 else g)
    new_type = refit_from_points(keep[["x", "y", "cluster_id"]])
    shared = [c for c in new_type if c in orig_type]
    a = [orig_type[c] for c in shared]
    b = [new_type[c] for c in shared]
    flip = np.mean([x != y for x, y in zip(a, b)])
    ari = adjusted_rand_score(a, b)
    print(f"  drop {int(frac*100)}%: clusters retained={len(shared):5d}  "
          f"type-flip rate={flip*100:5.1f}%  ARI={ari:.3f}")

print("\nINTERPRETATION: high detection-AUC, a strong PC1~confidence correlation, or a"
      "\nlarge type-flip rate under modest FN rates all indicate the typology cannot be"
      "\ncleanly separated from detector behaviour. Full validation still needs R1/R2"
      "\n(stratified precision/recall vs ground truth), which require external data.")
