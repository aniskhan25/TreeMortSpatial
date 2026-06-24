"""
M1 (discreteness) + M2 (covariance-preserving silhouette null) for the revision.

Question: is the k=2 typology a genuine dichotomy, or a cut through a single
correlated continuum? K-means features = (fractal_dim, clark_evans, aspect_ratio).
"""
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

SEED = 42
rng = np.random.default_rng(SEED)

df = pd.read_csv("data/clusters_typed.csv")
feats = ["fractal_dim", "clark_evans", "aspect_ratio"]
X = df[feats].to_numpy(float)
Xs = StandardScaler().fit_transform(X)
n = len(Xs)
print(f"n = {n} clusters; features = {feats}\n")

# --- reproduce the reported k=2 silhouette ---------------------------------
km2 = KMeans(n_clusters=2, n_init=10, random_state=SEED).fit(Xs)
sil_obs = silhouette_score(Xs, km2.labels_)
print(f"[repro] observed k=2 silhouette = {sil_obs:.4f}  (paper: 0.401)\n")

# --- M1a: GMM BIC, 1 vs k components --------------------------------------
print("[M1a] Gaussian-mixture BIC (lower = better):")
bics = {}
for k in range(1, 8):
    gm = GaussianMixture(n_components=k, covariance_type="full",
                         n_init=5, random_state=SEED).fit(Xs)
    bics[k] = gm.bic(Xs)
    print(f"   k={k}: BIC = {bics[k]:.1f}")
best_k = min(bics, key=bics.get)
print(f"   -> BIC-optimal k = {best_k}")
print(f"   -> 2-comp better than 1-comp? {'YES' if bics[2] < bics[1] else 'NO'} "
      f"(delta BIC 1->2 = {bics[1]-bics[2]:.1f})\n")

# --- M1b: bimodality on PC1 -----------------------------------------------
pc1 = PCA(n_components=1, random_state=SEED).fit_transform(Xs).ravel()
g1 = pd.Series(pc1).skew()
k_ex = pd.Series(pc1).kurt()  # excess kurtosis (Fisher)
# Sarle's bimodality coefficient; >0.555 (uniform) suggests bimodal
BC = (g1**2 + 1) / (k_ex + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
print(f"[M1b] PC1 explained var = {PCA(n_components=3).fit(Xs).explained_variance_ratio_[0]:.3f}")
print(f"      Sarle bimodality coefficient = {BC:.3f}  (>0.555 hints bimodal)")
try:
    import diptest
    d, pdip = diptest.diptest(pc1)
    print(f"      Hartigan dip test: D={d:.4f}, p={pdip:.4g} "
          f"({'multimodal' if pdip < 0.05 else 'unimodal not rejected'})")
except Exception as e:
    print(f"      (diptest unavailable: {e})")
print()

# --- M1c: how many clusters sit near the k=2 boundary ---------------------
d_cent = km2.transform(Xs)               # distance to each centroid
d_sorted = np.sort(d_cent, axis=1)
margin = (d_sorted[:, 1] - d_sorted[:, 0]) / d_sorted[:, 1]  # 0 = on boundary
for thr in (0.05, 0.10, 0.20):
    frac = np.mean(margin < thr)
    print(f"[M1c] clusters within {int(thr*100)}% of the k=2 boundary: "
          f"{frac*100:.1f}%")
print()

# --- M2: covariance-preserving silhouette null ----------------------------
mu = Xs.mean(axis=0)
cov = np.cov(Xs, rowvar=False)
B = 999
null = np.empty(B)
for b in range(B):
    Z = rng.multivariate_normal(mu, cov, size=n)
    lab = KMeans(n_clusters=2, n_init=5, random_state=SEED).fit_predict(Z)
    null[b] = silhouette_score(Z, lab)
p = (np.sum(null >= sil_obs) + 1) / (B + 1)
print("[M2] silhouette vs SINGLE multivariate-Gaussian null (preserves covariance):")
print(f"     observed = {sil_obs:.4f}")
print(f"     null mean = {null.mean():.4f}, 95th pct = {np.percentile(null,95):.4f}, "
      f"max = {null.max():.4f}")
print(f"     p(observed >= null) = {p:.4g}  "
      f"-> {'k=2 EXCEEDS unimodal null' if p < 0.05 else 'k=2 NOT distinguishable from a single Gaussian'}")
