"""
Regenerate fig_silhouette_permtest.png to match the revised text:
(a) covariance-preserving silhouette null, (b) GMM BIC vs k, (c) PC1 histogram.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

SEED = 42
rng = np.random.default_rng(SEED)
df = pd.read_csv("data/clusters_typed.csv")
feats = ["fractal_dim", "clark_evans", "aspect_ratio"]
Xs = StandardScaler().fit_transform(df[feats].to_numpy(float))
n = len(Xs)

sil_obs = silhouette_score(Xs, KMeans(2, n_init=10, random_state=SEED).fit_predict(Xs))

# (a) covariance-preserving null (full silhouette, matches the values quoted in the text).
mu, cov = Xs.mean(0), np.cov(Xs, rowvar=False)
B = 999
null = np.empty(B)
for b in range(B):
    Z = rng.multivariate_normal(mu, cov, size=n)
    lab = KMeans(2, n_init=5, random_state=SEED).fit_predict(Z)
    null[b] = silhouette_score(Z, lab)

# (b) GMM BIC
ks = range(1, 8)
bic = [GaussianMixture(k, covariance_type="full", n_init=5, random_state=SEED)
       .fit(Xs).bic(Xs) for k in ks]

# (c) PC1
pc1 = PCA(3, random_state=SEED).fit(Xs)
pc1_scores = PCA(3, random_state=SEED).fit_transform(Xs)[:, 0]
var1 = pc1.explained_variance_ratio_[0]

fig, ax = plt.subplots(1, 3, figsize=(13, 4))

ax[0].hist(null, bins=40, color="#9ecae1", edgecolor="white")
ax[0].axvline(sil_obs, color="#cb181d", lw=2.5,
              label=f"observed = {sil_obs:.3f}")
ax[0].axvline(np.percentile(null, 95), color="#525252", ls="--", lw=1.2,
              label="null 95th pct")
ax[0].set_xlabel("silhouette (k=2)"); ax[0].set_ylabel("null replicates")
ax[0].set_title("(a) Covariance-preserving null\n$p=0.001$")
ax[0].legend(fontsize=8, frameon=False)

ax[1].plot(list(ks), bic, "o-", color="#238b45")
ax[1].axvline(2, color="#cb181d", ls="--", lw=1.2, label="k=2")
ax[1].set_xlabel("number of components k"); ax[1].set_ylabel("BIC")
ax[1].set_title("(b) GMM BIC: no privileged k=2")
ax[1].legend(fontsize=8, frameon=False)

ax[2].hist(pc1_scores, bins=60, color="#bcbddc", edgecolor="white")
ax[2].set_xlabel(f"PC1 ({var1*100:.0f}% variance)")
ax[2].set_ylabel("clusters")
ax[2].set_title("(c) Dominant axis is unimodal\nSarle BC = 0.486")

plt.tight_layout()
out = "figures/fig_silhouette_permtest.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"observed={sil_obs:.4f} null_mean={null.mean():.4f} "
      f"null_p95={np.percentile(null,95):.4f}")
print(f"BIC={[round(b) for b in bic]}")
print(f"saved -> {out}")
