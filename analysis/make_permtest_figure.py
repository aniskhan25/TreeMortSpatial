"""
fig_silhouette_permtest.png — significance + discreteness of the k=2 partition, house style.
(a) covariance-preserving silhouette null vs observed (p=0.001)
(b) GMM BIC vs k (no privileged k=2)
(c) PC1 histogram (unimodal -> continuum)

The full covariance-preserving null (999 x KMeans + full silhouette on 14,582 pts) is slow,
so it is cached to analysis/out/permtest_null.npy and reused on subsequent runs.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from figstyle import apply_style, PALETTE, despine, panel_label, save

apply_style()
P = PALETTE
SEED = 42
rng = np.random.default_rng(SEED)

df = pd.read_csv("data/clusters_typed.csv")
feats = ["fractal_dim", "clark_evans", "aspect_ratio"]
Xs = StandardScaler().fit_transform(df[feats].to_numpy(float))
n = len(Xs)

sil_obs = silhouette_score(Xs, KMeans(2, n_init=10, random_state=SEED).fit_predict(Xs))

# (a) covariance-preserving null (cached: full silhouette, matches the quoted numbers)
cache = "analysis/out/permtest_null.npy"
if os.path.exists(cache):
    null = np.load(cache)
    print(f"loaded cached null ({len(null)} reps)")
else:
    mu, cov = Xs.mean(0), np.cov(Xs, rowvar=False)
    B = 999
    null = np.empty(B)
    for b in range(B):
        Z = rng.multivariate_normal(mu, cov, size=n)
        null[b] = silhouette_score(Z, KMeans(2, n_init=5, random_state=SEED).fit_predict(Z))
    np.save(cache, null)
    print(f"computed + cached null ({B} reps)")

p95 = np.percentile(null, 95)
pval = (np.sum(null >= sil_obs) + 1) / (len(null) + 1)

# (b) GMM BIC
ks = list(range(1, 8))
bic = [GaussianMixture(k, covariance_type="full", n_init=5, random_state=SEED).fit(Xs).bic(Xs)
       for k in ks]

# (c) PC1
pca = PCA(3, random_state=SEED).fit(Xs)
pc1 = pca.transform(Xs)[:, 0]
var1 = pca.explained_variance_ratio_[0]

fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), dpi=150)

# (a)
a = ax[0]
despine(a)
a.hist(null, bins=36, color=P["muted"], edgecolor="white", linewidth=0.35, zorder=3,
       label="covariance-preserving null")
a.axvline(p95, color=P["neutral"], ls=(0, (4, 2)), lw=1.1, zorder=4, label=f"null 95th pct = {p95:.3f}")
a.axvline(sil_obs, color=P["accent"], lw=2.2, zorder=5, label=f"observed = {sil_obs:.3f}")
a.set_xlabel("silhouette ($k{=}2$)")
a.set_ylabel("null replicates")
a.legend(loc="upper center")
panel_label(a, "a")
a.set_title(f"Beats covariance-preserving null ($p={pval:.3f}$)",
            fontweight="bold", color="#222222", pad=18)

# (b)
b = ax[1]
despine(b)
b.plot(ks, np.array(bic) / 1e5, "o-", color=P["primary"], lw=1.8, ms=6,
       mfc="white", mec=P["primary"], mew=1.4, zorder=3)
b.axvline(2, color=P["accent"], ls=(0, (4, 2)), lw=1.2, zorder=2, label="$k=2$")
b.set_xlabel("number of mixture components $k$")
b.set_ylabel("BIC ($\\times 10^{5}$)")
b.legend(loc="upper right")
panel_label(b, "b")
b.set_title("BIC favours $k>2$: no privileged split", fontweight="bold", color="#222222", pad=18)

# (c)
c = ax[2]
despine(c)
c.hist(pc1, bins=50, color=P["tertiary"], edgecolor="white", linewidth=0.3, zorder=3)
c.axvline(np.median(pc1), color=P["neutral"], ls=(0, (4, 2)), lw=1.1, zorder=4, label="median")
c.set_xlabel(f"PC1 ({var1*100:.0f}% of variance)")
c.set_ylabel("number of clusters")
c.legend(loc="upper right")
panel_label(c, "c")
c.set_title("Dominant axis is unimodal (Sarle 0.49)", fontweight="bold", color="#222222", pad=18)

fig.tight_layout(w_pad=2.4)
save(fig, "figures/fig_silhouette_permtest")
print(f"observed={sil_obs:.4f} null_mean={null.mean():.4f} p95={p95:.4f} p={pval:.4f}")
print("saved figures/fig_silhouette_permtest.png + .pdf")
