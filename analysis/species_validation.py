"""External validation that breaks the shape-only circularity.

The shape -> k-means "type" -> agent reading is circular: the evidence for the agent is the very
shape we clustered on. Here we bring in an EXTERNAL variable that lives outside the shape space and
outside the detector -- pre-mortality (2021) MS-NFI species composition (spruce/pine/birch volume) --
and ask whether cluster morphology covaries with the species mix of the stand that died.

  - If aspect ratio (the gradient axis) tracks species composition consistently across blocks, shape
    carries independent ecological signal -> external corroboration (non-circular).
  - If it does not, shape does not track species/host type -> strengthens "shape is not
    agent-diagnostic" (Contribution 2).

Inference respects the purposive 7-block design: the block is the replicate (effective n = 7), so we
report a within-block Spearman rho per block and combine by a sign test across blocks, never a
per-cluster p (which would be pseudo-replicated). A pooled estimate is shown only for reference.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "analysis")
from figstyle import apply_style, PALETTE, despine, panel_label, save
from host_pointprocess import _sample_raster

RAST = {"spruce": "data/kuusi_vmi1x_1721.tif", "pine": "data/manty_vmi1x_1721.tif",
        "birch": "data/koivu_vmi1x_1721.tif", "total": "data/tilavuus_vmi1x_1721.tif"}
OUT = "analysis/out"; os.makedirs(OUT, exist_ok=True)
SHAPE = "aspect_ratio"                 # the morphological-gradient axis (couples to compactness)
MINB = 60                              # minimum clusters for a block to enter the per-block test


def load():
    d = pd.read_csv("data/clusters_typed.csv")
    d["block"] = DBSCAN(eps=15000, min_samples=1).fit_predict(d[["cx", "cy"]].values)
    xy = d[["cx", "cy"]].values
    for k, p in RAST.items():
        d[k] = _sample_raster(p, xy)
    d = d[np.isfinite(d.total) & (d.total > 0)].copy()
    for k in ("spruce", "pine", "birch"):
        d[k + "_frac"] = np.clip(d[k].fillna(0) / d.total, 0, 1)
    d["conifer_frac"] = np.clip(d.spruce_frac + d.pine_frac, 0, 1)
    return d


def per_block_corr(d, species):
    """Within-block Spearman(shape, species_frac); combine across blocks by sign test."""
    rows = []
    for b, g in d.groupby("block"):
        if len(g) < MINB:
            continue
        rho, p = stats.spearmanr(g[SHAPE], g[species + "_frac"])
        rows.append((b, len(g), rho, p))
    R = pd.DataFrame(rows, columns=["block", "n", "rho", "p"])
    pos = int((R.rho > 0).sum()); k = len(R)
    sign_p = stats.binomtest(max(pos, k - pos), k, 0.5).pvalue   # two-sided sign test on rho direction
    pooled, _ = stats.spearmanr(d[SHAPE], d[species + "_frac"])
    return R, dict(species=species, n_blocks=k, pos=pos, mean_rho=R.rho.mean(),
                   sign_p=sign_p, pooled_rho=pooled)


def type_contrast(d):
    """Per-block difference in species fraction between the two morphological endpoints."""
    rows = []
    for sp in ("spruce", "pine", "birch"):
        diffs = []
        for b, g in d.groupby("block"):
            if len(g) < MINB:
                continue
            elong = g.loc[g.morph_type == "Type0", sp + "_frac"]   # Type0 = elongated/dispersed (23.5%)
            comp = g.loc[g.morph_type == "Type1", sp + "_frac"]    # Type1 = compact/aggregated (76.5%)
            if len(elong) >= 10 and len(comp) >= 10:
                diffs.append(elong.median() - comp.median())       # elongated - compact
        diffs = np.array(diffs)
        pos = int((diffs > 0).sum()); k = len(diffs)
        sp_p = stats.binomtest(max(pos, k - pos), k, 0.5).pvalue
        rows.append((sp, k, np.median(diffs), pos, sp_p))
    return pd.DataFrame(rows, columns=["species", "n_blocks", "median_diff_elong_minus_compact",
                                       "blocks_positive", "sign_p"])


def main():
    d = load()
    print(f"clusters with forest host: {len(d)} / 14582  (blocks: {sorted(d.block.unique())})")
    print(f"\n=== Q1. Does the morphological gradient ({SHAPE}) track species composition? ===")
    print("    within-block Spearman rho, combined across blocks by sign test (replicate = block)")
    summ, per = [], {}
    for sp in ("spruce", "pine", "birch", "conifer"):
        R, s = per_block_corr(d, sp); summ.append(s); per[sp] = R
        print(f"  {sp:8s}: mean rho={s['mean_rho']:+.3f}  ({s['pos']}/{s['n_blocks']} blocks +)  "
              f"sign-test p={s['sign_p']:.3f}   [pooled rho={s['pooled_rho']:+.3f}]")
    print(f"\n=== Q2. Do the two morphological endpoints sit in different species mixes? ===")
    print("    per-block median(elongated) - median(compact) species fraction, sign test across blocks")
    tc = type_contrast(d)
    for _, r in tc.iterrows():
        print(f"  {r.species:8s}: median diff={r.median_diff_elong_minus_compact:+.3f}  "
              f"({int(r.blocks_positive)}/{int(r.n_blocks)} +)  sign-test p={r.sign_p:.3f}")
    pd.DataFrame(summ).to_csv(f"{OUT}/species_corr.csv", index=False)
    tc.to_csv(f"{OUT}/species_typecontrast.csv", index=False)
    figure(d, per, tc)
    print(f"\nsaved {OUT}/species_corr.csv, {OUT}/species_typecontrast.csv, figures/fig_species_validation.*")


def figure(d, per, tc):
    apply_style()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    # (a) per-block Spearman rho of aspect ratio vs each species fraction
    species = ["spruce", "pine", "birch", "conifer"]
    cols = [PALETTE["tertiary"], PALETTE["secondary"], PALETTE["primary"], PALETTE["neutral"]]
    for i, (sp, c) in enumerate(zip(species, cols)):
        R = per[sp]
        ax[0].scatter(np.full(len(R), i) + np.linspace(-.15, .15, len(R)), R.rho,
                      s=24, color=c, zorder=3, edgecolor="white", linewidth=0.5)
        ax[0].plot([i - .25, i + .25], [R.rho.mean()] * 2, color=c, lw=2.2, zorder=4)
    ax[0].axhline(0, color=PALETTE["neutral"], lw=0.8, ls="--")
    ax[0].set_xticks(range(len(species))); ax[0].set_xticklabels(species)
    ax[0].set_ylabel(r"within-block Spearman $\rho$" + f"\n({SHAPE} vs species fraction)")
    ax[0].set_title("Shape vs species composition, per block")
    despine(ax[0]); panel_label(ax[0], "a")
    # (b) species mix by morphological endpoint (pooled distributions, for display)
    sps = ["spruce", "pine", "birch"]
    x = np.arange(len(sps)); w = 0.36
    m0 = [d.loc[d.morph_type == "Type0", sp + "_frac"].median() for sp in sps]   # elongated
    m1 = [d.loc[d.morph_type == "Type1", sp + "_frac"].median() for sp in sps]   # compact
    ax[1].bar(x - w/2, m0, w, label="elongated (Type0)", color=PALETTE["secondary"])
    ax[1].bar(x + w/2, m1, w, label="compact (Type1)", color=PALETTE["primary"])
    ax[1].set_xticks(x); ax[1].set_xticklabels(sps)
    ax[1].set_ylabel("median species fraction of total volume")
    ax[1].set_title("Species mix at the two morphological endpoints")
    ax[1].legend(loc="upper right")
    despine(ax[1]); panel_label(ax[1], "b")
    fig.tight_layout()
    save(fig, "figures/fig_species_validation")


if __name__ == "__main__":
    main()
