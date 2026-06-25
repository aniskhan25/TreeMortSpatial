"""
Host-conditioned mortality point-process engine (Ecology flagship, Phase 1).

Core idea (Vasquez/Lindqvist/Anand): test whether mortality clusters MORE than the forest
host already explains, by replacing the meaningless homogeneous-CSR null (which gave the
retired E4/E7 artefact g(r)=400-525 over heterogeneous forest) with an INHOMOGENEOUS null
whose first-order intensity is driven by host availability.

Estimators (computed on cluster centroids, per super-site block, the 7 PSUs):
  - K_inhom(r): inhomogeneous Ripley's K with translation edge correction,
        K_inhom(r) = (1/|W|) * sum_{i!=j, d<=r} w_ij / (lam_i lam_j);  E[K]=pi r^2 under inhom. Poisson.
  - g_inhom(r): ring derivative of K_inhom; g=1 under inhom. Poisson, >1 aggregation, <1 inhibition.

Intensity lam(u):
  - Phase 2: sampled from the MS-NFI host raster (pass --raster).
  - Phase 1 (no raster): a kernel forest-presence proxy, used ONLY to exercise the code path;
    NOT an ecological substitute for real host density.

Run `python analysis/host_pointprocess.py --selftest` for the simulation unit-test that proves
g_inhom~1 on an inhomogeneous Poisson pattern while the CSR-PCF spuriously reports aggregation.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde


# --------------------------------------------------------------------------- #
# Inhomogeneous second-order estimators (rectangular window, translation edge) #
# --------------------------------------------------------------------------- #
def _translation_weights(coords, W):
    """Pairwise translation edge-correction weights for a rectangular window W=(x0,x1,y0,y1).
    w_ij = |W| / area(W intersect W shifted by (x_j-x_i)).  Returns (n,n) matrix."""
    x0, x1, y0, y1 = W
    Lx, Ly = x1 - x0, y1 - y0
    dx = np.abs(coords[:, 0][:, None] - coords[:, 0][None, :])
    dy = np.abs(coords[:, 1][:, None] - coords[:, 1][None, :])
    overlap = np.clip(Lx - dx, 1e-9, None) * np.clip(Ly - dy, 1e-9, None)
    return (Lx * Ly) / overlap


def k_num(coords, lam, r_vals, W):
    """Numerator of the inhomogeneous K on rectangular window W (NOT divided by area):
    N(r) = sum_{i!=j, d<=r} w_ij / (lam_i lam_j). Used for pooling across tiles."""
    if len(coords) < 2:
        return np.zeros(len(r_vals))
    d = np.hypot(coords[:, 0][:, None] - coords[:, 0][None, :],
                 coords[:, 1][:, None] - coords[:, 1][None, :])
    w = _translation_weights(coords, W)
    inv = 1.0 / np.outer(lam, lam)
    np.fill_diagonal(d, np.inf)
    contrib = w * inv
    return np.array([contrib[d <= r].sum() for r in r_vals])


def k_inhom(coords, lam, r_vals, W):
    """Inhomogeneous K on a single rectangular window = numerator / area."""
    x0, x1, y0, y1 = W
    return k_num(coords, lam, r_vals, W) / ((x1 - x0) * (y1 - y0))


def k_pooled(tiles, r_vals):
    """Pool the inhomogeneous K across independent rectangular tile-windows.
    tiles: list of dicts {coords, lam, W}.  K(r) = sum_t num_t(r) / sum_t area_t.
    Under inhomogeneous Poisson (lam = true intensity) E[K]=pi r^2, so g=1."""
    num = np.zeros(len(r_vals)); denom = 0.0
    for t in tiles:
        x0, x1, y0, y1 = t["W"]
        num += k_num(t["coords"], t["lam"], r_vals, t["W"])
        denom += (x1 - x0) * (y1 - y0)
    return num / denom


def g_from_k(r_vals, K):
    """Pair-correlation g(r) from K via ring derivative: g=ΔK/(π Δ(r²)). g=1 ~ Poisson."""
    r_mid = 0.5 * (r_vals[:-1] + r_vals[1:])
    dK = np.diff(K)
    dr2 = np.pi * np.diff(r_vals ** 2)
    return r_mid, dK / dr2


def g_csr(coords, r_vals, W):
    """Homogeneous (CSR) PCF: same ring estimator but constant intensity = n/|W|.
    This is the WRONG null over heterogeneous forest; included to show the artefact."""
    x0, x1, y0, y1 = W
    area = (x1 - x0) * (y1 - y0)
    lam0 = len(coords) / area
    return g_from_k(r_vals, k_inhom(coords, np.full(len(coords), lam0), r_vals, W))


# --------------------------------------------------------------------------- #
# Intensity models                                                            #
# --------------------------------------------------------------------------- #
def lam_from_raster(coords, raster_path):
    """Phase 2: sample host (e.g. spruce volume, m3/ha) at coords from an MS-NFI GeoTIFF.
    coords are EPSG:3067; reprojects to the raster CRS if it differs. Returns a
    strictly-positive intensity proportional to host density (floored where host~0)."""
    import rasterio
    with rasterio.open(raster_path) as src:
        pts = coords
        if src.crs and src.crs.to_epsg() not in (3067, None):
            from pyproj import Transformer
            tr = Transformer.from_crs(3067, src.crs.to_epsg(), always_xy=True)
            xs, ys = tr.transform(coords[:, 0], coords[:, 1])
            pts = np.column_stack([xs, ys])
        nod = src.nodata
        vals = np.array([v[0] for v in src.sample(pts)], dtype=float)
    vals = np.where(np.isfinite(vals), vals, 0.0)
    if nod is not None:
        vals = np.where(vals == nod, 0.0, vals)
    # floor so non-host (volume~0) locations get a small positive intensity, not div-by-zero
    floor = np.nanpercentile(vals[vals > 0], 5) if (vals > 0).any() else 1.0
    return np.clip(vals + floor, 1e-9, None)


def lam_proxy(coords, bw=None):
    """Phase 1 placeholder: kernel forest-presence proxy from the points themselves.
    NOT a substitute for host density (it is partly circular); only exercises the code path."""
    kde = gaussian_kde(coords.T, bw_method=bw)
    lam = kde(coords.T)
    return np.clip(lam, 1e-12, None)


def _scale_intensity(lam, n, W):
    """Scale a relative intensity so it integrates to n over W (mean-matches point count)."""
    x0, x1, y0, y1 = W
    area = (x1 - x0) * (y1 - y0)
    return lam * (n / (lam.mean() * area)) * (area / area)  # mean(lam)*area == n


# --------------------------------------------------------------------------- #
# Simulation self-test                                                        #
# --------------------------------------------------------------------------- #
def selftest(seed=42):
    """Inhomogeneous Poisson with a known east-west intensity gradient.
    Expect: g_inhom(r) ~ 1 (no true clustering); g_csr(r) > 1 (spurious, from the gradient)."""
    rng = np.random.default_rng(seed)
    W = (0.0, 1000.0, 0.0, 1000.0)
    x0, x1, y0, y1 = W
    # true intensity: exp gradient in x, ~2000 expected points
    def lam_true_fn(xy):
        return np.exp(2.0 * (xy[:, 0] - x0) / (x1 - x0))   # relative
    # sample by thinning a homogeneous Poisson
    n_prop = 9000
    pts = np.column_stack([rng.uniform(x0, x1, n_prop), rng.uniform(y0, y1, n_prop)])
    lam_rel = lam_true_fn(pts)
    keep = rng.random(n_prop) < lam_rel / lam_rel.max()
    pts = pts[keep]
    n = len(pts)
    # intensity at points, scaled to integrate to n
    lam_pts = lam_true_fn(pts)
    lam_pts = lam_pts * n / (np.mean(lam_true_fn(
        np.column_stack([rng.uniform(x0, x1, 20000), rng.uniform(y0, y1, 20000)]))) *
        (x1 - x0) * (y1 - y0))
    r_vals = np.linspace(0, 150, 16)
    Kin = k_inhom(pts, lam_pts, r_vals, W)
    rm, gin = g_from_k(r_vals, Kin)
    _, gcsr = g_csr(pts, r_vals, W)
    print(f"[selftest] inhomogeneous Poisson, n={n}")
    print(f"  r(m):        {np.round(rm[1:],0)}")
    print(f"  g_inhom(r):  {np.round(gin[1:],2)}   (expect ~1.0 — host-conditioned null is correct)")
    print(f"  g_csr(r):    {np.round(gcsr[1:],2)}   (expect >1 — CSR mistakes the gradient for clustering)")
    ok = np.nanmedian(np.abs(gin[1:] - 1)) < 0.25 and np.nanmedian(gcsr[1:]) > np.nanmedian(gin[1:])
    print(f"  RESULT: {'PASS' if ok else 'CHECK'} — "
          f"median|g_inhom-1|={np.nanmedian(np.abs(gin[1:]-1)):.2f}, "
          f"median g_csr={np.nanmedian(gcsr[1:]):.2f} vs g_inhom={np.nanmedian(gin[1:]):.2f}")
    return ok


# --------------------------------------------------------------------------- #
# Per-block real-data driver                                                  #
# --------------------------------------------------------------------------- #
def _window_mean_host(raster, W, ngrid=150):
    """Mean host (spruce volume) over rectangular window W, from a grid sample of the raster.
    Needed to normalise the inhomogeneous intensity correctly (g scales with this constant)."""
    x0, x1, y0, y1 = W
    gx = np.linspace(x0, x1, ngrid); gy = np.linspace(y0, y1, ngrid)
    XX, YY = np.meshgrid(gx, gy)
    h = lam_from_raster(np.column_stack([XX.ravel(), YY.ravel()]), raster)
    return float(np.mean(h))


def run_blocks(raster=None, r_max=2000.0, nr=20, save="analysis/out/host_pcf.json"):
    import json, os
    from scipy.stats import spearmanr
    from sklearn.cluster import DBSCAN
    d = pd.read_csv("data/clusters_typed.csv")
    d["block"] = DBSCAN(eps=15000, min_samples=1).fit_predict(d[["cx", "cy"]].values)
    lat = lambda y: 60 + (y - 6660000) / 111320.0
    r_vals = np.linspace(0, r_max, nr + 1)
    mode = "MS-NFI spruce volume" if raster else "forest-presence PROXY (placeholder, not host)"
    print(f"\n=== per-block host-conditioned inhomogeneous PCF  [intensity: {mode}] ===")
    out = {}
    rows = []
    for b, g in d.groupby("block"):
        xy = g[["cx", "cy"]].to_numpy(float)
        if len(xy) < 50:
            continue
        W = (xy[:, 0].min(), xy[:, 0].max(), xy[:, 1].min(), xy[:, 1].max())
        area = (W[1] - W[0]) * (W[3] - W[2])
        n = len(xy)
        if raster:
            host_pts = lam_from_raster(xy, raster)
            mean_host = _window_mean_host(raster, W)           # window integral of host
            lam = host_pts * n / (mean_host * area)            # c = n / ∫host ; λ_i = c·host_i
            host_frac = float(np.mean(host_pts > np.nanpercentile(host_pts, 5) + 1))
        else:
            lam = lam_proxy(xy); lam = lam * n / (lam.mean() * area); host_frac = np.nan
        rm, gin = g_from_k(r_vals, k_inhom(xy, lam, r_vals, W))
        _, gcsr = g_from_k(r_vals, k_inhom(xy, np.full(n, n / area), r_vals, W))
        rec = dict(block=int(b), lat=round(lat(g.cy.mean()), 1), n=n,
                   r=rm.tolist(), g_inhom=gin.tolist(), g_csr=gcsr.tolist(),
                   g_inhom_200=round(float(np.interp(200, rm, gin)), 2),
                   g_inhom_500=round(float(np.interp(500, rm, gin)), 2),
                   g_csr_500=round(float(np.interp(500, rm, gcsr)), 2),
                   mean_spruce=round(float(np.nanmean(host_pts)), 0) if raster else None)
        out[str(b)] = rec
        rows.append(rec)
    print(f"{'blk':>3} {'lat':>5} {'n':>5} {'spruce':>7} {'g_in(200m)':>10} {'g_in(500m)':>10} {'g_csr(500m)':>11}")
    for r in sorted(rows, key=lambda z: z["lat"]):
        print(f"{r['block']:>3} {r['lat']:>5.1f} {r['n']:>5} {str(r['mean_spruce']):>7} "
              f"{r['g_inhom_200']:>10} {r['g_inhom_500']:>10} {r['g_csr_500']:>11}")
    # latitudinal trend in host-conditioned aggregation (7 blocks; honest small-n)
    la = [r["lat"] for r in rows]; gi500 = [r["g_inhom_500"] for r in rows]
    if len(rows) >= 4:
        rho, p = spearmanr(la, gi500)
        print(f"\nlatitudinal trend in g_inhom(500m): Spearman rho={rho:+.2f}, p={p:.3f} (n={len(rows)} blocks)")
        print("(rho>0 => host-conditioned aggregation strengthens northward — the leading-edge signal)")
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        json.dump(out, open(save, "w"), indent=2)
        print(f"saved {save}")
    print("\ng_inhom > 1 = mortality clusters MORE than spruce host explains (genuine aggregation);"
          "\nthe g_csr column is the discredited homogeneous null, inflated by forest heterogeneity.")
    return out


# --------------------------------------------------------------------------- #
# Pooled, tile-windowed host-conditioned analysis (the correct observation window) #
# --------------------------------------------------------------------------- #
def build_tiles(raster, min_clusters=10):
    """Map filtered clusters -> 6 km aerial tiles (the true observation windows), assign blocks,
    and attach host (spruce) intensity. Returns list of tile dicts and a global rate rho."""
    from sklearn.cluster import DBSCAN
    typed = pd.read_csv("data/clusters_typed.csv")
    tr = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id", "image_name"])
    tr = tr[tr.cluster_id != -1]
    # tile window = bbox of a tile's detections (~6 km imaged rectangle)
    tw = tr.groupby("image_name").agg(x0=("x", "min"), x1=("x", "max"),
                                      y0=("y", "min"), y1=("y", "max"))
    cl_tile = tr.groupby("cluster_id").image_name.agg(lambda s: s.mode().iloc[0])
    typed = typed[typed.cluster_id.isin(cl_tile.index)].copy()
    typed["tile"] = typed.cluster_id.map(cl_tile)
    # blocks from tile-window centres
    tw["cx"] = (tw.x0 + tw.x1) / 2; tw["cy"] = (tw.y0 + tw.y1) / 2
    tw["block"] = DBSCAN(eps=15000, min_samples=1).fit_predict(tw[["cx", "cy"]].values)
    # build per-tile records; global host integral for rho
    tiles, host_integral, N = [], 0.0, 0
    for tile, sub in typed.groupby("tile"):
        if len(sub) < min_clusters or tile not in tw.index:
            continue
        W = (tw.loc[tile, "x0"], tw.loc[tile, "x1"], tw.loc[tile, "y0"], tw.loc[tile, "y1"])
        area = (W[1] - W[0]) * (W[3] - W[2])
        if area <= 0:
            continue
        xy = sub[["cx", "cy"]].to_numpy(float)
        host_pts = lam_from_raster(xy, raster)
        mean_host = _window_mean_host(raster, W, ngrid=40)
        tiles.append(dict(tile=tile, block=int(tw.loc[tile, "block"]), W=W, area=area,
                          coords=xy, host=host_pts, n=len(xy),
                          cy=sub.cy.mean(), mean_host=mean_host))
        host_integral += mean_host * area
        N += len(xy)
    for t in tiles:
        # PER-TILE inhomogeneous null: lambda(u) = spruce(u) * n_t / integral_t(spruce),
        # so it integrates to n_t within the tile. This frees each tile's overall rate, so the
        # pooled g measures WITHIN-tile aggregation beyond host structure (not between-tile rate
        # variation, which is the separate first-order leading-edge signal reported elsewhere).
        t["lam"] = t["host"] * t["n"] / (t["mean_host"] * t["area"])
        t["lam0"] = np.full(t["n"], t["n"] / t["area"])   # per-tile homogeneous (CSR) null
    rho_global = N / host_integral               # first-order mortality per unit spruce (for trend)
    return tiles, rho_global


def run_pooled(raster, r_max=1500.0, nr=15, save="analysis/out/host_pcf_pooled.json"):
    import json, os
    from scipy.stats import spearmanr
    lat = lambda y: 60 + (y - 6660000) / 111320.0
    r_vals = np.linspace(0, r_max, nr + 1)
    tiles, rho = build_tiles(raster)
    print(f"\n=== POOLED host-conditioned inhomogeneous PCF (tile windows) ===")
    print(f"tiles used: {len(tiles)}  | clusters: {sum(t['n'] for t in tiles):,}  | "
          f"global rate rho={rho:.3e} clusters per (m3/ha * m2)")

    def g_of(tile_subset, null="host"):
        key = "lam" if null == "host" else "lam0"
        K = k_pooled([{"coords": t["coords"], "lam": t[key], "W": t["W"]} for t in tile_subset], r_vals)
        return g_from_k(r_vals, K)

    rm, g_all = g_of(tiles, "host")
    _, g_all_csr = g_of(tiles, "csr")
    print(f"\nALL tiles pooled:")
    print(f"  r(m):        {np.round(rm[::2],0)}")
    print(f"  g_inhom(r):  {np.round(g_all[::2],2)}   (>1 = clusters beyond spruce host)")
    print(f"  g_csr(r):    {np.round(g_all_csr[::2],2)}   (homogeneous null, for contrast)")

    print(f"\nper block (tiles grouped):")
    print(f"{'blk':>3} {'lat':>5} {'tiles':>5} {'clus':>5} {'g_in(200)':>9} {'g_in(500)':>9} {'g_csr(500)':>10}")
    out = {"all": {"r": rm.tolist(), "g_inhom": g_all.tolist(), "g_csr": g_all_csr.tolist()}, "blocks": {}}
    rows = []
    for b in sorted({t["block"] for t in tiles}):
        ts = [t for t in tiles if t["block"] == b]
        if sum(t["n"] for t in ts) < 100:
            continue
        rmb, gb = g_of(ts, "host"); _, gbc = g_of(ts, "csr")
        la = lat(np.mean([t["cy"] for t in ts]))
        g2, g5, gc5 = (float(np.interp(x, rmb, arr)) for x, arr in [(200, gb), (500, gb), (500, gbc)])
        print(f"{b:>3} {la:>5.1f} {len(ts):>5} {sum(t['n'] for t in ts):>5} "
              f"{g2:>9.2f} {g5:>9.2f} {gc5:>10.2f}")
        out["blocks"][str(b)] = dict(lat=round(la, 1), tiles=len(ts), n=sum(t["n"] for t in ts),
                                     g_inhom_200=round(g2, 2), g_inhom_500=round(g5, 2),
                                     g_csr_500=round(gc5, 2), r=rmb.tolist(), g_inhom=gb.tolist())
        rows.append((la, g5))
    if len(rows) >= 4:
        la_, g5_ = zip(*rows)
        rho_s, p_s = spearmanr(la_, g5_)
        print(f"\n[second-order] latitudinal trend in g_inhom(500m): Spearman rho={rho_s:+.2f}, "
              f"p={p_s:.3f} (n={len(rows)} blocks; rho>0 => within-stand contagion strengthens northward)")
        out["lat_trend_g"] = {"rho": round(float(rho_s), 3), "p": round(float(p_s), 3), "n": len(rows)}

    # first-order LEADING-EDGE signal: mortality per unit spruce, by block vs latitude
    print(f"\n[first-order] mortality per unit spruce (clusters per m3/ha*km2), by block:")
    print(f"{'blk':>3} {'lat':>5} {'mort/spruce':>12}")
    fo = []
    for b in sorted({t["block"] for t in tiles}):
        ts = [t for t in tiles if t["block"] == b]
        n_b = sum(t["n"] for t in ts)
        host_int = sum(t["mean_host"] * t["area"] for t in ts)
        if host_int <= 0:
            continue
        rate = n_b / host_int * 1e6              # per m3/ha per km2
        la = lat(np.mean([t["cy"] for t in ts]))
        fo.append((la, rate))
        print(f"{b:>3} {la:>5.1f} {rate:>12.3f}")
    if len(fo) >= 4:
        la2, rt2 = zip(*fo)
        rs, ps = spearmanr(la2, rt2)
        print(f"latitudinal trend in mortality/spruce: Spearman rho={rs:+.2f}, p={ps:.3f} "
              f"(rho>0 => host-relative mortality rises northward = leading-edge signal)")
        out["lat_trend_rate"] = {"rho": round(float(rs), 3), "p": round(float(ps), 3), "n": len(fo)}
    os.makedirs(os.path.dirname(save), exist_ok=True)
    json.dump(out, open(save, "w"), indent=2)
    print(f"saved {save}")
    return out


def selftest_pooled(seed=1):
    """Pooled estimator on many small tiles, each an inhomogeneous Poisson with a host gradient.
    Expect pooled g_inhom~1 (correct) and g_csr>1 (homogeneous null mistakes the gradient)."""
    rng = np.random.default_rng(seed)
    tiles = []
    L = 6000.0
    Ntiles = 60
    for _ in range(Ntiles):
        ox, oy = rng.uniform(0, 5e5), rng.uniform(0, 5e5)
        W = (ox, ox + L, oy, oy + L)
        host_fn = lambda xy, a=rng.uniform(0.5, 2.5): np.exp(a * (xy[:, 0] - ox) / L)
        prop = np.column_stack([rng.uniform(ox, ox + L, 1500), rng.uniform(oy, oy + L, 1500)])
        hr = host_fn(prop); keep = rng.random(len(prop)) < hr / hr.max()
        xy = prop[keep]
        if len(xy) < 10:
            continue
        host = host_fn(xy)
        tiles.append(dict(coords=xy, W=W, host=host, area=L * L, n=len(xy),
                          mean_host=host_fn(np.column_stack([rng.uniform(ox, ox+L, 4000),
                                                             rng.uniform(oy, oy+L, 4000)])).mean()))
    for t in tiles:  # per-tile normalisation (matches build_tiles)
        t["lam"] = t["host"] * t["n"] / (t["mean_host"] * t["area"])
        t["lam0"] = np.full(t["n"], t["n"] / t["area"])
    N = sum(t["n"] for t in tiles)
    r_vals = np.linspace(0, 800, 17)
    rm, gin = g_from_k(r_vals, k_pooled([{"coords": t["coords"], "lam": t["lam"], "W": t["W"]} for t in tiles], r_vals))
    _, gcsr = g_from_k(r_vals, k_pooled([{"coords": t["coords"], "lam": t["lam0"], "W": t["W"]} for t in tiles], r_vals))
    print(f"[selftest_pooled] {len(tiles)} tiles, N={N}")
    print(f"  g_inhom: {np.round(gin[1:],2)}  (expect ~1)")
    print(f"  g_csr:   {np.round(gcsr[1:],2)}  (expect >1)")
    ok = np.nanmedian(np.abs(gin[1:] - 1)) < 0.25 and np.nanmedian(gcsr[1:]) > 1.1
    print(f"  RESULT: {'PASS' if ok else 'CHECK'} (median|g_inhom-1|={np.nanmedian(np.abs(gin[1:]-1)):.2f}, "
          f"median g_csr={np.nanmedian(gcsr[1:]):.2f})")
    return ok


# --------------------------------------------------------------------------- #
# Forest-masked, pre-mortality-host run with host-conditioned envelopes        #
# (the two decisive fixes: land-class forest window + 2021 spruce host)        #
# --------------------------------------------------------------------------- #
def _sample_raster(path, xy):
    import rasterio
    with rasterio.open(path) as s:
        nod = s.nodata
        v = np.array([t[0] for t in s.sample(xy)], dtype=float)
    return np.where((v == nod) if nod is not None else False, np.nan, v)


def run_masked(spruce_path, landclass_path, forest_classes=(1,),
               r_max=1500.0, nr=15, nsim=99, grid=60, save="analysis/out/host_pcf_masked.json"):
    """Tree/centroid inhomogeneous PCF on FOREST-MASKED tile windows with PRE-mortality spruce host
    and host-conditioned simulation ENVELOPES (the spruce-proportional inhomogeneous null)."""
    import json, os
    from sklearn.cluster import DBSCAN
    rng = np.random.default_rng(42)
    d = pd.read_csv("data/clusters_typed.csv")
    tr = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id", "image_name"])
    tr = tr[tr.cluster_id != -1]
    cl_tile = tr.groupby("cluster_id").image_name.agg(lambda s: s.mode().iloc[0])
    tw = tr.groupby("image_name").agg(x0=("x", "min"), x1=("x", "max"), y0=("y", "min"), y1=("y", "max"))
    d = d[d.cluster_id.isin(cl_tile.index)].copy(); d["tile"] = d.cluster_id.map(cl_tile)
    r_vals = np.linspace(0, r_max, nr + 1)

    obs_tiles, sim_pools = [], [[] for _ in range(nsim)]
    used = 0
    for tile, sub in d.groupby("tile"):
        if len(sub) < 15 or tile not in tw.index:
            continue
        W = (tw.loc[tile, "x0"], tw.loc[tile, "x1"], tw.loc[tile, "y0"], tw.loc[tile, "y1"])
        # forest-masked grid over the tile: cells that are forest land, with their spruce value
        gx = np.linspace(W[0], W[1], grid); gy = np.linspace(W[2], W[3], grid)
        GX, GY = np.meshgrid(gx, gy); cells = np.column_stack([GX.ravel(), GY.ravel()])
        lc = _sample_raster(landclass_path, cells)
        forest = np.isin(np.nan_to_num(lc, nan=-1).astype(int), forest_classes)
        if forest.sum() < 25:
            continue
        fcells = cells[forest]
        spr_cells = np.nan_to_num(_sample_raster(spruce_path, fcells), nan=0.0)
        if spr_cells.sum() <= 0:
            continue
        cell_area = ((W[1]-W[0])/(grid-1)) * ((W[3]-W[2])/(grid-1))
        forest_area = forest.sum() * cell_area
        # host = spruce + a robust floor (10th pct of forest spruce) so 1/lambda is bounded
        floor = np.percentile(spr_cells[spr_cells > 0], 10) if (spr_cells > 0).any() else 1.0
        host_cells = spr_cells + floor
        mean_host = host_cells.mean()
        norm = len(sub) / (mean_host * forest_area)             # per-tile, integrates to n over forest
        n = len(sub)
        lam_floor = 0.15 * n / forest_area                       # symmetric lower bound on lambda
        # observed points: host sampled at their locations
        xy = sub[["cx", "cy"]].to_numpy(float)
        spr_pts = np.nan_to_num(_sample_raster(spruce_path, xy), nan=0.0) + floor
        lam = np.maximum(spr_pts * norm, lam_floor)
        obs_tiles.append({"coords": xy, "lam": lam, "W": W, "area": forest_area})
        # host-conditioned sims: n points placed on forest cells with prob ∝ host; lambda = the
        # CHOSEN cell's host (consistent with placement, no re-sample -> fast). Same floor as obs.
        p = host_cells / host_cells.sum()
        cellsize = np.sqrt(cell_area)
        for s in range(nsim):
            idx = rng.choice(len(fcells), size=n, p=p)
            sxy = fcells[idx] + (rng.random((n, 2)) - 0.5) * cellsize
            sim_pools[s].append({"coords": sxy, "lam": np.maximum(host_cells[idx] * norm, lam_floor),
                                 "W": W, "area": forest_area})
        used += 1
    print(f"forest-masked tiles used: {used}; clusters: {sum(len(t['coords']) for t in obs_tiles):,}")

    # ring-difference g (fine-scale localization, noisy at large r)
    rm, g_obs = g_from_k(r_vals, k_pooled(obs_tiles, r_vals))
    g_sims = np.array([g_from_k(r_vals, k_pooled(sp, r_vals))[1] for sp in sim_pools])
    g_lo, g_hi = np.nanpercentile(g_sims, 2.5, axis=0), np.nanpercentile(g_sims, 97.5, axis=0)
    # PRIMARY object: cumulative L_inhom = sqrt(K/pi) (variance-stabilized, stable at all r)
    K_obs = k_pooled(obs_tiles, r_vals); L_obs = np.sqrt(np.clip(K_obs, 0, None) / np.pi)
    L_sims = np.array([np.sqrt(np.clip(k_pooled(sp, r_vals), 0, None) / np.pi) for sp in sim_pools])
    L_lo, L_hi = np.nanpercentile(L_sims, 2.5, axis=0), np.nanpercentile(L_sims, 97.5, axis=0)
    print(f"\n=== PRIMARY: L_inhom(r) - r  vs host-conditioned envelope (clustering if above) ===")
    print(f"{'r(m)':>6} {'Lobs-r':>9} {'env_lo-r':>9} {'env_hi-r':>9}  {'verdict':>22}")
    for i, r in enumerate(r_vals):
        if r == 0:
            continue
        lo_, hi_ = L_lo[i] - r, L_hi[i] - r
        v = ("ABOVE (residual foci)" if (L_obs[i] - r) > hi_
             else "below" if (L_obs[i] - r) < lo_ else "within (host-tracking)")
        print(f"{r:>6.0f} {L_obs[i]-r:>9.1f} {lo_:>9.1f} {hi_:>9.1f}  {v:>22}")
    out = {"r": rm.tolist(), "g_obs": g_obs.tolist(), "g_lo": g_lo.tolist(), "g_hi": g_hi.tolist(),
           "L_r": r_vals.tolist(), "L_obs": L_obs.tolist(), "L_lo": L_lo.tolist(), "L_hi": L_hi.tolist()}
    os.makedirs(os.path.dirname(save), exist_ok=True); json.dump(out, open(save, "w"), indent=2)
    print(f"\nsaved {save}")
    print("L_inhom-r ABOVE the host-conditioned envelope => residual aggregation beyond spruce host;"
          "\nWITHIN => mortality is host-tracking at that scale. (g table above localizes the scale.)")
    return out


def run_masked_ppm(spruce_path, landclass_path, totvol_path, forest_classes=(1,),
                   r_max=1500.0, nr=15, nsim=99, grid=60,
                   save="analysis/out/host_pcf_ppm.json"):
    """ppm DISCRIMINATOR (Vasquez): fit lambda ~ spruce + total-volume + stand-edge by a pixel-Poisson
    GLM, then regenerate the L_inhom envelope from THAT richer fitted null. If observed L_inhom still
    exceeds the envelope -> residual foci survive the richer covariates (genuine). If it collapses ->
    the apparent foci were unresolved host structure captured by total volume / edge."""
    import json, os
    from scipy import ndimage
    from sklearn.linear_model import PoissonRegressor
    from sklearn.preprocessing import StandardScaler
    rng = np.random.default_rng(42)
    d = pd.read_csv("data/clusters_typed.csv")
    tr = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id", "image_name"])
    tr = tr[tr.cluster_id != -1]
    cl_tile = tr.groupby("cluster_id").image_name.agg(lambda s: s.mode().iloc[0])
    tw = tr.groupby("image_name").agg(x0=("x","min"), x1=("x","max"), y0=("y","min"), y1=("y","max"))
    d = d[d.cluster_id.isin(cl_tile.index)].copy(); d["tile"] = d.cluster_id.map(cl_tile)
    r_vals = np.linspace(0, r_max, nr + 1)

    # ---- pass 1: per-tile forest-cell covariates + observed counts (global design matrix) ----
    tile_data, X_rows, y_rows = [], [], []
    for tile, sub in d.groupby("tile"):
        if len(sub) < 15 or tile not in tw.index:
            continue
        W = (tw.loc[tile,"x0"], tw.loc[tile,"x1"], tw.loc[tile,"y0"], tw.loc[tile,"y1"])
        gx = np.linspace(W[0], W[1], grid); gy = np.linspace(W[2], W[3], grid)
        GX, GY = np.meshgrid(gx, gy); cells = np.column_stack([GX.ravel(), GY.ravel()])
        lc = np.nan_to_num(_sample_raster(landclass_path, cells), nan=-1).astype(int)
        fmask2d = np.isin(lc.reshape(grid, grid), forest_classes)
        if fmask2d.sum() < 30:
            continue
        edge2d = ndimage.distance_transform_edt(fmask2d) * ((W[1]-W[0])/(grid-1))  # m to forest edge
        fidx = np.where(fmask2d.ravel())[0]
        fcells = cells[fidx]
        spr = np.nan_to_num(_sample_raster(spruce_path, fcells), nan=0.0)
        tot = np.nan_to_num(_sample_raster(totvol_path, fcells), nan=0.0)
        edg = edge2d.ravel()[fidx]
        # observed cluster counts per forest cell
        from scipy.spatial import cKDTree
        xy = sub[["cx","cy"]].to_numpy(float)
        nn = cKDTree(fcells).query(xy)[1]
        cnt = np.bincount(nn, minlength=len(fcells)).astype(float)
        cell_area = ((W[1]-W[0])/(grid-1)) * ((W[3]-W[2])/(grid-1))
        tile_data.append(dict(tile=tile, W=W, fcells=fcells, spr=spr, tot=tot, edg=edg,
                              cell_area=cell_area, xy=xy, nn=nn, n=len(xy)))
        X_rows.append(np.column_stack([spr, tot, edg])); y_rows.append(cnt)
    X = np.vstack(X_rows); ycount = np.concatenate(y_rows)
    Xs = StandardScaler().fit_transform(X)
    glm = PoissonRegressor(alpha=1e-6, max_iter=500).fit(Xs, ycount)
    print(f"ppm intensity model  lambda ~ spruce + total + edge  (n_cells={len(ycount):,})")
    print(f"  standardized coefficients: spruce={glm.coef_[0]:+.3f}  total={glm.coef_[1]:+.3f}  edge={glm.coef_[2]:+.3f}")
    sc = StandardScaler().fit(X)

    # ---- pass 2: fitted lambda per tile; observed + host(fitted)-conditioned sims ----
    obs_tiles, sim_pools = [], [[] for _ in range(nsim)]
    for t in tile_data:
        rel = glm.predict(sc.transform(np.column_stack([t["spr"], t["tot"], t["edg"]])))  # exp(Xb), per cell
        rel = np.clip(rel, 1e-9, None)
        n, A = t["n"], t["fcells"].shape[0] * t["cell_area"]
        norm = n / (rel.mean() * A)                       # integrate to n over forest
        lam_cell = rel * norm
        lam_floor = 0.15 * n / A
        lam_obs = np.maximum(lam_cell[t["nn"]], lam_floor)
        obs_tiles.append({"coords": t["xy"], "lam": lam_obs, "W": t["W"], "area": A})
        p = rel / rel.sum(); cs = np.sqrt(t["cell_area"])
        for s in range(nsim):
            idx = rng.choice(len(t["fcells"]), size=n, p=p)
            sxy = t["fcells"][idx] + (rng.random((n, 2)) - 0.5) * cs
            sim_pools[s].append({"coords": sxy, "lam": np.maximum(lam_cell[idx], lam_floor),
                                 "W": t["W"], "area": A})
    L_obs = np.sqrt(np.clip(k_pooled(obs_tiles, r_vals), 0, None) / np.pi)
    L_sims = np.array([np.sqrt(np.clip(k_pooled(sp, r_vals), 0, None)/np.pi) for sp in sim_pools])
    L_lo, L_hi = np.nanpercentile(L_sims, 2.5, axis=0), np.nanpercentile(L_sims, 97.5, axis=0)
    print(f"\n=== L_inhom(r)-r vs RICHER (spruce+total+edge) host-conditioned envelope ===")
    print(f"{'r(m)':>6} {'Lobs-r':>9} {'env_lo-r':>9} {'env_hi-r':>9}  {'verdict':>22}")
    out = {"r": r_vals.tolist(), "L_obs": L_obs.tolist(), "L_lo": L_lo.tolist(), "L_hi": L_hi.tolist(),
           "coef": {"spruce": glm.coef_[0], "total": glm.coef_[1], "edge": glm.coef_[2]}}
    for i, r in enumerate(r_vals):
        if r == 0: continue
        lo_, hi_ = L_lo[i]-r, L_hi[i]-r
        v = ("ABOVE (residual foci)" if (L_obs[i]-r) > hi_ else
             "below" if (L_obs[i]-r) < lo_ else "within (host-tracking)")
        print(f"{r:>6.0f} {L_obs[i]-r:>9.1f} {lo_:>9.1f} {hi_:>9.1f}  {v:>22}")
    os.makedirs(os.path.dirname(save), exist_ok=True); json.dump(out, open(save, "w"), indent=2, default=float)
    print(f"\nsaved {save}")
    print("If observed L still ABOVE this richer envelope -> residual FOCI beyond spruce+total+edge;"
          "\nif it collapses INTO it -> the apparent foci were unresolved host structure.")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true", help="single-window simulation test")
    ap.add_argument("--selftest-pooled", action="store_true", help="pooled multi-tile simulation test")
    ap.add_argument("--raster", default=None, help="path to MS-NFI spruce GeoTIFF (EPSG:3067)")
    ap.add_argument("--masked", action="store_true", help="forest-masked + pre-mortality host + envelopes")
    ap.add_argument("--ppm", action="store_true", help="ppm discriminator: lambda ~ spruce+total+edge")
    ap.add_argument("--spruce", default="data/kuusi_vmi1x_1721.tif")
    ap.add_argument("--landclass", default="data/maaluokka_vmi1x_1721.tif")
    ap.add_argument("--totvol", default="data/tilavuus_vmi1x_1721.tif")
    ap.add_argument("--bbox-mode", action="store_true", help="(deprecated) per-block bbox window")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    elif a.selftest_pooled:
        selftest_pooled()
    elif a.ppm:
        run_masked_ppm(a.spruce, a.landclass, a.totvol)
    elif a.masked:
        run_masked(a.spruce, a.landclass)
    elif a.raster and a.bbox_mode:
        run_blocks(raster=a.raster)
    elif a.raster:
        run_pooled(a.raster)
    else:
        run_blocks(raster=None)
