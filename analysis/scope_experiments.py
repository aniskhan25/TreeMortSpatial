"""
Scope confirmation + the two reviewer experiments on the corrected national scale.

(0) CONFIRM sampling spans Finland: tile/cluster locations, lat/lon ranges, per-band
    counts, and a national-extent map (figures/fig_sampling_extent.png).
(1) Tanaka R4, latitude-conditioned: does the detector confound (AUC 0.727) survive once
    latitude is partialled out? Report AUC(detection only), AUC(latitude only),
    AUC(detection | latitude residualised).
(2) Vasquez: per-threshold neighbour / island counts for the Moran subsample, to scope the
    distance-band weights (esp. whether the 5 km column is island-dominated).
"""
import sys, warnings, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import cross_val_score
warnings.filterwarnings("ignore")
sys.path.insert(0, "analysis")
from figstyle import apply_style, PALETTE, save
apply_style(); P = PALETTE
SEED = 42

typed = pd.read_csv("data/clusters_typed.csv")
trees = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id", "max_value"])
trees = trees[trees.cluster_id.isin(set(typed.cluster_id))]

def lon(x): return 27 + (x - 500000) / (111320 * math.cos(math.radians(63)))
def lat(y): return 60 + (y - 6660000) / 111320.0

# ---- (0) CONFIRM spatial extent ----------------------------------------------
print("=" * 60, "\n(0) SAMPLING EXTENT CONFIRMATION\n", "=" * 60)
print(f"distinct aerial tiles: {pd.read_csv('data/trees_clustered.csv', usecols=['image_name']).image_name.nunique()}")
print(f"clusters: {len(typed):,}")
print(f"longitude: {lon(typed.cx.min()):.1f}E .. {lon(typed.cx.max()):.1f}E (full national width)")
print(f"latitude : {lat(typed.cy.min()):.1f}N .. {lat(typed.cy.max()):.1f}N (south boreal -> Lapland)")
def band(cy): return "S (~60-62N)" if cy < 7e6 else ("C (~64-65N)" if cy < 7.35e6 else "N (~67-68N)")
typed["band"] = typed.cy.map(band)
print("\nper-band cluster counts:"); print(typed.band.value_counts().to_string())
occ = sorted(set((typed.cy // 50000 * 50).astype(int)))
print(f"\noccupied 50km northing bands (km): {occ}")

# national-extent map
try:
    import geopandas as gpd, contextily as ctx
    g = gpd.GeoDataFrame(typed, geometry=gpd.points_from_xy(typed.cx, typed.cy),
                         crs="EPSG:3067").to_crs(epsg=3857)
    fig, ax = plt.subplots(figsize=(6.5, 9), dpi=150)
    g.plot(ax=ax, color=P["accent"], markersize=3, alpha=0.35, edgecolors="none")
    # pad the view so Finland's outline is visible around the points
    xmin, ymin, xmax, ymax = g.total_bounds
    padx = (xmax - xmin) * 0.25; pady = (ymax - ymin) * 0.08
    ax.set_xlim(xmin - padx, xmax + padx); ax.set_ylim(ymin - pady, ymax + pady)
    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, attribution_size=5)
    # scale bar: data are EPSG:3857 (Web Mercator), so 1 map unit = cos(lat) ground metres;
    # correct at the view's centre latitude so "150 km" is true ground distance.
    from matplotlib_scalebar.scalebar import ScaleBar
    yc = sum(ax.get_ylim()) / 2
    latc = math.degrees(math.atan(math.sinh(yc / 6378137.0)))
    ax.add_artist(ScaleBar(math.cos(math.radians(latc)), units="m", location="lower right",
                           fixed_value=150, fixed_units="km", box_alpha=0.6, color="#333",
                           frameon=True, font_properties={"size": 8}))
    ax.set_axis_off()
    ax.set_title("Sampled cluster locations across Finland\n"
                 "($n$=14,582 clusters; 312 aerial tiles; ~60.4–67.8°N)",
                 fontsize=11, pad=8)
    save(fig, "figures/fig_sampling_extent")
    plt.close(fig)
    print("\nmap -> figures/fig_sampling_extent.png (+ .pdf)")
except Exception as e:
    print(f"\n(map skipped: {e})")

# ---- (1) Tanaka R4 latitude-conditioned --------------------------------------
print("\n" + "=" * 60, "\n(1) DETECTOR CONFOUND vs LATITUDE (Tanaka R4)\n", "=" * 60)
conf = trees.groupby("cluster_id").max_value.agg(["mean", "std"]).rename(
    columns={"mean": "conf_mean", "std": "conf_std"})
d = typed.merge(conf, on="cluster_id", how="left")
d["density"] = d.n_points / d.area_m2
y = (d.morph_type == "Type0").astype(int).to_numpy()
det_cols = ["conf_mean", "conf_std", "n_points", "density"]
det = d[det_cols].fillna(d[det_cols].median())
latc = d[["cx", "cy"]].fillna(d[["cx", "cy"]].median())

def auc(X):
    return cross_val_score(LogisticRegression(max_iter=2000),
                           StandardScaler().fit_transform(X), y, cv=5, scoring="roc_auc")

a_det = auc(det); a_lat = auc(latc)
# residualise each detection covariate on latitude+longitude, then AUC of residuals
res = det.copy()
XY = StandardScaler().fit_transform(latc)
for c in det_cols:
    res[c] = det[c].values - LinearRegression().fit(XY, det[c].values).predict(XY)
a_res = auc(res)
a_both = auc(pd.concat([det, latc], axis=1))
print(f"AUC(detection covariates only)        = {a_det.mean():.3f} +/- {a_det.std():.3f}  [the 0.727 confound]")
print(f"AUC(latitude+longitude only)          = {a_lat.mean():.3f} +/- {a_lat.std():.3f}  [does location alone predict type?]")
print(f"AUC(detection | location residualised)= {a_res.mean():.3f} +/- {a_res.std():.3f}  [confound AFTER removing latitude]")
print(f"AUC(detection + location together)    = {a_both.mean():.3f} +/- {a_both.std():.3f}")
print("Type0 prevalence by band:")
print((100 * typed.groupby('band', observed=True).apply(lambda g:(g.morph_type=='Type0').mean())).round(1).to_string())
verdict = ("confound is LARGELY latitudinal" if a_res.mean() < 0.6
           else "confound SURVIVES latitude conditioning (not just a latitudinal gradient)")
print(f"-> {verdict}")

# ---- (2) Vasquez: Moran per-threshold neighbour/island counts ----------------
print("\n" + "=" * 60, "\n(2) MORAN DISTANCE-BAND SUPPORT (Vasquez)\n", "=" * 60)
from libpysal.weights import DistanceBand
t0 = typed[typed.morph_type == "Type0"]; t1 = typed[typed.morph_type == "Type1"]
rng = np.random.default_rng(SEED)
n0 = int(2000 * len(t0) / len(typed)); n1 = 2000 - n0
sub = pd.concat([t0.iloc[rng.choice(len(t0), n0, replace=False)],
                 t1.iloc[rng.choice(len(t1), n1, replace=False)]])
coords = sub[["cx", "cy"]].values
print(f"subsample n = {len(coords)}")
print(f"{'thr':>7s} {'islands':>9s} {'with_nbr':>9s} {'%island':>8s} {'mean_deg':>9s}")
for th in [200, 500, 1000, 2000, 5000]:
    w = DistanceBand.from_array(coords, threshold=th, silence_warnings=True)
    isl = len(w.islands); withn = len(coords) - isl
    degs = [len(w.neighbors[i]) for i in w.neighbors]
    print(f"{th:7d} {isl:9d} {withn:9d} {100*isl/len(coords):7.1f}% {np.mean(degs):9.1f}")
print("(islands DROP as threshold grows; if the 5 km row is NOT island-dominated, the column is usable)")
