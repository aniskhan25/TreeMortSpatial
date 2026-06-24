"""
Restyle the notebook-generated MAP figures to the house style: clean CartoDB Positron
basemap, palette type-colours, gradient-framing labels, tidy legends/scalebars/insets.

Figures: fig1_spatial_extent, fig6_typology_map, fig8_local_morans,
         fig_contagion_a, fig_contagion_b.  Needs network tiles + geopandas.
"""
import sys, warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import contextily as ctx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib_scalebar.scalebar import ScaleBar
from figstyle import apply_style, PALETTE, save

sys.path.insert(0, "src")
from metrics.spatial import local_morans_i  # noqa: E402

apply_style()
P = PALETTE
warnings.filterwarnings("ignore")

TYPE_COLOR  = {"Type0": P["secondary"], "Type1": P["primary"]}
TYPE_MARKER = {"Type0": "^",            "Type1": "o"}
TYPE_LABEL  = {"Type0": "Type 0 (elongated/dispersed)",
               "Type1": "Type 1 (compact/aggregated)"}
BASEMAP = ctx.providers.CartoDB.Positron
BBOX = dict(xmin=397000, xmax=421000, ymin=6715000, ymax=6740750)

df = pd.read_csv("data/clusters_typed.csv")


def to3857(d):
    return gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d.cx, d.cy),
                            crs="EPSG:3067").to_crs(epsg=3857)


def in_bbox(d):
    return d[(d.cx >= BBOX["xmin"]) & (d.cx <= BBOX["xmax"]) &
             (d.cy >= BBOX["ymin"]) & (d.cy <= BBOX["ymax"])]


def add_locator(ax, all_gdf):
    ins = ax.inset_axes([0.66, 0.66, 0.32, 0.32])
    ins.set_facecolor("white")
    all_gdf.plot(ax=ins, color="#9AB0BE", markersize=0.3, alpha=0.6)
    ins.set_aspect("equal"); ins.set_axis_off()
    for s in ins.spines.values():
        s.set_visible(True); s.set_edgecolor("#bbb"); s.set_linewidth(0.6)


def basemap(ax):
    ctx.add_basemap(ax, source=BASEMAP, attribution_size=5)


def scalebar(ax, loc="lower right"):
    ax.add_artist(ScaleBar(1, units="m", location=loc, box_alpha=0.6,
                           color="#333", frameon=True, font_properties={"size": 8}))


# ---- fig1: study-region extent (all 14,582 filtered clusters) -----------------
def fig1():
    g = to3857(df)
    fig, ax = plt.subplots(figsize=(8, 8.4), dpi=150)
    g.plot(ax=ax, color=P["primary"], markersize=2, alpha=0.35, edgecolors="none")
    basemap(ax); scalebar(ax)
    ax.set_axis_off()
    ax.set_title(f"Dead-tree cluster centroids  ($n$={len(g):,})", fontsize=11, pad=8)
    fig.tight_layout(); save(fig, "figures/fig1_spatial_extent"); plt.close(fig)


# ---- fig6: typology map -------------------------------------------------------
def fig6():
    z = to3857(in_bbox(df.dropna(subset=["morph_type"])))
    allg = to3857(df.dropna(subset=["morph_type"]))
    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
    for t in ("Type1", "Type0"):
        sub = z[z.morph_type == t]
        sub.plot(ax=ax, color=TYPE_COLOR[t], markersize=10, alpha=0.7,
                 marker=TYPE_MARKER[t], edgecolors="white", linewidth=0.2)
    basemap(ax); scalebar(ax)
    handles = [mlines.Line2D([], [], color=TYPE_COLOR[t], marker=TYPE_MARKER[t],
               ls="", markersize=8, label=TYPE_LABEL[t]) for t in ("Type0", "Type1")]
    ax.legend(handles=handles, loc="lower left", framealpha=0.9)
    ax.set_axis_off()
    ax.set_title(f"Gradient endpoints — study region  ($n$={len(z):,})", fontsize=11, pad=8)
    add_locator(ax, allg)
    fig.tight_layout(); save(fig, "figures/fig6_typology_map"); plt.close(fig)


# ---- fig8: LISA hotspots (fractal dimension) ----------------------------------
def fig8():
    sub = df.dropna(subset=["fractal_dim"]).copy()
    lm = local_morans_i(sub[["cx", "cy"]].values, sub["fractal_dim"].values, threshold=1000.0)
    sub["lm"] = lm["label"].values
    z = to3857(in_bbox(sub)); allg = to3857(sub)
    col = {"HH": P["accent"], "LL": P["primary"], "HL": "#E8A33D",
           "LH": P["muted"], "NS": "#D5D5D5"}
    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
    for lbl in ["NS", "LH", "HL", "LL", "HH"]:
        pts = z[z.lm == lbl]
        if len(pts):
            pts.plot(ax=ax, color=col[lbl], markersize=3 if lbl == "NS" else 9,
                     alpha=0.3 if lbl == "NS" else 0.85, marker="o",
                     edgecolors="none", label=f"{lbl}  ($n$={len(pts):,})")
    basemap(ax); scalebar(ax, "lower left")
    ax.set_axis_off()
    ax.set_title("Local Moran's $I$ — fractal-dimension hot/cold spots", fontsize=11, pad=8)
    ax.legend(title="LISA class", loc="lower right", framealpha=0.9)
    add_locator(ax, allg)
    fig.tight_layout(); save(fig, "figures/fig8_local_morans"); plt.close(fig)


# ---- fig_contagion_a / _b -----------------------------------------------------
def fig_contagion():
    z = to3857(in_bbox(df.dropna(subset=["morph_type"])))
    allg = to3857(df.dropna(subset=["morph_type"]))
    thr = df["betweenness_200m"].quantile(0.99)

    # (a) reach choropleth
    fig, ax = plt.subplots(figsize=(8, 7.4), dpi=150)
    z.plot(ax=ax, column="reach_200m", cmap="YlOrRd", markersize=10, alpha=0.8,
           legend=True, legend_kwds={"label": "Reach at 200 m (component size)", "shrink": 0.55},
           edgecolors="none")
    basemap(ax); scalebar(ax, "lower left")
    ax.set_axis_off()
    ax.set_title("Connectivity reach at 200 m", fontsize=11, pad=8)
    add_locator(ax, allg)
    fig.tight_layout(); save(fig, "figures/fig_contagion_a"); plt.close(fig)

    # (b) high-betweenness bridges (no outbreak/sentinel framing)
    bridges = z[z.betweenness_200m >= thr]
    fig, ax = plt.subplots(figsize=(8, 7.4), dpi=150)
    z.plot(ax=ax, color="#C9C9C9", markersize=5, alpha=0.5, edgecolors="none")
    bridges.plot(ax=ax, color=P["accent"], markersize=16, alpha=0.9, marker="D",
                 edgecolors="white", linewidth=0.3,
                 label=f"Top 1% betweenness ($n$={len(bridges):,})")
    basemap(ax); scalebar(ax, "lower left")
    ax.set_axis_off()
    ax.legend(loc="lower right", framealpha=0.9)
    ax.set_title("High-betweenness bridge clusters\n(geometric articulation points; no temporal validation)",
                 fontsize=10.5, pad=8)
    add_locator(ax, allg)
    fig.tight_layout(); save(fig, "figures/fig_contagion_b"); plt.close(fig)


if __name__ == "__main__":
    fig1();         print("fig1_spatial_extent")
    fig6();         print("fig6_typology_map")
    fig8();         print("fig8_local_morans")
    fig_contagion(); print("fig_contagion_a / _b")
    print("maps done")
