"""
Bridge for the spatstat leg: export per-block covariate grids + points + tile windows to CSV,
so the R script needs only spatstat (no terra/GDAL). Uses rasterio windowed reads (fast).

Per block <b> writes to analysis/out/rblocks/:
  block<b>_grid.csv   x, y, spr, tot, lc     (RES-m grid over the block bbox; lc = land class)
  block<b>_pts.csv    cx, cy                 (cluster centroids in the block)
  block<b>_tiles.csv  x0, x1, y0, y1         (6 km tile windows)
"""
import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
from sklearn.cluster import DBSCAN

RES = 100  # export grid resolution (m); 16 m native aggregated to this
NODATA = 32767
SPRUCE, TOTAL, LAND = ("data/kuusi_vmi1x_1721.tif", "data/tilavuus_vmi1x_1721.tif",
                       "data/maaluokka_vmi1x_1721.tif")
OUT = "analysis/out/rblocks"
os.makedirs(OUT, exist_ok=True)


def block_grid(path, x0, x1, y0, y1, categorical=False):
    """Windowed read of [x0,x1]x[y0,y1], aggregated to RES-m; returns (xs, ys, Z) with NaN nodata."""
    with rasterio.open(path) as s:
        win = from_bounds(x0, y0, x1, y1, s.transform)
        arr = s.read(1, window=win).astype(float)
        arr[arr == NODATA] = np.nan
        wt = s.window_transform(win)
        nyr, nxr = arr.shape
        xs0 = wt.c + wt.a / 2 + np.arange(nxr) * wt.a
        ys0 = wt.f + wt.e / 2 + np.arange(nyr) * wt.e   # wt.e < 0 (north-up)
    fac = max(1, round(RES / 16))
    ny, nx = nyr // fac, nxr // fac
    arr = arr[:ny * fac, :nx * fac].reshape(ny, fac, nx, fac)
    if categorical:
        # modal over the block (ignore nan): use rounded mean as a cheap mode proxy for class 1 dominance
        Z = np.nanmedian(arr, axis=(1, 3))
    else:
        Z = np.nanmean(arr, axis=(1, 3))
    xs = xs0[:nx * fac].reshape(nx, fac).mean(1)
    ys = ys0[:ny * fac].reshape(ny, fac).mean(1)
    return xs, ys, Z


def main():
    ct = pd.read_csv("data/clusters_typed.csv")[["cluster_id", "cx", "cy"]]
    ct["block"] = DBSCAN(eps=15000, min_samples=1).fit_predict(ct[["cx", "cy"]].values)
    tr = pd.read_csv("data/trees_clustered.csv", usecols=["x", "y", "cluster_id", "image_name"])
    tr = tr[tr.cluster_id != -1]
    cl_tile = tr.groupby("cluster_id").image_name.agg(lambda s: s.mode().iloc[0])
    tw = tr.groupby("image_name").agg(x0=("x", "min"), x1=("x", "max"),
                                      y0=("y", "min"), y1=("y", "max"))
    ct = ct[ct.cluster_id.isin(cl_tile.index)].copy()
    ct["tile"] = ct.cluster_id.map(cl_tile)

    for b, sub in ct.groupby("block"):
        if len(sub) < 60:
            continue
        tiles = tw.loc[sorted(sub.tile.unique())]
        x0, x1 = tiles.x0.min() - 200, tiles.x1.max() + 200
        y0, y1 = tiles.y0.min() - 200, tiles.y1.max() + 200
        xs, ys, spr = block_grid(SPRUCE, x0, x1, y0, y1)
        _, _, tot = block_grid(TOTAL, x0, x1, y0, y1)
        _, _, lc = block_grid(LAND, x0, x1, y0, y1, categorical=True)
        ny, nx = spr.shape
        GX, GY = np.meshgrid(xs[:nx], ys[:ny])
        grid = pd.DataFrame({"x": GX.ravel(), "y": GY.ravel(),
                             "spr": np.nan_to_num(spr).ravel(),
                             "tot": np.nan_to_num(tot).ravel(),
                             "lc": np.nan_to_num(lc, nan=-1).ravel()})
        grid.to_csv(f"{OUT}/block{b}_grid.csv", index=False)
        sub[["cx", "cy"]].to_csv(f"{OUT}/block{b}_pts.csv", index=False)
        tiles.reset_index()[["x0", "x1", "y0", "y1"]].to_csv(f"{OUT}/block{b}_tiles.csv", index=False)
        print(f"block {b}: pts={len(sub)} tiles={len(tiles)} grid={ny}x{nx} "
              f"forest%={100*np.mean(np.round(np.nan_to_num(lc,nan=-1))==1):.0f}")
    print(f"exported to {OUT}/")


if __name__ == "__main__":
    main()
