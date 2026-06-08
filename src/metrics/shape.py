"""Shape and morphology metrics for individual dead-tree clusters."""

import warnings
import numpy as np
import pandas as pd
from scipy.stats import linregress
from shapely.geometry import box as shapely_box
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Fractal dimension (box-counting)
# ---------------------------------------------------------------------------

def fractal_dimension(geom, min_scale=0.1, max_scale=100.0, num_scales=10):
    """Box-counting fractal dimension of a Shapely geometry."""
    if geom is None or geom.is_empty:
        return np.nan
    scales = np.logspace(np.log10(min_scale), np.log10(max_scale), num_scales)
    x_min, y_min, x_max, y_max = geom.bounds
    counts = []
    for s in scales:
        x_bins = np.arange(x_min, x_max + s, s)
        y_bins = np.arange(y_min, y_max + s, s)
        count = sum(
            1
            for i in range(len(x_bins) - 1)
            for j in range(len(y_bins) - 1)
            if geom.intersects(shapely_box(x_bins[i], y_bins[j], x_bins[i + 1], y_bins[j + 1]))
        )
        counts.append(count)
    valid = np.array(counts) > 0
    if valid.sum() < 2:
        return np.nan
    slope, *_ = linregress(np.log(scales[valid]), np.log(np.array(counts)[valid]))
    return -slope


# ---------------------------------------------------------------------------
# Per-cluster shape metrics
# ---------------------------------------------------------------------------

def _alpha_shape(points, alpha=0.5):
    """Return an alphashape Polygon for an (n,2) array of points."""
    try:
        import alphashape as _as
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shape = _as.alphashape(points, alpha)
        if shape is None or shape.is_empty:
            return None
        if shape.geom_type == "MultiPolygon":
            shape = unary_union(shape)
        return shape if shape.geom_type == "Polygon" else None
    except Exception:
        return None


def _mbr_props(geom):
    """Aspect ratio and orientation (degrees) from minimum rotated rectangle."""
    rect = geom.minimum_rotated_rectangle
    coords = list(rect.exterior.coords)
    sides = [
        np.hypot(coords[i + 1][0] - coords[i][0], coords[i + 1][1] - coords[i][1])
        for i in range(4)
    ]
    length = max(sides)
    width = max(min(sides), 1.0)   # guard against degenerate rectangles
    aspect_ratio = length / width
    dx = coords[1][0] - coords[0][0]
    dy = coords[1][1] - coords[0][1]
    orientation = np.degrees(np.arctan2(dy, dx)) % 180
    return aspect_ratio, orientation


def compute_cluster_shape(points, alpha=0.5, fd_min_scale=0.1, fd_max_scale=100.0):
    """Compute all shape metrics for one cluster given its (n,2) point array.

    Returns a dict with: area_m2, perimeter_m, compactness, convexity,
    aspect_ratio, elongation, orientation, fractal_dim, alpha_compactness,
    core_to_edge_ratio.
    """
    if len(points) < 3:
        return None

    shape = _alpha_shape(points, alpha)
    if shape is None:
        return None

    area = shape.area
    perimeter = shape.length
    if area <= 0 or perimeter <= 0:
        return None

    convex = shape.convex_hull
    compactness = (4 * np.pi * area) / (perimeter ** 2)
    convexity = area / convex.area if convex.area > 0 else np.nan
    aspect_ratio, orientation = _mbr_props(shape)
    elongation = 1.0 - (1.0 / aspect_ratio) if aspect_ratio > 0 else np.nan

    # alpha_compactness = concavity_ratio = alpha_shape_area / convex_hull_area
    alpha_compactness = area / convex.area if convex.area > 0 else np.nan

    core_area = max(shape.buffer(-1.0).area, 0.0)
    core_to_edge = core_area / area

    fd = fractal_dimension(shape, min_scale=fd_min_scale, max_scale=fd_max_scale)

    return {
        "area_m2": area,
        "perimeter_m": perimeter,
        "compactness": compactness,
        "convexity": convexity,
        "aspect_ratio": aspect_ratio,
        "elongation": elongation,
        "orientation": orientation,
        "fractal_dim": fd,
        "alpha_compactness": alpha_compactness,
        "core_to_edge_ratio": core_to_edge,
    }


def compute_all_shape_metrics(df, cluster_col="cluster_id", alpha=0.5):
    """Apply compute_cluster_shape to every cluster in df.

    df must have columns: x, y, and cluster_col.
    Returns a DataFrame indexed by cluster id.
    """
    rows = []
    for cid, grp in df.groupby(cluster_col):
        pts = grp[["x", "y"]].values
        m = compute_cluster_shape(pts, alpha=alpha)
        if m is not None:
            m[cluster_col] = cid
            m["n_points"] = len(pts)
            m["cx"] = pts[:, 0].mean()
            m["cy"] = pts[:, 1].mean()
            rows.append(m)
    return pd.DataFrame(rows).set_index(cluster_col)
