import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN


def run_dbscan(df, eps=20, min_samples=3):
    coords = df[["x", "y"]].values
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coords)
    df = df.copy()
    df["cluster_id"] = labels
    return df


def filter_clusters(df, min_trees=5, min_area=100.0, max_elongation=10.0):
    """Apply the paper's filtering criteria to a cluster-level metrics DataFrame.

    Expects columns: n_points, area_m2, aspect_ratio.
    Returns the filtered DataFrame.
    """
    mask = (
        (df["n_points"] >= min_trees)
        & (df["area_m2"] >= min_area)
        & (df["aspect_ratio"] <= max_elongation)
    )
    return df[mask].copy()


def compute_nnd(coords):
    """Return per-point nearest-neighbour distance for an (n,2) array."""
    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=2)
    return distances[:, 1]


def compute_density_within_radius(coords, radius):
    """Return point density (trees/m²) within `radius` for each point."""
    tree = cKDTree(coords)
    area = np.pi * radius ** 2
    counts = np.array([len(tree.query_ball_point(p, radius)) for p in coords])
    return counts / area


def cluster_nnd_stats(coords):
    """Mean, SD, and skewness of NND for a single cluster's point array."""
    from scipy.stats import skew
    nnd = compute_nnd(coords)
    return {"nnd_mean": nnd.mean(), "nnd_std": nnd.std(), "nnd_skew": float(skew(nnd))}
