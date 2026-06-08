"""Spatial statistics for dead-tree cluster centroids."""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import skew
from joblib import Parallel, delayed
from tqdm import tqdm


# ---------------------------------------------------------------------------
# NND-based metrics
# ---------------------------------------------------------------------------

def nnd_array(coords):
    """Per-point nearest-neighbour distances for an (n,2) centroid array."""
    tree = cKDTree(coords)
    d, _ = tree.query(coords, k=2)
    return d[:, 1]


def nnd_stats(coords):
    """Mean, SD, and skewness of NND between cluster centroids."""
    d = nnd_array(coords)
    return {"nnd_mean": d.mean(), "nnd_std": d.std(), "nnd_skew": float(skew(d))}


# ---------------------------------------------------------------------------
# Clark-Evans Index
# ---------------------------------------------------------------------------

def clark_evans_index(coords, area):
    """Clark-Evans aggregation index R for a set of points in a study area.

    R < 1 → clustering, R > 1 → dispersion, R ≈ 1 → random.

    Parameters
    ----------
    coords : (n, 2) array  — point coordinates
    area   : float         — area of the study window in the same units²
    """
    n = len(coords)
    if n < 2:
        return np.nan
    mean_nnd = nnd_array(coords).mean()
    expected = 0.5 * np.sqrt(area / n)
    return mean_nnd / expected if expected > 0 else np.nan


# ---------------------------------------------------------------------------
# Moran's I
# ---------------------------------------------------------------------------

def local_morans_i(coords, values, threshold=1000.0):
    """Local Moran's I using libpysal / esda.

    Returns a DataFrame with columns: Ii, p_sim, label (HH/LL/HL/LH/NS).
    """
    from libpysal.weights import DistanceBand
    from esda.moran import Moran_Local

    w = DistanceBand.from_array(coords, threshold=threshold)
    w.transform = "r"
    lm = Moran_Local(values, w)
    labels = np.full(len(values), "NS", dtype=object)
    sig = lm.p_sim < 0.05
    labels[sig & (lm.q == 1)] = "HH"
    labels[sig & (lm.q == 3)] = "LL"
    labels[sig & (lm.q == 2)] = "HL"
    labels[sig & (lm.q == 4)] = "LH"
    return pd.DataFrame({"Ii": lm.Is, "p_sim": lm.p_sim, "label": labels})


# ---------------------------------------------------------------------------
# Ripley's K function
# ---------------------------------------------------------------------------

def _k_for_plot(group, support, plot_size):
    """Compute Ripley's K statistic for one spatial plot (group of centroids)."""
    from pointpats.distance_statistics import k_test
    from shapely.geometry import Point, Polygon

    if len(group) < 3:
        return None
    xs, ys = group[:, 0], group[:, 1]
    cx, cy = xs.mean(), ys.mean()
    half = plot_size / 2
    window = Polygon([
        (cx - half, cy - half), (cx - half, cy + half),
        (cx + half, cy + half), (cx + half, cy - half),
    ])
    result = k_test(group, support=support, hull=window)
    return result.statistic


def ripley_k(coords, plot_ids, max_dist=1500, n_intervals=50,
             plot_size=6000, n_jobs=-1):
    """Compute aggregated Ripley's K across spatial plots.

    Parameters
    ----------
    coords    : (n, 2) array of centroid coordinates
    plot_ids  : (n,) array of integer plot labels
    max_dist  : maximum distance for K function (metres)
    n_intervals : number of distance bins
    plot_size : side length of the square study window per plot (metres)
    n_jobs    : joblib parallel workers (-1 = all CPUs)

    Returns
    -------
    d_vals    : (n_intervals,) distance values
    k_obs_mean: (n_intervals,) mean observed K across plots
    """
    support = (0, max_dist, n_intervals)
    groups = {pid: coords[plot_ids == pid] for pid in np.unique(plot_ids)}

    k_list = Parallel(n_jobs=n_jobs)(
        delayed(_k_for_plot)(g, support, plot_size)
        for g in tqdm(groups.values(), desc="Ripley K plots")
    )
    k_list = [k for k in k_list if k is not None]
    if not k_list:
        raise ValueError("No plots had enough points to compute K.")

    d_vals = np.linspace(0, max_dist, n_intervals)
    k_obs_mean = np.mean(k_list, axis=0)
    return d_vals, k_obs_mean


# ---------------------------------------------------------------------------
# Mark Correlation Function
# ---------------------------------------------------------------------------

def mark_correlation_function(coords, marks, r_values, dr=None):
    """Mark Correlation Function k_m(r) for a set of labelled points.

    k_m(r) = E[m_i * m_j | dist(i,j) ≈ r] / E[m]²

    Parameters
    ----------
    coords   : (n, 2) centroid array
    marks    : (n,) mark values (e.g., fractal_dim per cluster)
    r_values : 1-D array of distances at which to evaluate k_m
    dr       : half-width of distance bin (default: half the smallest r step)

    Returns
    -------
    km : (len(r_values),) MCF values; NaN where no pairs fall in a bin
    """
    if dr is None:
        dr = (r_values[1] - r_values[0]) / 2.0 if len(r_values) > 1 else r_values[0] / 2.0

    marks = np.asarray(marks, dtype=float)
    mean_mark_sq = marks.mean() ** 2
    if mean_mark_sq == 0:
        return np.full(len(r_values), np.nan)

    tree = cKDTree(coords)
    km = np.full(len(r_values), np.nan)

    for k, r in enumerate(r_values):
        pairs = tree.query_pairs(r + dr)
        if not pairs:
            continue
        close = [
            (i, j) for i, j in pairs
            if np.linalg.norm(coords[i] - coords[j]) >= r - dr
        ]
        if not close:
            continue
        products = np.array([marks[i] * marks[j] for i, j in close])
        km[k] = products.mean() / mean_mark_sq

    return km
