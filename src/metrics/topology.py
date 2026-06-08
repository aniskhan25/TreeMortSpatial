"""Topological metrics for dead-tree clusters."""

import numpy as np


def euler_number(coords):
    """Euler characteristic χ = V - E + F via Delaunay triangulation.

    For a simply-connected cluster χ = 1; holes decrease it.

    Parameters
    ----------
    coords : (n, 2) array of point positions within a cluster

    Returns
    -------
    int or nan
    """
    from scipy.spatial import Delaunay

    if len(coords) < 3:
        return np.nan

    try:
        tri = Delaunay(coords)
    except Exception:
        return np.nan

    V = len(coords)
    F = len(tri.simplices)

    # Collect unique edges from the triangulation
    edges = set()
    for simplex in tri.simplices:
        for i in range(3):
            edge = tuple(sorted((simplex[i], simplex[(i + 1) % 3])))
            edges.add(edge)
    E = len(edges)

    return V - E + F
