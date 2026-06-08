"""Composite Severity Index (CSI) for dead-tree clusters."""

import numpy as np


def composite_severity_index(fractal_dim, clark_evans, n_points):
    """CSI = fractal_dim × (1 / clark_evans) × n_points.

    Higher values indicate greater ecological impact: complex shape,
    tight internal aggregation, and large cluster size.

    Parameters
    ----------
    fractal_dim  : float or array
    clark_evans  : float or array  (must be > 0)
    n_points     : float or array

    Returns
    -------
    float or array — NaN where clark_evans ≤ 0
    """
    ce = np.asarray(clark_evans, dtype=float)
    fd = np.asarray(fractal_dim, dtype=float)
    n = np.asarray(n_points, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        csi = np.where(ce > 0, fd * (1.0 / ce) * n, np.nan)
    return float(csi) if csi.ndim == 0 else csi
