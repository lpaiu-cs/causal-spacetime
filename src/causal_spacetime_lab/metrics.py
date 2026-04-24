"""Metric calculations and estimators used by reconstruction experiments."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.constants import ETA_1P1_CAUSAL_DIAMOND


def minkowski_tau_1p1(
    p: ArrayLike,
    q: ArrayLike,
    *,
    atol: float = 1e-12,
) -> float:
    """Return 1+1D Minkowski proper time between timelike or null events."""

    p_event = np.asarray(p, dtype=float)
    q_event = np.asarray(q, dtype=float)
    if p_event.shape != (2,) or q_event.shape != (2,):
        raise ValueError("p and q must each have shape (2,), representing (t, x)")

    dt = float(q_event[0] - p_event[0])
    dx = float(q_event[1] - p_event[1])
    interval_squared = dt * dt - dx * dx
    if interval_squared < -atol:
        raise ValueError("events are spacelike-separated; proper time is undefined")
    return math.sqrt(max(interval_squared, 0.0))


def causal_diamond_volume_1p1(T: float) -> float:
    """Return the 1+1D continuum volume of a causal diamond of duration ``T``."""

    if T <= 0:
        raise ValueError("T must be positive")
    return 0.5 * T * T


def estimate_tau_from_interval_count(
    count: int | float,
    rho: float,
    eta_d: float = ETA_1P1_CAUSAL_DIAMOND,
    d: int | float = 2,
) -> float:
    """Estimate proper time from causal interval cardinality and event density.

    ``rho`` supplies metric scale information. The causal order alone provides
    the interval count, but not the density scale.
    """

    if count < 0:
        raise ValueError("count must be non-negative")
    if rho <= 0:
        raise ValueError("rho must be positive")
    if eta_d <= 0:
        raise ValueError("eta_d must be positive")
    if d <= 0:
        raise ValueError("d must be positive")

    volume = float(count) / rho
    return (volume / eta_d) ** (1.0 / float(d))

