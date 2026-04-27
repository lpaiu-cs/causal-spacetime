"""Minkowski Alexandrov interval geometry utilities."""

from __future__ import annotations

import math


def _validate_spacetime_dim(spacetime_dim: int) -> int:
    dimension = int(spacetime_dim)
    if dimension < 2:
        raise ValueError("spacetime_dim must be at least 2")
    return dimension


def alexandrov_eta(spacetime_dim: int) -> float:
    """Return ``eta_D`` for ``V_D(tau) = eta_D * tau**D`` in Minkowski space."""

    dimension = _validate_spacetime_dim(spacetime_dim)
    spatial_dim = dimension - 1
    numerator = math.pi ** (spatial_dim / 2.0)
    denominator = (
        dimension
        * (2.0 ** spatial_dim)
        * math.gamma((dimension + 1.0) / 2.0)
    )
    return numerator / denominator


def alexandrov_volume(tau: float, spacetime_dim: int) -> float:
    """Return D-dimensional Alexandrov interval volume for proper time ``tau``."""

    tau_value = float(tau)
    if tau_value < 0.0:
        raise ValueError("tau must be non-negative")
    dimension = _validate_spacetime_dim(spacetime_dim)
    return alexandrov_eta(dimension) * (tau_value**dimension)


def estimate_tau_from_interval_count(
    interval_count: int,
    rho: float,
    spacetime_dim: int,
) -> float:
    """Estimate proper time from interval cardinality in D-dimensional space."""

    count = int(interval_count)
    if count < 0:
        raise ValueError("interval_count must be non-negative")
    if rho <= 0:
        raise ValueError("rho must be positive")
    dimension = _validate_spacetime_dim(spacetime_dim)
    eta = alexandrov_eta(dimension)
    return (count / (rho * eta)) ** (1.0 / dimension)
