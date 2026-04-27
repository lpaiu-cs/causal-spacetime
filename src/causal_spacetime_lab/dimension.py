"""Order-statistical dimension estimators for finite causal sets."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike


def relation_fraction(causal_matrix: ArrayLike) -> float:
    """Return ordered relation fraction ``R / (N * (N - 1))``."""

    matrix = np.asarray(causal_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("causal_matrix must be square")
    n_events = matrix.shape[0]
    if n_events < 2:
        raise ValueError("at least two events are required")
    return float(np.count_nonzero(matrix) / (n_events * (n_events - 1)))


def myrheim_meyer_relation_fraction(spacetime_dim: float) -> float:
    """Return the expected Myrheim-Meyer ordered relation fraction curve."""

    dimension = float(spacetime_dim)
    if dimension <= 0.0:
        raise ValueError("spacetime_dim must be positive")
    log_value = (
        math.lgamma(dimension + 1.0)
        + math.lgamma(dimension / 2.0)
        - math.log(4.0)
        - math.lgamma(1.5 * dimension)
    )
    return math.exp(log_value)


def estimate_myrheim_meyer_dimension(
    observed_relation_fraction: float,
    min_dim: float = 1.01,
    max_dim: float = 10.0,
    tol: float = 1e-6,
) -> float:
    """Invert the Myrheim-Meyer relation fraction curve by bisection."""

    observed = float(observed_relation_fraction)
    if not 0.0 < observed < 1.0:
        raise ValueError("observed_relation_fraction must be between 0 and 1")
    lower = float(min_dim)
    upper = float(max_dim)
    if lower <= 0.0 or upper <= lower:
        raise ValueError("dimension bounds must satisfy 0 < min_dim < max_dim")
    if tol <= 0.0:
        raise ValueError("tol must be positive")

    lower_value = myrheim_meyer_relation_fraction(lower)
    upper_value = myrheim_meyer_relation_fraction(upper)
    if not upper_value <= observed <= lower_value:
        raise ValueError(
            "observed relation fraction is outside the requested dimension bounds"
        )

    while upper - lower > tol:
        mid = 0.5 * (lower + upper)
        mid_value = myrheim_meyer_relation_fraction(mid)
        if mid_value > observed:
            lower = mid
        else:
            upper = mid
    return 0.5 * (lower + upper)

