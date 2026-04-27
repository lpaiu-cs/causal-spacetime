"""Finite-sampling error models for 1+1D interval-count reconstruction."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.estimators import bias, mae, relative_rmse, rmse
from causal_spacetime_lab.metrics import causal_diamond_volume_1p1

__all__ = [
    "bias",
    "binomial_interval_mean_1p1",
    "binomial_interval_var_1p1",
    "binomial_tau_abs_std_delta_1p1",
    "binomial_tau_rel_std_delta_1p1",
    "binned_summary",
    "bootstrap_mean_ci",
    "interval_volume_1p1",
    "mae",
    "poisson_interval_mean_1p1",
    "poisson_tau_abs_std_delta_1p1",
    "poisson_tau_rel_std_delta_1p1",
    "relative_rmse",
    "rmse",
]


def interval_volume_1p1(tau: float) -> float:
    """Return the 1+1D Alexandrov interval volume ``V = tau**2 / 2``."""

    tau_value = float(tau)
    if tau_value < 0.0:
        raise ValueError("tau must be non-negative")
    return 0.5 * tau_value * tau_value


def poisson_interval_mean_1p1(tau: float, rho: float) -> float:
    """Return Poisson mean interval count for proper time ``tau``."""

    if rho <= 0:
        raise ValueError("rho must be positive")
    return rho * interval_volume_1p1(tau)


def poisson_tau_abs_std_delta_1p1(rho: float) -> float:
    """Delta-method absolute standard deviation for ``sqrt(2K / rho)``.

    For Poisson ``K`` with mean ``rho * tau**2 / 2``, this leading-order
    expression is independent of ``tau`` away from the very-small-interval
    regime.
    """

    if rho <= 0:
        raise ValueError("rho must be positive")
    return 1.0 / math.sqrt(2.0 * rho)


def poisson_tau_rel_std_delta_1p1(tau: float, rho: float) -> float:
    """Delta-method relative standard deviation for the Poisson model."""

    tau_value = float(tau)
    if tau_value <= 0.0:
        raise ValueError("tau must be positive")
    return poisson_tau_abs_std_delta_1p1(rho) / tau_value


def _binomial_interval_probability_1p1(tau: float, T: float) -> float:
    interval_volume = interval_volume_1p1(tau)
    global_volume = causal_diamond_volume_1p1(T)
    probability = interval_volume / global_volume
    if probability > 1.0 + 1e-12:
        raise ValueError("interval volume cannot exceed global diamond volume")
    return min(max(probability, 0.0), 1.0)


def binomial_interval_mean_1p1(tau: float, T: float, n_support: int) -> float:
    """Return fixed-N binomial mean interval count."""

    count = int(n_support)
    if count < 0:
        raise ValueError("n_support must be non-negative")
    return count * _binomial_interval_probability_1p1(tau, T)


def binomial_interval_var_1p1(tau: float, T: float, n_support: int) -> float:
    """Return fixed-N binomial interval-count variance."""

    count = int(n_support)
    if count < 0:
        raise ValueError("n_support must be non-negative")
    probability = _binomial_interval_probability_1p1(tau, T)
    return count * probability * (1.0 - probability)


def binomial_tau_abs_std_delta_1p1(
    tau: float,
    T: float,
    n_support: int,
    rho: float,
) -> float:
    """Delta-method absolute tau standard deviation for fixed-N sprinkling."""

    tau_value = float(tau)
    if tau_value <= 0.0:
        raise ValueError("tau must be positive")
    if rho <= 0:
        raise ValueError("rho must be positive")

    variance = binomial_interval_var_1p1(tau_value, T, n_support)
    derivative = 1.0 / (rho * tau_value)
    return derivative * math.sqrt(variance)


def binomial_tau_rel_std_delta_1p1(
    tau: float,
    T: float,
    n_support: int,
    rho: float,
) -> float:
    """Delta-method relative tau standard deviation for fixed-N sprinkling."""

    tau_value = float(tau)
    if tau_value <= 0.0:
        raise ValueError("tau must be positive")
    return binomial_tau_abs_std_delta_1p1(tau_value, T, n_support, rho) / tau_value


def binned_summary(
    x: ArrayLike,
    y: ArrayLike,
    bins: ArrayLike,
) -> list[dict[str, float]]:
    """Summarize ``y`` values in half-open bins of ``x``.

    The last bin includes its right edge. Empty bins are retained with ``count``
    zero and ``NaN`` summary statistics.
    """

    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    bin_edges = np.asarray(bins, dtype=float)
    if x_values.shape != y_values.shape:
        raise ValueError("x and y must have the same shape")
    if x_values.size == 0:
        raise ValueError("x and y must not be empty")
    if bin_edges.ndim != 1 or bin_edges.size < 2:
        raise ValueError("bins must contain at least two edges")
    if np.any(np.diff(bin_edges) <= 0):
        raise ValueError("bins must be strictly increasing")

    rows: list[dict[str, float]] = []
    for index, (left, right) in enumerate(
        zip(bin_edges[:-1], bin_edges[1:], strict=True)
    ):
        if index == bin_edges.size - 2:
            mask = (x_values >= left) & (x_values <= right)
        else:
            mask = (x_values >= left) & (x_values < right)
        values = y_values[mask]
        count = int(values.size)
        if count == 0:
            mean = float("nan")
            std = float("nan")
            root_mean_square = float("nan")
        else:
            mean = float(np.mean(values))
            std = float(np.std(values, ddof=1)) if count > 1 else 0.0
            root_mean_square = float(np.sqrt(np.mean(values * values)))
        rows.append(
            {
                "bin_left": float(left),
                "bin_right": float(right),
                "bin_center": float(0.5 * (left + right)),
                "count": float(count),
                "mean": mean,
                "std": std,
                "rmse": root_mean_square,
            }
        )
    return rows


def bootstrap_mean_ci(
    values: ArrayLike,
    confidence: float = 0.95,
    n_resamples: int = 1000,
    seed: int | None = None,
) -> tuple[float, float, float]:
    """Return mean and percentile bootstrap confidence interval."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError("values must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(array)):
        raise ValueError("values must be finite")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must be between 0 and 1")
    if n_resamples <= 0:
        raise ValueError("n_resamples must be positive")

    rng = np.random.default_rng(seed)
    sample_indices = rng.integers(0, array.size, size=(int(n_resamples), array.size))
    means = np.mean(array[sample_indices], axis=1)
    alpha = 1.0 - confidence
    lower, upper = np.quantile(means, [alpha / 2.0, 1.0 - alpha / 2.0])
    return float(np.mean(array)), float(lower), float(upper)
