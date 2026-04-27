"""Estimators and error metrics for causal-set reconstruction experiments."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.constants import ETA_1P1_CAUSAL_DIAMOND
from causal_spacetime_lab.geometry import estimate_tau_from_interval_count
from causal_spacetime_lab.metrics import causal_diamond_volume_1p1


def global_density_1p1(n_events: int, T: float) -> float:
    """Return event density in a 1+1D causal diamond of proper duration ``T``.

    The density is supplied scale information for reconstruction experiments,
    not information derived from causal order alone.
    """

    count = int(n_events)
    if count < 0:
        raise ValueError("n_events must be non-negative")
    return count / causal_diamond_volume_1p1(T)


def estimate_tau_from_interval_count_1p1(
    interval_count: int,
    rho: float,
    eta_d: float = ETA_1P1_CAUSAL_DIAMOND,
) -> float:
    """Estimate timelike proper time from Alexandrov interval cardinality.

    ``interval_count`` is expected to count sampled events strictly between the
    two endpoint events. ``rho`` supplies event density as additional metric
    scale information.
    """

    count = int(interval_count)
    if count < 0:
        raise ValueError("interval_count must be non-negative")
    if rho <= 0:
        raise ValueError("rho must be positive")
    if eta_d <= 0:
        raise ValueError("eta_d must be positive")

    if eta_d == ETA_1P1_CAUSAL_DIAMOND:
        return estimate_tau_from_interval_count(count, rho, spacetime_dim=2)

    volume_hat = count / rho
    return math.sqrt(volume_hat / eta_d)


def estimate_tau_from_longest_chain_1p1(
    chain_length: int,
    rho: float,
    chain_counts_endpoints: bool = True,
) -> float:
    """Estimate proper time from a longest-chain length in 1+1D.

    This uses the asymptotic relation ``L ~ sqrt(2 * rho) * tau``. When
    ``chain_counts_endpoints`` is true, ``chain_length`` is interpreted as the
    convention used by :func:`causal_spacetime_lab.chains.longest_chain_length`:
    both requested endpoint events are included in the chain. The estimator then
    subtracts two endpoint contributions before applying the asymptotic
    normalization. This correction is simple and finite-size sensitive, so
    chain estimates should be treated as diagnostic rather than calibrated.
    """

    length = int(chain_length)
    if length < 0:
        raise ValueError("chain_length must be non-negative")
    if rho <= 0:
        raise ValueError("rho must be positive")

    effective_length = max(length - 2, 0) if chain_counts_endpoints else length
    return effective_length / math.sqrt(2.0 * rho)


def _paired_arrays(
    y_true: ArrayLike,
    y_pred: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    true = np.asarray(y_true, dtype=float)
    pred = np.asarray(y_pred, dtype=float)
    if true.shape != pred.shape:
        raise ValueError("y_true and y_pred must have the same shape")
    if true.size == 0:
        raise ValueError("metric inputs must not be empty")
    if not np.all(np.isfinite(true)) or not np.all(np.isfinite(pred)):
        raise ValueError("metric inputs must be finite")
    return true, pred


def mae(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """Return mean absolute error."""

    true, pred = _paired_arrays(y_true, y_pred)
    return float(np.mean(np.abs(pred - true)))


def rmse(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """Return root mean squared error."""

    true, pred = _paired_arrays(y_true, y_pred)
    return float(np.sqrt(np.mean((pred - true) ** 2)))


def bias(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """Return mean signed error ``prediction - truth``."""

    true, pred = _paired_arrays(y_true, y_pred)
    return float(np.mean(pred - true))


def relative_rmse(y_true: ArrayLike, y_pred: ArrayLike) -> float:
    """Return root mean squared relative error.

    The denominator is the true value for each sample, so zero true values are
    not accepted.
    """

    true, pred = _paired_arrays(y_true, y_pred)
    if np.any(true == 0.0):
        raise ValueError("relative_rmse requires nonzero true values")
    return float(np.sqrt(np.mean(((pred - true) / true) ** 2)))


def error_metric_summary(y_true: ArrayLike, y_pred: ArrayLike) -> dict[str, float]:
    """Return standard reconstruction error metrics."""

    return {
        "mae": mae(y_true, y_pred),
        "rmse": rmse(y_true, y_pred),
        "bias": bias(y_true, y_pred),
        "relative_rmse": relative_rmse(y_true, y_pred),
    }
