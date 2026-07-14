"""Shared descriptive diagnostics for finite 2D-order MCMC runs."""

from __future__ import annotations

import numpy as np

IAT_ESTIMATOR = "Geyer initial-positive-pair sequence"


def _initial_positive_pair_tau(correlations: list[float]) -> float:
    """Return Geyer's initial-positive-pair IAT, floored at one."""

    pair_sum = 0.0
    for start in range(0, len(correlations) - 1, 2):
        gamma = correlations[start] + correlations[start + 1]
        if gamma <= 0.0:
            break
        pair_sum += gamma
    return max(1.0, -1.0 + 2.0 * pair_sum)


def integrated_autocorrelation(values: list[float]) -> tuple[float, float]:
    """Return Geyer initial-positive-pair IAT and ESS for a scalar chain."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError("values must be a non-empty finite one-dimensional list")
    centered = array - array.mean()
    variance = float(np.dot(centered, centered) / array.size)
    if variance == 0.0:
        return 1.0, float(array.size)
    correlations = [1.0]
    for lag in range(1, array.size):
        covariance = float(np.dot(centered[:-lag], centered[lag:]) / array.size)
        correlations.append(covariance / variance)
    tau = _initial_positive_pair_tau(correlations)
    return tau, float(array.size / tau)


def classify_phase(samples: list[dict], reference: dict) -> str:
    """Apply the P4/P5-derived phase profile rules frozen for the caller."""

    mean_n0 = float(np.mean([row["n0"] for row in samples]))
    mean_n12 = float(np.mean([row["n1"] + row["n2"] for row in samples]))
    mean_height = float(np.mean([row["height"] for row in samples]))
    reference_n12 = reference["n1"] + reference["n2"]
    if (
        abs(mean_n0 / reference["n0"] - 1.0) <= 0.5
        and mean_height >= 0.7 * reference["height"]
    ):
        return "continuum"
    if mean_n12 <= 0.5 * reference_n12 and mean_height <= 0.25 * reference["height"]:
        return "crystal"
    return "intermediate"
