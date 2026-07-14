"""Shared descriptive diagnostics for finite 2D-order MCMC runs."""

from __future__ import annotations

import numpy as np


def integrated_autocorrelation(values: list[float]) -> tuple[float, float]:
    """Return initial-positive-sequence IAT and ESS for a scalar chain."""

    array = np.asarray(values, dtype=float)
    centered = array - array.mean()
    variance = float(np.dot(centered, centered) / array.size)
    if variance == 0.0:
        return 1.0, float(array.size)
    tau = 1.0
    for lag in range(1, array.size // 2 + 1):
        covariance = float(np.dot(centered[:-lag], centered[lag:]) / (array.size - lag))
        correlation = covariance / variance
        if correlation <= 0.0:
            break
        tau += 2.0 * correlation
    return tau, float(array.size / tau)


def classify_phase(samples: list[dict], reference: dict) -> str:
    """Apply the frozen P5 profile rules to retained chain samples."""

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
