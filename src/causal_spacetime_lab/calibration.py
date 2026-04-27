"""Calibration diagnostics for order-first representation layers."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _as_vector(values: ArrayLike, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def equal_step_deviation(values: ArrayLike) -> float:
    """Return standard deviation of adjacent step sizes."""

    array = _as_vector(values, "values")
    if array.size < 2:
        return float("nan")
    return float(np.std(np.diff(array)))


def ratio_error_under_transform(
    values: ArrayLike,
    transformed: ArrayLike,
    intervals: ArrayLike,
) -> float:
    """Return mean absolute interval-ratio error under a transformation."""

    source = _as_vector(values, "values")
    target = _as_vector(transformed, "transformed")
    if source.shape != target.shape:
        raise ValueError("values and transformed must have equal shape")
    interval_array = np.asarray(intervals, dtype=int)
    if interval_array.ndim != 2 or interval_array.shape[1] != 4:
        raise ValueError("intervals must have shape (m, 4)")
    if np.any((interval_array < 0) | (interval_array >= source.size)):
        raise IndexError("intervals contain an index out of range")

    i = interval_array[:, 0]
    j = interval_array[:, 1]
    k = interval_array[:, 2]
    ell = interval_array[:, 3]
    source_denominator = source[ell] - source[k]
    target_denominator = target[ell] - target[k]
    valid = (np.abs(source_denominator) > 1e-12) & (
        np.abs(target_denominator) > 1e-12
    )
    if not np.any(valid):
        return float("nan")
    source_ratio = (source[j] - source[i]) / source_denominator
    target_ratio = (target[j] - target[i]) / target_denominator
    return float(np.mean(np.abs(source_ratio[valid] - target_ratio[valid])))


def affine_fit_1d(source: ArrayLike, target: ArrayLike) -> tuple[float, float, float]:
    """Fit ``target ~= scale * source + offset`` and return scale, offset, RMSE."""

    x = _as_vector(source, "source")
    y = _as_vector(target, "target")
    if x.shape != y.shape:
        raise ValueError("source and target must have equal shape")
    if x.size == 0:
        raise ValueError("source and target must be nonempty")
    design = np.column_stack((x, np.ones_like(x)))
    scale, offset = np.linalg.lstsq(design, y, rcond=None)[0]
    prediction = scale * x + offset
    rmse = float(np.sqrt(np.mean((prediction - y) ** 2)))
    return float(scale), float(offset), rmse
