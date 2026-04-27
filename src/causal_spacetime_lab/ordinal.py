"""Ordinal comparison utilities for order-first reconstruction diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_value_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def order_matrix_from_values(
    values: ArrayLike,
    strict: bool = True,
) -> NDArray[np.bool_]:
    """Return the pairwise order matrix induced by scalar values."""

    array = _as_value_vector(values, "values")
    if strict:
        return array[:, None] < array[None, :]
    return array[:, None] <= array[None, :]


def _pair_order_signs(
    values: NDArray[np.float64],
    tolerance: float,
) -> NDArray[np.float64]:
    diffs = values[None, :] - values[:, None]
    signs = np.zeros_like(diffs, dtype=np.float64)
    signs[diffs > tolerance] = 1.0
    signs[diffs < -tolerance] = -1.0
    iu, ju = np.triu_indices(values.size, k=1)
    return signs[iu, ju]


def order_inversion_rate(
    true_values: ArrayLike,
    estimated_values: ArrayLike,
    *,
    ignore_ties: bool = True,
    tolerance: float = 0.0,
) -> float:
    """Return the fraction of comparable scalar pairs whose order is not preserved.

    When true ties are ignored, any estimated reversal or finite-resolution tie
    for a truly ordered pair is counted as an ordinal failure. This keeps the
    diagnostic useful for coarse order reconstructions where ties represent lost
    order information.
    """

    truth = _as_value_vector(true_values, "true_values")
    estimate = _as_value_vector(estimated_values, "estimated_values")
    if truth.shape != estimate.shape:
        raise ValueError("true_values and estimated_values must have equal shape")
    if truth.size < 2:
        return float("nan")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")

    true_signs = _pair_order_signs(truth, tol)
    estimated_signs = _pair_order_signs(estimate, tol)
    if ignore_ties:
        comparable = true_signs != 0.0
    else:
        comparable = np.ones(true_signs.shape, dtype=bool)
    if not np.any(comparable):
        return float("nan")
    inverted = true_signs[comparable] != estimated_signs[comparable]
    return float(np.count_nonzero(inverted) / np.count_nonzero(comparable))


def order_agreement_rate(
    true_values: ArrayLike,
    estimated_values: ArrayLike,
    *,
    ignore_ties: bool = True,
    tolerance: float = 0.0,
) -> float:
    """Return one minus the order inversion rate on comparable scalar pairs."""

    inversion = order_inversion_rate(
        true_values,
        estimated_values,
        ignore_ties=ignore_ties,
        tolerance=tolerance,
    )
    return float("nan") if not np.isfinite(inversion) else 1.0 - inversion


def pair_distance_values_1d(
    coords: ArrayLike,
    pairs: ArrayLike,
) -> NDArray[np.float64]:
    """Return absolute 1D pair separations for index pairs."""

    coordinate_values = np.asarray(coords, dtype=float)
    if coordinate_values.ndim == 2 and coordinate_values.shape[1] == 1:
        coordinate_values = coordinate_values[:, 0]
    if coordinate_values.ndim != 1:
        raise ValueError("coords must be one-dimensional or have shape (n, 1)")
    pair_array = np.asarray(pairs, dtype=int)
    if pair_array.ndim != 2 or pair_array.shape[1] != 2:
        raise ValueError("pairs must have shape (m, 2)")
    if np.any((pair_array < 0) | (pair_array >= coordinate_values.size)):
        raise IndexError("pairs contain an index out of range")
    return np.abs(
        coordinate_values[pair_array[:, 0]] - coordinate_values[pair_array[:, 1]]
    ).astype(np.float64, copy=False)


def pair_distance_order_inversion_rate(
    true_coords: ArrayLike,
    estimated_coords: ArrayLike,
    pairs: ArrayLike,
    *,
    tolerance: float = 0.0,
) -> float:
    """Compare ordinal pair-distance order between true and estimated 1D coords."""

    true_distances = pair_distance_values_1d(true_coords, pairs)
    estimated_distances = pair_distance_values_1d(estimated_coords, pairs)
    return order_inversion_rate(
        true_distances,
        estimated_distances,
        ignore_ties=True,
        tolerance=tolerance,
    )
