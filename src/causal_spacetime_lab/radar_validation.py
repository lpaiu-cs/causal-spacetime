"""Validation helpers for discrete radar reconstruction experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def true_stationary_radar_coordinates_1p1(
    events: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return analytic stationary-observer radar coordinates for validation."""

    event_array = np.asarray(events, dtype=float)
    if event_array.ndim != 2 or event_array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")
    return event_array[:, 0].copy(), np.abs(event_array[:, 1])


def radar_time_error(
    reconstructed_time: ArrayLike,
    true_time: ArrayLike,
) -> NDArray[np.float64]:
    """Return signed radar-time error ``reconstructed - true``."""

    reconstructed = np.asarray(reconstructed_time, dtype=float)
    true = np.asarray(true_time, dtype=float)
    if reconstructed.shape != true.shape:
        raise ValueError("reconstructed_time and true_time must have the same shape")
    return reconstructed - true


def radar_distance_error(
    reconstructed_distance: ArrayLike,
    true_distance: ArrayLike,
) -> NDArray[np.float64]:
    """Return signed radar-distance error ``reconstructed - true``."""

    reconstructed = np.asarray(reconstructed_distance, dtype=float)
    true = np.asarray(true_distance, dtype=float)
    if reconstructed.shape != true.shape:
        raise ValueError(
            "reconstructed_distance and true_distance must have the same shape"
        )
    return reconstructed - true


def accessible_fraction(accessible: ArrayLike) -> float:
    """Return the fraction of targets with accessible radar ticks."""

    values = np.asarray(accessible, dtype=bool)
    if values.size == 0:
        raise ValueError("accessible must not be empty")
    return float(np.mean(values))

