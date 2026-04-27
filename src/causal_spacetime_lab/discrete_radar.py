"""Discrete observer-chain radar reconstruction from causal order."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.observer import validate_observer_chain_clock_order


@dataclass(frozen=True)
class DiscreteRadarCoordinates:
    """Radar coordinates reconstructed from order relations and clock labels."""

    target_index: int
    accessible: bool
    tau_minus: float | None
    tau_plus: float | None
    radar_time: float | None
    radar_distance: float | None


def _as_square_bool_matrix(matrix: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(matrix, dtype=bool)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("causal_matrix must be a square boolean matrix")
    return array


def _validate_indices(indices: ArrayLike, n: int, name: str) -> NDArray[np.int_]:
    array = np.asarray(indices, dtype=int)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if np.any((array < 0) | (array >= n)):
        raise IndexError(f"{name} contains an index out of range")
    return array


def find_radar_ticks_from_order(
    causal_matrix: ArrayLike,
    observer_indices: ArrayLike,
    target_index: int,
    clock_times: ArrayLike,
) -> tuple[float, float] | None:
    """Find ``tau_minus`` and ``tau_plus`` using only order and observer clocks."""

    causal_order = _as_square_bool_matrix(causal_matrix)
    n = causal_order.shape[0]
    observers = _validate_indices(observer_indices, n, "observer_indices")
    target = int(target_index)
    if target < 0 or target >= n:
        raise IndexError("target_index out of range")
    clocks = validate_observer_chain_clock_order(clock_times)
    if clocks.shape != observers.shape:
        raise ValueError("clock_times and observer_indices must have the same length")

    predecessor_mask = causal_order[observers, target]
    successor_mask = causal_order[target, observers]
    if not np.any(predecessor_mask) or not np.any(successor_mask):
        return None

    tau_minus = float(np.max(clocks[predecessor_mask]))
    tau_plus = float(np.min(clocks[successor_mask]))
    return tau_minus, tau_plus


def discrete_radar_coordinates_from_order(
    causal_matrix: ArrayLike,
    observer_indices: ArrayLike,
    target_indices: ArrayLike,
    clock_times: ArrayLike,
) -> list[DiscreteRadarCoordinates]:
    """Reconstruct radar coordinates for targets from causal order and clocks."""

    causal_order = _as_square_bool_matrix(causal_matrix)
    targets = _validate_indices(target_indices, causal_order.shape[0], "target_indices")
    results: list[DiscreteRadarCoordinates] = []

    for target in targets:
        ticks = find_radar_ticks_from_order(
            causal_order,
            observer_indices,
            int(target),
            clock_times,
        )
        if ticks is None:
            results.append(
                DiscreteRadarCoordinates(
                    target_index=int(target),
                    accessible=False,
                    tau_minus=None,
                    tau_plus=None,
                    radar_time=None,
                    radar_distance=None,
                )
            )
            continue

        tau_minus, tau_plus = ticks
        radar_time = 0.5 * (tau_plus + tau_minus)
        radar_distance = 0.5 * (tau_plus - tau_minus)
        results.append(
            DiscreteRadarCoordinates(
                target_index=int(target),
                accessible=True,
                tau_minus=tau_minus,
                tau_plus=tau_plus,
                radar_time=radar_time,
                radar_distance=radar_distance,
            )
        )

    return results

