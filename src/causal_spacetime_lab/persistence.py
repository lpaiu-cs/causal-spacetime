"""Supplied object-persistence utilities for protocol-relative velocity tests."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def generate_persistent_object_events_1p1(
    slice_time_values: ArrayLike,
    initial_positions: ArrayLike,
    velocities: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.int_], NDArray[np.int_]]:
    """Generate synthetic persistent object events in 1+1D."""

    times = np.asarray(slice_time_values, dtype=float)
    initial = np.asarray(initial_positions, dtype=float)
    velocity_values = np.asarray(velocities, dtype=float)
    if times.ndim != 1 or initial.ndim != 1 or velocity_values.ndim != 1:
        raise ValueError("all inputs must be one-dimensional")
    if initial.shape != velocity_values.shape:
        raise ValueError("initial_positions and velocities must have same shape")
    events: list[tuple[float, float]] = []
    slice_ids: list[int] = []
    object_ids: list[int] = []
    t0 = float(times[0]) if times.size else 0.0
    for slice_index, time in enumerate(times):
        for object_id, (x0, velocity) in enumerate(
            zip(initial, velocity_values, strict=True)
        ):
            events.append((float(time), float(x0 + velocity * (time - t0))))
            slice_ids.append(slice_index)
            object_ids.append(object_id)
    return (
        np.asarray(events, dtype=np.float64),
        np.asarray(slice_ids, dtype=int),
        np.asarray(object_ids, dtype=int),
    )


def object_observations_by_id(
    event_indices: ArrayLike,
    object_ids: ArrayLike,
) -> dict[int, NDArray[np.int_]]:
    """Group event indices by supplied persistent object id."""

    indices = np.asarray(event_indices, dtype=int)
    objects = np.asarray(object_ids, dtype=int)
    if indices.ndim != 1 or objects.ndim != 1 or indices.shape != objects.shape:
        raise ValueError("event_indices and object_ids must align")
    return {
        int(object_id): indices[objects == object_id].astype(int, copy=False)
        for object_id in sorted(set(int(value) for value in objects))
    }


def finite_difference_velocity_by_object(
    positions: ArrayLike,
    slice_labels: ArrayLike,
    object_ids: ArrayLike,
    time_values_by_slice: dict[int, float],
) -> dict[int, float]:
    """Estimate mean finite-difference velocity for each persistent object."""

    pos = np.asarray(positions, dtype=float)
    if pos.ndim == 2 and pos.shape[1] == 1:
        pos = pos[:, 0]
    if pos.ndim != 1:
        raise ValueError("positions must be one-dimensional or shape (n, 1)")
    slices = np.asarray(slice_labels, dtype=int)
    objects = np.asarray(object_ids, dtype=int)
    if slices.shape != pos.shape or objects.shape != pos.shape:
        raise ValueError("positions, slice_labels, and object_ids must align")
    result: dict[int, float] = {}
    for object_id in sorted(set(int(value) for value in objects)):
        mask = objects == object_id
        valid = mask & np.isfinite(pos)
        valid &= np.asarray([int(s) in time_values_by_slice for s in slices])
        if np.count_nonzero(valid) < 2:
            continue
        times = np.asarray([time_values_by_slice[int(s)] for s in slices[valid]])
        values = pos[valid]
        order = np.argsort(times)
        times = times[order]
        values = values[order]
        dt = np.diff(times)
        finite = dt != 0.0
        if not np.any(finite):
            continue
        result[object_id] = float(np.mean(np.diff(values)[finite] / dt[finite]))
    return result
