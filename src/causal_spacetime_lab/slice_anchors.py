"""Supplied anchor structures for cross-slice transport diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def generate_stationary_anchor_events_1p1(
    slice_time_values: ArrayLike,
    anchor_positions: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.int_], NDArray[np.int_]]:
    """Generate supplied persistent stationary anchor events in 1+1D."""

    times = np.asarray(slice_time_values, dtype=float)
    positions = np.asarray(anchor_positions, dtype=float)
    if times.ndim != 1 or positions.ndim != 1:
        raise ValueError("slice_time_values and anchor_positions must be 1D")
    events: list[tuple[float, float]] = []
    slice_ids: list[int] = []
    anchor_ids: list[int] = []
    for slice_index, time in enumerate(times):
        for anchor_id, position in enumerate(positions):
            events.append((float(time), float(position)))
            slice_ids.append(slice_index)
            anchor_ids.append(anchor_id)
    return (
        np.asarray(events, dtype=np.float64),
        np.asarray(slice_ids, dtype=int),
        np.asarray(anchor_ids, dtype=int),
    )


def anchor_indices_by_slice(
    anchor_event_indices: ArrayLike,
    anchor_slice_ids: ArrayLike,
) -> dict[int, NDArray[np.int_]]:
    """Group anchor event indices by slice id."""

    indices = np.asarray(anchor_event_indices, dtype=int)
    slices = np.asarray(anchor_slice_ids, dtype=int)
    if indices.ndim != 1 or slices.ndim != 1 or indices.shape != slices.shape:
        raise ValueError("anchor_event_indices and anchor_slice_ids must align")
    return {
        int(slice_id): indices[slices == slice_id].astype(int, copy=False)
        for slice_id in sorted(set(int(value) for value in slices))
    }


def anchor_reference_positions_by_slice(
    anchor_ids: ArrayLike,
    anchor_positions: ArrayLike,
    anchor_slice_ids: ArrayLike,
) -> dict[int, NDArray[np.float64]]:
    """Return anchor reference positions grouped by slice id."""

    ids = np.asarray(anchor_ids, dtype=int)
    positions = np.asarray(anchor_positions, dtype=float)
    slices = np.asarray(anchor_slice_ids, dtype=int)
    if ids.ndim != 1 or slices.ndim != 1 or ids.shape != slices.shape:
        raise ValueError("anchor_ids and anchor_slice_ids must align")
    if positions.ndim != 1:
        raise ValueError("anchor_positions must be one-dimensional")
    if np.any((ids < 0) | (ids >= positions.size)):
        raise IndexError("anchor_ids contain an index out of range")
    return {
        int(slice_id): positions[ids[slices == slice_id]].astype(np.float64)
        for slice_id in sorted(set(int(value) for value in slices))
    }
