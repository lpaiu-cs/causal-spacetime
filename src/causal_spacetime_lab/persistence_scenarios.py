"""Synthetic unlabeled persistence scenarios."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _emit_unlabeled(
    positions_by_object: np.ndarray,
    rng: np.random.Generator,
    random_per_slice_permutation: bool,
) -> tuple[np.ndarray, dict[int, int]]:
    n = positions_by_object.size
    if random_per_slice_permutation:
        local_to_object = rng.permutation(n)
    else:
        local_to_object = np.arange(n)
    local_positions = positions_by_object[local_to_object]
    track = {
        int(object_id): int(np.flatnonzero(local_to_object == object_id)[0])
        for object_id in range(n)
    }
    return local_positions.astype(float), track


def generate_unlabeled_static_history(
    slice_count: int,
    object_positions: ArrayLike,
    *,
    random_per_slice_permutation: bool = True,
    seed: int | None = None,
) -> tuple[dict[int, np.ndarray], dict[int, dict[int, int]]]:
    """Generate an unlabeled static persistent-object history."""

    positions = np.asarray(object_positions, dtype=float)
    if slice_count <= 0 or positions.ndim != 1:
        raise ValueError("slice_count must be positive and object_positions 1D")
    rng = np.random.default_rng(seed)
    unlabeled: dict[int, np.ndarray] = {}
    tracks: dict[int, dict[int, int]] = {}
    for label in range(slice_count):
        unlabeled[label], tracks[label] = _emit_unlabeled(
            positions,
            rng,
            random_per_slice_permutation,
        )
    return unlabeled, tracks


def generate_unlabeled_small_motion_history(
    slice_count: int,
    object_positions: ArrayLike,
    velocities: ArrayLike,
    *,
    random_per_slice_permutation: bool = True,
    seed: int | None = None,
) -> tuple[dict[int, np.ndarray], dict[int, dict[int, int]]]:
    """Generate an unlabeled small-motion persistent-object history."""

    positions = np.asarray(object_positions, dtype=float)
    velocity_values = np.asarray(velocities, dtype=float)
    if (
        slice_count <= 0
        or positions.ndim != 1
        or velocity_values.shape != positions.shape
    ):
        raise ValueError("slice_count must be positive and inputs must align")
    rng = np.random.default_rng(seed)
    unlabeled: dict[int, np.ndarray] = {}
    tracks: dict[int, dict[int, int]] = {}
    for label in range(slice_count):
        current = positions + velocity_values * float(label)
        unlabeled[label], tracks[label] = _emit_unlabeled(
            current,
            rng,
            random_per_slice_permutation,
        )
    return unlabeled, tracks


def generate_unlabeled_crossing_history_1d(
    slice_count: int,
    object_positions: ArrayLike,
    moving_object_ids: tuple[int, int],
    *,
    random_per_slice_permutation: bool = True,
    seed: int | None = None,
) -> tuple[dict[int, np.ndarray], dict[int, dict[int, int]]]:
    """Generate an unlabeled 1D crossing history."""

    positions = np.asarray(object_positions, dtype=float)
    if slice_count <= 1 or positions.ndim != 1:
        raise ValueError("slice_count > 1 and 1D object_positions are required")
    first, second = int(moving_object_ids[0]), int(moving_object_ids[1])
    if first < 0 or second < 0 or first >= positions.size or second >= positions.size:
        raise IndexError("moving_object_ids are out of range")
    rng = np.random.default_rng(seed)
    start_a, start_b = positions[first], positions[second]
    unlabeled: dict[int, np.ndarray] = {}
    tracks: dict[int, dict[int, int]] = {}
    for label in range(slice_count):
        alpha = label / float(slice_count - 1)
        current = positions.copy()
        current[first] = (1.0 - alpha) * start_a + alpha * start_b
        current[second] = (1.0 - alpha) * start_b + alpha * start_a
        unlabeled[label], tracks[label] = _emit_unlabeled(
            current,
            rng,
            random_per_slice_permutation,
        )
    return unlabeled, tracks


def generate_symmetric_configuration_history(
    slice_count: int,
    object_count: int,
    spacing: float = 1.0,
    *,
    random_per_slice_permutation: bool = True,
    seed: int | None = None,
) -> tuple[dict[int, np.ndarray], dict[int, dict[int, int]]]:
    """Generate an unlabeled symmetric static configuration."""

    if object_count <= 0:
        raise ValueError("object_count must be positive")
    positions = (np.arange(object_count) - (object_count - 1) / 2.0) * float(spacing)
    return generate_unlabeled_static_history(
        slice_count,
        positions,
        random_per_slice_permutation=random_per_slice_permutation,
        seed=seed,
    )
