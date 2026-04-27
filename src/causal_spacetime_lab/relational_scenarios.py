"""Synthetic persistent-object histories for relational evolution diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _labels_and_positions(
    slice_labels: ArrayLike,
    object_positions: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.asarray(slice_labels, dtype=int)
    positions = np.asarray(object_positions, dtype=float)
    if labels.ndim != 1 or positions.ndim != 1:
        raise ValueError("slice_labels and object_positions must be one-dimensional")
    return labels, positions


def generate_static_configuration_history(
    slice_labels: ArrayLike,
    object_positions: ArrayLike,
) -> dict[int, dict[int, float]]:
    """Generate a persistent static object configuration."""

    labels, positions = _labels_and_positions(slice_labels, object_positions)
    return {
        int(label): {int(i): float(position) for i, position in enumerate(positions)}
        for label in labels
    }


def generate_expanding_configuration_history(
    slice_labels: ArrayLike,
    object_positions: ArrayLike,
    expansion_factors: ArrayLike,
) -> dict[int, dict[int, float]]:
    """Generate a persistent uniformly expanding 1D configuration."""

    labels, positions = _labels_and_positions(slice_labels, object_positions)
    factors = np.asarray(expansion_factors, dtype=float)
    if factors.ndim != 1 or factors.shape[0] != labels.shape[0]:
        raise ValueError("expansion_factors must match slice_labels")
    if np.any(factors <= 0.0):
        raise ValueError("expansion_factors must be positive")
    return {
        int(label): {
            int(i): float(factor * position) for i, position in enumerate(positions)
        }
        for label, factor in zip(labels, factors, strict=True)
    }


def generate_shear_or_reordering_history_1d(
    slice_labels: ArrayLike,
    object_positions: ArrayLike,
    moving_object_id: int,
    displacement_by_slice: ArrayLike,
) -> dict[int, dict[int, float]]:
    """Generate a 1D history where one persistent object moves across slices."""

    labels, positions = _labels_and_positions(slice_labels, object_positions)
    displacements = np.asarray(displacement_by_slice, dtype=float)
    if displacements.ndim != 1 or displacements.shape[0] != labels.shape[0]:
        raise ValueError("displacement_by_slice must match slice_labels")
    moving_id = int(moving_object_id)
    if moving_id < 0 or moving_id >= positions.size:
        raise IndexError("moving_object_id is out of range")
    output: dict[int, dict[int, float]] = {}
    for label, displacement in zip(labels, displacements, strict=True):
        current = positions.copy()
        current[moving_id] += displacement
        output[int(label)] = {
            int(i): float(position) for i, position in enumerate(current)
        }
    return output


def add_position_noise_to_history(
    history: dict[int, dict[int, float]],
    noise_scale: float,
    seed: int | None = None,
) -> dict[int, dict[int, float]]:
    """Add independent position noise while preserving slice and object keys."""

    if noise_scale < 0.0:
        raise ValueError("noise_scale must be nonnegative")
    rng = np.random.default_rng(seed)
    return {
        int(label): {
            int(object_id): float(value + rng.normal(scale=noise_scale))
            for object_id, value in positions.items()
        }
        for label, positions in sorted(history.items())
    }
