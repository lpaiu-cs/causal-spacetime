"""Validation utilities for effective metric-representation diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.ordinal_embedding import (
    quadruplet_hinge_loss,
    quadruplet_violation_rate,
)


def _as_constraints(constraints: ArrayLike) -> NDArray[np.int_]:
    array = np.asarray(constraints, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    if np.any(array < 0):
        raise IndexError("constraints contain negative indices")
    return array


def split_constraints(
    constraints: ArrayLike,
    train_fraction: float = 0.7,
    seed: int | None = None,
) -> tuple[NDArray[np.int_], NDArray[np.int_]]:
    """Split quadruplet constraints into deterministic train/test subsets."""

    constraint_array = _as_constraints(constraints)
    fraction = float(train_fraction)
    if not 0.0 < fraction < 1.0:
        raise ValueError("train_fraction must satisfy 0 < train_fraction < 1")
    rng = np.random.default_rng(seed)
    order = rng.permutation(constraint_array.shape[0])
    train_count = int(round(fraction * constraint_array.shape[0]))
    train_count = min(max(train_count, 1), max(1, constraint_array.shape[0] - 1))
    train = constraint_array[order[:train_count]]
    test = constraint_array[order[train_count:]]
    return train.astype(int, copy=False), test.astype(int, copy=False)


def shuffle_constraint_sides(
    constraints: ArrayLike,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Randomly swap left/right pair in each constraint with probability 0.5."""

    shuffled = _as_constraints(constraints).copy()
    rng = np.random.default_rng(seed)
    swap = rng.random(shuffled.shape[0]) < 0.5
    left = shuffled[swap, 0:2].copy()
    shuffled[swap, 0:2] = shuffled[swap, 2:4]
    shuffled[swap, 2:4] = left
    return shuffled


def random_quadruplet_constraints(
    n_points: int,
    num_constraints: int,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Return unconstrained random quadruplets with valid non-self pairs."""

    n = int(n_points)
    count = int(num_constraints)
    if n < 2:
        raise ValueError("n_points must be at least 2")
    if count < 0:
        raise ValueError("num_constraints must be nonnegative")
    rng = np.random.default_rng(seed)
    rows = np.empty((count, 4), dtype=int)
    accepted = 0
    while accepted < count:
        batch = max(256, 2 * (count - accepted))
        candidate = rng.integers(0, n, size=(batch, 4))
        left_ok = candidate[:, 0] != candidate[:, 1]
        right_ok = candidate[:, 2] != candidate[:, 3]
        left_pair = np.sort(candidate[:, 0:2], axis=1)
        right_pair = np.sort(candidate[:, 2:4], axis=1)
        different_pairs = np.any(left_pair != right_pair, axis=1)
        valid = candidate[left_ok & right_ok & different_pairs]
        take = min(valid.shape[0], count - accepted)
        if take:
            rows[accepted : accepted + take] = valid[:take]
            accepted += take
    return rows


def flip_constraint_fraction(
    constraints: ArrayLike,
    flip_fraction: float,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Swap left/right pairs for a controlled fraction of constraints."""

    flipped = _as_constraints(constraints).copy()
    fraction = float(flip_fraction)
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("flip_fraction must satisfy 0 <= flip_fraction <= 1")
    count = int(round(fraction * flipped.shape[0]))
    if count == 0:
        return flipped
    rng = np.random.default_rng(seed)
    selected = rng.choice(flipped.shape[0], size=count, replace=False)
    left = flipped[selected, 0:2].copy()
    flipped[selected, 0:2] = flipped[selected, 2:4]
    flipped[selected, 2:4] = left
    return flipped


def evaluate_embedding_on_constraints(
    coords: ArrayLike,
    constraints: ArrayLike,
    margin: float = 1e-3,
) -> dict[str, float]:
    """Evaluate coordinates against quadruplet order constraints."""

    constraint_array = _as_constraints(constraints)
    return {
        "violation_rate": quadruplet_violation_rate(coords, constraint_array),
        "hinge_loss": quadruplet_hinge_loss(coords, constraint_array, margin=margin),
        "constraint_count": float(constraint_array.shape[0]),
    }


def representation_generalization_report(
    coords: ArrayLike,
    train_constraints: ArrayLike,
    test_constraints: ArrayLike,
    margin: float = 1e-3,
) -> dict[str, float]:
    """Report train/test order-validation diagnostics for one embedding."""

    train = evaluate_embedding_on_constraints(coords, train_constraints, margin=margin)
    test = evaluate_embedding_on_constraints(coords, test_constraints, margin=margin)
    return {
        "train_violation_rate": train["violation_rate"],
        "test_violation_rate": test["violation_rate"],
        "train_hinge_loss": train["hinge_loss"],
        "test_hinge_loss": test["hinge_loss"],
        "generalization_gap": test["violation_rate"] - train["violation_rate"],
        "train_constraint_count": train["constraint_count"],
        "test_constraint_count": test["constraint_count"],
    }
