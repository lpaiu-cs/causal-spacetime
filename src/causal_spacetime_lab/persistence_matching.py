"""Finite persistence matching utilities for unlabeled slice histories."""

from __future__ import annotations

from itertools import permutations
from math import factorial

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.relational_evolution import (
    pair_distance_order_signature,
    signature_order_disagreement,
)


def all_permutations(n: int, max_n: int = 8) -> NDArray[np.int_]:
    """Return all permutations of ``range(n)`` for small finite matching tests."""

    if n < 0:
        raise ValueError("n must be nonnegative")
    if n > max_n:
        raise ValueError(f"n={n} exceeds max_n={max_n}; refusing factorial search")
    if n == 0:
        return np.empty((1, 0), dtype=int)
    return np.asarray(list(permutations(range(n))), dtype=int).reshape(factorial(n), n)


def apply_permutation_to_positions(
    positions: ArrayLike,
    permutation: ArrayLike,
) -> NDArray[np.float64]:
    """Return current positions ordered by previous object ids.

    ``permutation[i]`` is the current local index matched to previous object id
    ``i``.
    """

    values = np.asarray(positions, dtype=float)
    perm = np.asarray(permutation, dtype=int)
    if values.ndim != 1 or perm.ndim != 1:
        raise ValueError("positions and permutation must be one-dimensional")
    if perm.size != values.size:
        raise ValueError("permutation length must match positions length")
    if sorted(perm.tolist()) != list(range(values.size)):
        raise ValueError("permutation must contain each local index exactly once")
    return values[perm].astype(np.float64)


def signature_disagreement_for_permutation(
    previous_positions: ArrayLike,
    current_positions: ArrayLike,
    permutation: ArrayLike,
    tolerance: float = 0.0,
) -> float:
    """Compute relational signature disagreement for one persistence hypothesis."""

    previous = np.asarray(previous_positions, dtype=float)
    current_ordered = apply_permutation_to_positions(current_positions, permutation)
    if previous.ndim != 1 or previous.shape != current_ordered.shape:
        raise ValueError("previous and current positions must have equal 1D shape")
    previous_signature = pair_distance_order_signature(
        {int(i): float(value) for i, value in enumerate(previous)},
        tolerance=tolerance,
    )
    current_signature = pair_distance_order_signature(
        {int(i): float(value) for i, value in enumerate(current_ordered)},
        tolerance=tolerance,
    )
    return signature_order_disagreement(previous_signature, current_signature)


def best_persistence_matchings_by_relational_order(
    previous_positions: ArrayLike,
    current_positions: ArrayLike,
    *,
    max_n: int = 8,
    tolerance: float = 0.0,
    top_k: int = 5,
) -> list[dict[str, object]]:
    """Rank persistence hypotheses by slice-local relational-order continuity."""

    previous = np.asarray(previous_positions, dtype=float)
    current = np.asarray(current_positions, dtype=float)
    if previous.ndim != 1 or current.ndim != 1 or previous.shape != current.shape:
        raise ValueError("previous_positions and current_positions must align")
    rows: list[dict[str, object]] = []
    for perm in all_permutations(previous.size, max_n=max_n):
        cost = signature_disagreement_for_permutation(
            previous,
            current,
            perm,
            tolerance=tolerance,
        )
        rows.append({"permutation": perm.copy(), "cost": float(cost)})
    rows.sort(key=lambda row: (float(row["cost"]), tuple(row["permutation"])))
    for rank, row in enumerate(rows[:top_k], start=1):
        row["rank"] = float(rank)
    return rows[:top_k]


def matching_ambiguity_gap(matches: list[dict[str, object]]) -> float:
    """Return second-best minus best matching cost."""

    if len(matches) < 2:
        return float("nan")
    return float(matches[1]["cost"]) - float(matches[0]["cost"])


def matching_accuracy(
    inferred_permutation: ArrayLike,
    true_permutation: ArrayLike,
) -> float:
    """Return fraction of object matches equal to the validation permutation."""

    inferred = np.asarray(inferred_permutation, dtype=int)
    truth = np.asarray(true_permutation, dtype=int)
    if inferred.shape != truth.shape:
        raise ValueError("permutations must have the same shape")
    if inferred.size == 0:
        return 1.0
    return float(np.mean(inferred == truth))


def compose_permutations(
    first: ArrayLike,
    second: ArrayLike,
) -> NDArray[np.int_]:
    """Compose adjacent matchings: apply ``first`` then ``second``."""

    first_perm = np.asarray(first, dtype=int)
    second_perm = np.asarray(second, dtype=int)
    if first_perm.ndim != 1 or second_perm.ndim != 1:
        raise ValueError("permutations must be one-dimensional")
    if first_perm.size != second_perm.size:
        raise ValueError("permutations must have equal length")
    return second_perm[first_perm].astype(int, copy=False)
