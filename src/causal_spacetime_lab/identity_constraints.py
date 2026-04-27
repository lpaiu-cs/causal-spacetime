"""Partial-label and anchor constraints for persistence matching."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.persistence_matching import (
    all_permutations,
    signature_disagreement_for_permutation,
)


def filter_permutations_by_fixed_points(
    permutations: ArrayLike,
    fixed_matches: dict[int, int],
) -> NDArray[np.int_]:
    """Filter permutations by fixed previous-id -> current-index matches."""

    perms = np.asarray(permutations, dtype=int)
    if perms.ndim != 2:
        raise ValueError("permutations must have shape (m, n)")
    keep = np.ones(perms.shape[0], dtype=bool)
    for previous_id, current_index in fixed_matches.items():
        if previous_id < 0 or previous_id >= perms.shape[1]:
            raise IndexError("fixed previous object id is out of range")
        keep &= perms[:, int(previous_id)] == int(current_index)
    return perms[keep].astype(int, copy=False)


def best_persistence_matchings_with_fixed_points(
    previous_positions: ArrayLike,
    current_positions: ArrayLike,
    fixed_matches: dict[int, int],
    *,
    max_n: int = 8,
    tolerance: float = 0.0,
    top_k: int = 5,
) -> list[dict[str, object]]:
    """Rank persistence matchings subject to supplied fixed-point constraints."""

    previous = np.asarray(previous_positions, dtype=float)
    current = np.asarray(current_positions, dtype=float)
    if previous.ndim != 1 or current.ndim != 1 or previous.shape != current.shape:
        raise ValueError("previous_positions and current_positions must align")
    perms = filter_permutations_by_fixed_points(
        all_permutations(previous.size, max_n=max_n),
        fixed_matches,
    )
    rows: list[dict[str, object]] = []
    for perm in perms:
        rows.append(
            {
                "permutation": perm.copy(),
                "cost": float(
                    signature_disagreement_for_permutation(
                        previous,
                        current,
                        perm,
                        tolerance=tolerance,
                    )
                ),
            }
        )
    rows.sort(key=lambda row: (float(row["cost"]), tuple(row["permutation"])))
    for rank, row in enumerate(rows[:top_k], start=1):
        row["rank"] = float(rank)
    return rows[:top_k]


def generate_partial_identity_constraints(
    true_permutation: ArrayLike,
    known_fraction: float,
    seed: int | None = None,
) -> dict[int, int]:
    """Select a controlled fraction of true fixed identity matches."""

    truth = np.asarray(true_permutation, dtype=int)
    if truth.ndim != 1:
        raise ValueError("true_permutation must be one-dimensional")
    if not 0.0 <= known_fraction <= 1.0:
        raise ValueError("known_fraction must be in [0, 1]")
    count = int(round(float(known_fraction) * truth.size))
    if count == 0:
        return {}
    rng = np.random.default_rng(seed)
    previous_ids = np.sort(rng.choice(truth.size, size=count, replace=False))
    return {int(previous_id): int(truth[previous_id]) for previous_id in previous_ids}
