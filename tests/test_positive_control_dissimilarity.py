"""Exact tests for dissimilarity, margin, and leak-free constraint splits."""

from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.positive_control.dissimilarity import (
    build_constraint_split,
    flip_constraint_sides,
    margin_from_probe_quantile,
    pair_split_matrix,
    parallax_profiles,
    profile_dissimilarity_matrix,
)
from causal_spacetime_lab.positive_control.echo_profiles import EchoProfileMatrix


def _profiles(delays: np.ndarray, reachable: np.ndarray) -> EchoProfileMatrix:
    return EchoProfileMatrix(
        delay_ranks=delays.astype(float),
        reachable=reachable.astype(bool),
        target_indices=np.arange(delays.shape[0]),
    )


def test_raw_rms_dissimilarity_exact_values() -> None:
    delays = np.array([[2.0, 4.0], [4.0, 8.0], [6.0, 2.0]])
    reachable = np.ones_like(delays, dtype=bool)
    matrix = profile_dissimilarity_matrix(
        _profiles(delays, reachable), 2, center_references=False
    )
    assert matrix[0, 1] == pytest.approx(np.sqrt(10.0))
    assert matrix[0, 2] == pytest.approx(np.sqrt(10.0))
    assert matrix[1, 2] == pytest.approx(np.sqrt(20.0))
    assert np.allclose(matrix, matrix.T, equal_nan=True)
    assert np.all(np.diag(matrix) == 0.0)


def test_centered_parallax_dissimilarity_exact_values() -> None:
    # Rows centered across columns: [-1,1], [-2,2], [2,-2].
    delays = np.array([[2.0, 4.0], [4.0, 8.0], [6.0, 2.0]])
    reachable = np.ones_like(delays, dtype=bool)
    matrix = profile_dissimilarity_matrix(_profiles(delays, reachable), 2)
    assert matrix[0, 1] == pytest.approx(1.0)
    assert matrix[0, 2] == pytest.approx(3.0)
    assert matrix[1, 2] == pytest.approx(4.0)


def test_shared_scalar_profile_carries_no_structure() -> None:
    # Every reference agrees per target (a pure shared scalar / common mode):
    # centering collapses each row to zero, so no distance structure exists.
    shared = np.array([[3.0, 3.0, 3.0, 3.0], [7.0, 7.0, 7.0, 7.0]])
    reachable = np.ones_like(shared, dtype=bool)
    centered = parallax_profiles(_profiles(shared, reachable))
    assert np.allclose(centered, 0.0)
    matrix = profile_dissimilarity_matrix(_profiles(shared, reachable), 4)
    off_diagonal = matrix[~np.eye(2, dtype=bool)]
    assert np.allclose(off_diagonal, 0.0)


def test_dissimilarity_undefined_below_min_common_columns() -> None:
    delays = np.array([[1.0, 2.0], [3.0, 4.0]])
    reachable = np.array([[True, False], [False, True]])
    matrix = profile_dissimilarity_matrix(_profiles(delays, reachable), 1)
    assert np.isnan(matrix[0, 1])


def test_margin_probe_is_deterministic_and_positive() -> None:
    rng = np.random.default_rng(0)
    coords = rng.uniform(size=20)
    matrix = np.abs(coords[:, None] - coords[None, :])
    first = margin_from_probe_quantile(matrix, seed=4)
    second = margin_from_probe_quantile(matrix, seed=4)
    assert first == second
    assert first > 0.0


def test_pair_split_is_symmetric_deterministic_and_balanced() -> None:
    split = pair_split_matrix(200, train_fraction=0.8, seed=9)
    assert np.array_equal(split, split.T)
    assert not np.any(np.diag(split))
    assert np.array_equal(split, pair_split_matrix(200, 0.8, seed=9))
    upper = split[np.triu_indices(200, k=1)]
    assert 0.75 < float(np.mean(upper)) < 0.85


def test_constraint_split_has_no_pair_leakage_and_respects_margin() -> None:
    rng = np.random.default_rng(1)
    coords = rng.uniform(size=30)
    matrix = np.abs(coords[:, None] - coords[None, :])
    margin = margin_from_probe_quantile(matrix, seed=2)
    split = build_constraint_split(matrix, 300, 80, margin, 0.8, seed=7)
    train_pairs = pair_split_matrix(30, 0.8, seed=7)

    for rows, expected_side in ((split.train, True), (split.heldout, False)):
        assert rows.shape[1] == 4
        left = train_pairs[rows[:, 0], rows[:, 1]]
        right = train_pairs[rows[:, 2], rows[:, 3]]
        assert bool(np.all(left == expected_side))
        assert bool(np.all(right == expected_side))
        dij = matrix[rows[:, 0], rows[:, 1]]
        dkl = matrix[rows[:, 2], rows[:, 3]]
        assert bool(np.all(dij + margin < dkl))


def test_flip_control_swaps_sides() -> None:
    constraints = np.array([[0, 1, 2, 3], [4, 5, 6, 7]])
    flipped = flip_constraint_sides(constraints, flip_probability=1.0, seed=0)
    assert np.array_equal(flipped, constraints[:, [2, 3, 0, 1]])
