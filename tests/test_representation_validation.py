from __future__ import annotations

import numpy as np

from causal_spacetime_lab.representation_validation import (
    evaluate_embedding_on_constraints,
    flip_constraint_fraction,
    random_quadruplet_constraints,
    representation_generalization_report,
    shuffle_constraint_sides,
    split_constraints,
)


def test_split_constraints_deterministic_sizes() -> None:
    constraints = np.arange(40, dtype=int).reshape(10, 4)

    train, test = split_constraints(constraints, train_fraction=0.6, seed=1)

    assert train.shape == (6, 4)
    assert test.shape == (4, 4)
    assert {tuple(row) for row in train} | {tuple(row) for row in test} == {
        tuple(row) for row in constraints
    }


def test_shuffle_constraint_sides_preserves_shape_and_values() -> None:
    constraints = np.asarray([[0, 1, 2, 3], [1, 2, 3, 4]], dtype=int)

    shuffled = shuffle_constraint_sides(constraints, seed=2)

    assert shuffled.shape == constraints.shape
    assert sorted(map(tuple, np.sort(shuffled.reshape(-1, 2), axis=1))) == sorted(
        map(tuple, np.sort(constraints.reshape(-1, 2), axis=1))
    )


def test_random_quadruplet_constraints_returns_valid_nonself_pairs() -> None:
    constraints = random_quadruplet_constraints(5, 20, seed=3)

    assert constraints.shape == (20, 4)
    assert np.all(constraints[:, 0] != constraints[:, 1])
    assert np.all(constraints[:, 2] != constraints[:, 3])


def test_flip_constraint_fraction_zero_and_one() -> None:
    constraints = np.asarray([[0, 1, 2, 3], [1, 2, 3, 4]], dtype=int)

    assert np.array_equal(
        flip_constraint_fraction(constraints, 0.0, seed=4),
        constraints,
    )
    flipped = flip_constraint_fraction(constraints, 1.0, seed=4)

    assert np.array_equal(flipped[:, 0:2], constraints[:, 2:4])
    assert np.array_equal(flipped[:, 2:4], constraints[:, 0:2])


def test_evaluate_embedding_on_constraints_returns_finite_fields() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0]])
    constraints = np.asarray([[0, 1, 0, 2]], dtype=int)

    report = evaluate_embedding_on_constraints(coords, constraints)

    assert np.isfinite(report["violation_rate"])
    assert np.isfinite(report["hinge_loss"])
    assert report["constraint_count"] == 1.0


def test_representation_generalization_report_returns_train_test_fields() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0]])
    constraints = np.asarray([[0, 1, 0, 2]], dtype=int)

    report = representation_generalization_report(coords, constraints, constraints)

    assert "train_violation_rate" in report
    assert "test_violation_rate" in report
    assert "generalization_gap" in report
