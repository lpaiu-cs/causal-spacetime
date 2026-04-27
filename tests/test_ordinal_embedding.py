from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    ordinal_embedding_dimension_curve,
    procrustes_align,
    quadruplet_hinge_loss,
    quadruplet_violation_rate,
    sample_quadruplet_constraints_from_coords,
    squared_distance_matrix,
)


def test_squared_distance_matrix_hand_coded_points() -> None:
    coords = np.asarray([[0.0, 0.0], [3.0, 4.0], [1.0, 0.0]])

    distances = squared_distance_matrix(coords)

    assert distances == pytest.approx(
        np.asarray(
            [
                [0.0, 25.0, 1.0],
                [25.0, 0.0, 20.0],
                [1.0, 20.0, 0.0],
            ]
        )
    )


def test_sample_quadruplet_constraints_are_satisfied_by_true_coords() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0], [7.0]])

    constraints = sample_quadruplet_constraints_from_coords(coords, 20, seed=5)

    assert constraints.shape == (20, 4)
    assert quadruplet_violation_rate(coords, constraints) == pytest.approx(0.0)


def test_quadruplet_hinge_loss_is_finite_and_nonnegative() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0]])
    constraints = np.asarray([[0, 1, 0, 2]], dtype=int)

    loss = quadruplet_hinge_loss(coords, constraints)

    assert np.isfinite(loss)
    assert loss >= 0.0


def test_fit_ordinal_embedding_gradient_descent_tiny_example() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0], [6.0]])
    constraints = sample_quadruplet_constraints_from_coords(coords, 30, seed=7)

    embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
        4,
        1,
        constraints,
        steps=20,
        restarts=1,
        seed=11,
    )

    assert embedding.shape == (4, 1)
    assert np.isfinite(diagnostics["final_loss"])
    assert np.isfinite(diagnostics["violation_rate"])


def test_procrustes_align_recovers_similarity_transform() -> None:
    reference = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, 2.0]])
    theta = np.pi / 4.0
    rotation = np.asarray(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )
    estimated = 3.0 * reference @ rotation + np.asarray([5.0, -2.0])

    _, diagnostics = procrustes_align(estimated, reference)

    assert diagnostics["rmse"] < 1e-12


def test_embedding_distance_order_error_identical_coords_is_zero() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0], [6.0]])

    error = embedding_distance_order_error(coords, coords, seed=13)

    assert error == pytest.approx(0.0)


def test_ordinal_embedding_dimension_curve_returns_candidate_rows() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0], [6.0]])
    constraints = sample_quadruplet_constraints_from_coords(coords, 20, seed=17)

    rows = ordinal_embedding_dimension_curve(
        4,
        constraints,
        [1, 2],
        steps=10,
        restarts=1,
        seed=19,
    )

    assert [row["embedding_dim"] for row in rows] == [1.0, 2.0]
    assert all(np.isfinite(row["final_loss"]) for row in rows)
