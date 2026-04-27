from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.metric_representation import (
    euclidean_embedding_diagnostics,
    has_order_cycle,
    topological_rank_representation,
    triangle_inequality_violations,
    unordered_pair_indices,
)


def test_unordered_pair_indices() -> None:
    pairs = unordered_pair_indices(4)

    assert pairs.tolist() == [
        [0, 1],
        [0, 2],
        [0, 3],
        [1, 2],
        [1, 3],
        [2, 3],
    ]


def test_has_order_cycle_positive_and_negative_examples() -> None:
    assert not has_order_cycle(3, np.asarray([[0, 1], [1, 2]], dtype=int))
    assert has_order_cycle(3, np.asarray([[0, 1], [1, 2], [2, 0]], dtype=int))


def test_topological_rank_representation_satisfies_constraints() -> None:
    constraints = np.asarray([[0, 1], [1, 2], [0, 3]], dtype=int)

    ranks = topological_rank_representation(4, constraints)

    for source, target in constraints:
        assert ranks[source] < ranks[target]


def test_triangle_inequality_violations_metric_and_nonmetric() -> None:
    assert triangle_inequality_violations(np.asarray([1.0, 1.0, 1.0]), 3) == 0
    assert triangle_inequality_violations(np.asarray([1.0, 1.0, 3.0]), 3) == 1


def test_euclidean_embedding_diagnostics_simple_distance_matrix() -> None:
    matrix = np.asarray([[0.0, 1.0], [1.0, 0.0]])

    diagnostics = euclidean_embedding_diagnostics(matrix)

    assert diagnostics["is_symmetric"] is True
    assert diagnostics["has_zero_diagonal"] is True
    assert diagnostics["is_euclidean_candidate"] is True
    assert diagnostics["min_gram_eigenvalue"] == pytest.approx(0.0)
