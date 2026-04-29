from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_response_representability import (
    has_response_order_cycle,
    response_order_directed_edges,
    response_order_topological_ranks,
    scalar_rank_representation_error,
    scalar_representability_report,
)


def test_response_order_directed_edges() -> None:
    signs = np.asarray([[0, -1, 1], [1, 0, 0], [-1, 0, 0]], dtype=int)

    edges = response_order_directed_edges(signs)

    assert {tuple(edge) for edge in edges.tolist()} == {(0, 1), (2, 0)}


def test_has_response_order_cycle() -> None:
    acyclic = np.asarray([[0, -1, -1], [1, 0, -1], [1, 1, 0]], dtype=int)
    cyclic = np.asarray([[0, -1, 1], [1, 0, -1], [-1, 1, 0]], dtype=int)

    assert not has_response_order_cycle(acyclic)
    assert has_response_order_cycle(cyclic)


def test_response_order_topological_ranks() -> None:
    signs = np.asarray([[0, -1, -1], [1, 0, -1], [1, 1, 0]], dtype=int)

    ranks = response_order_topological_ranks(signs)

    assert np.array_equal(ranks, np.asarray([0, 1, 2]))


def test_response_order_topological_ranks_rejects_cycle() -> None:
    cyclic = np.asarray([[0, -1, 1], [1, 0, -1], [-1, 1, 0]], dtype=int)

    with pytest.raises(ValueError):
        response_order_topological_ranks(cyclic)


def test_scalar_rank_representation_error() -> None:
    signs = np.asarray([[0, -1], [1, 0]], dtype=int)

    good = scalar_rank_representation_error(signs, np.asarray([0, 1]))
    bad = scalar_rank_representation_error(signs, np.asarray([1, 0]))

    assert good["violation_count"] == 0.0
    assert bad["violation_fraction"] == 1.0


def test_scalar_representability_report() -> None:
    signs = np.asarray([[0, -1, 0], [1, 0, -1], [0, 1, 0]], dtype=int)

    report = scalar_representability_report(signs)

    assert report["target_count"] == 3.0
    assert report["nonzero_pair_count"] == 2.0
    assert report["scalar_representable"] == 1.0
