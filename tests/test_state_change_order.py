from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_order import (
    causal_interval_indices,
    immediate_trigger_adjacency,
    is_irreflexive,
    is_transitive,
    local_finiteness_report,
    topological_order_from_adjacency,
    transitive_closure_dag,
)


def test_immediate_trigger_adjacency_shape() -> None:
    network = generate_state_change_network(3, 15, seed=21)

    adjacency = immediate_trigger_adjacency(network)

    assert adjacency.shape == (15, 15)


def test_transitive_closure_dag_hand_coded_example() -> None:
    adjacency = np.zeros((4, 4), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 2] = True
    adjacency[0, 3] = True

    closure = transitive_closure_dag(adjacency)

    assert closure[0, 2]
    assert closure[0, 3]
    assert not closure[2, 0]


def test_is_irreflexive_and_transitive_for_closure() -> None:
    adjacency = np.zeros((3, 3), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 2] = True
    closure = transitive_closure_dag(adjacency)

    assert is_irreflexive(closure)
    assert is_transitive(closure)


def test_causal_interval_indices_hand_coded_example() -> None:
    adjacency = np.zeros((5, 5), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 3] = True
    adjacency[0, 2] = True
    adjacency[2, 3] = True
    adjacency[3, 4] = True
    closure = transitive_closure_dag(adjacency)

    interval = causal_interval_indices(closure, 0, 4)

    assert np.array_equal(interval, np.asarray([1, 2, 3]))


def test_topological_order_from_adjacency_rejects_cycle() -> None:
    adjacency = np.zeros((2, 2), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 0] = True

    with pytest.raises(ValueError):
        topological_order_from_adjacency(adjacency)


def test_local_finiteness_report_returns_finite_intervals() -> None:
    network = generate_state_change_network(4, 30, seed=22)
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    report = local_finiteness_report(closure)

    assert report["finite_interval_fraction"] == pytest.approx(1.0)
    assert report["comparable_pair_count"] > 0.0
    assert report["max_interval_size"] < 30.0
