from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_echo import (
    compare_echo_orders,
    echo_delay_rank_for_emission,
    echo_order_matrix_from_delay_ranks,
    echo_reachability_from_emission,
    echo_return_position_for_emission,
)
from causal_spacetime_lab.state_change_echo_selection import (
    emission_positions_for_reference_chain,
    select_echo_reachable_targets,
    select_targets_after_emission,
)
from causal_spacetime_lab.state_change_order import transitive_closure_dag


def _hand_closure() -> np.ndarray:
    adjacency = np.zeros((5, 5), dtype=bool)
    adjacency[0, 2] = True
    adjacency[1, 2] = True
    adjacency[1, 3] = True
    adjacency[3, 4] = True
    adjacency[2, 4] = True
    return transitive_closure_dag(adjacency)


def test_echo_return_position_for_emission_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([1, 3, 4])

    assert (
        echo_return_position_for_emission(
            closure,
            reference,
            target_index=2,
            emission_position=0,
        )
        == 2
    )


def test_echo_delay_rank_for_emission_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([1, 3, 4])

    assert (
        echo_delay_rank_for_emission(
            closure,
            reference,
            target_index=2,
            emission_position=0,
        )
        == 2
    )


def test_echo_reachability_from_emission_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([1, 3, 4])

    returns, delays, reachable = echo_reachability_from_emission(
        closure,
        reference,
        emission_position=0,
    )

    assert np.array_equal(returns, np.asarray([-1, -1, 2, -1, -1]))
    assert np.array_equal(delays, np.asarray([-1, -1, 2, -1, -1]))
    assert np.array_equal(reachable, np.asarray([False, False, True, False, False]))


def test_echo_reachability_handles_no_return() -> None:
    closure = _hand_closure()
    reference = np.asarray([1, 3, 4])

    assert (
        echo_return_position_for_emission(
            closure,
            reference,
            target_index=0,
            emission_position=0,
        )
        is None
    )


def test_echo_order_matrix_from_delay_ranks() -> None:
    delays = np.asarray([-1, 2, 3, -1])
    reachable = np.asarray([False, True, True, False])

    matrix = echo_order_matrix_from_delay_ranks(delays, reachable)

    assert matrix[1, 2]
    assert not matrix[2, 1]
    assert not matrix[0, 1]


def test_compare_echo_orders_identical_delays_zero_inversion() -> None:
    delays = np.asarray([-1, 2, 3, 4])
    reachable = np.asarray([False, True, True, True])

    comparison = compare_echo_orders(delays, delays, reachable, reachable)

    assert comparison["common_reachable_count"] == 3.0
    assert comparison["order_inversion_rate"] == 0.0
    assert comparison["order_agreement_rate"] == 1.0


def test_emission_positions_for_reference_chain_strategies() -> None:
    reference = np.arange(6)

    assert np.array_equal(
        emission_positions_for_reference_chain(reference, "all_interior"),
        np.asarray([1, 2, 3, 4]),
    )
    quantiles = emission_positions_for_reference_chain(
        reference,
        "interior_quantiles",
        count=3,
    )
    early_middle_late = emission_positions_for_reference_chain(
        reference,
        "early_middle_late",
    )

    assert quantiles.size == 3
    assert quantiles[0] >= 1 and quantiles[-1] <= 4
    assert early_middle_late.size <= 3
    assert np.all((early_middle_late >= 1) & (early_middle_late <= 4))


def test_select_targets_after_emission_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([1, 3, 4])

    targets = select_targets_after_emission(
        closure,
        reference,
        emission_position=0,
        target_indices=np.asarray([0, 2]),
    )

    assert np.array_equal(targets, np.asarray([2]))


def test_select_echo_reachable_targets_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([1, 3, 4])

    targets = select_echo_reachable_targets(
        closure,
        reference,
        emission_position=0,
        target_indices=np.asarray([0, 2]),
    )

    assert np.array_equal(targets, np.asarray([2]))
