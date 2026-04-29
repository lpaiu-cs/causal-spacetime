from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_response_signature import (
    echo_response_signature_from_motif_rows,
    response_order_sign_matrix,
    signature_reachable_fraction,
    signature_strict_pair_fraction,
    signature_tie_fraction,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    consensus_response_order_matrix,
    response_order_cycle_count,
    stable_response_order_core,
)


def test_response_order_sign_matrix_deterministic_values() -> None:
    signs = response_order_sign_matrix([2, 5, 5, -1], [True, True, True, False])

    assert signs[0, 1] == -1
    assert signs[1, 0] == 1
    assert signs[1, 2] == 0
    assert signs[0, 3] == 0


def test_echo_response_signature_from_motif_rows() -> None:
    signature = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 10.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 11.0, "recovered_delay_rank": float("nan")},
        ]
    )

    assert np.array_equal(signature.target_event_ids, np.asarray([10, 11]))
    assert np.array_equal(signature.delay_ranks, np.asarray([2, -1]))
    assert np.array_equal(signature.reachable_mask, np.asarray([True, False]))


def test_signature_tie_and_strict_pair_fractions() -> None:
    signature = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 5.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 5.0},
        ]
    )

    assert signature_reachable_fraction(signature) == 1.0
    assert signature_tie_fraction(signature) == pytest.approx(1.0 / 3.0)
    assert signature_strict_pair_fraction(signature) == pytest.approx(2.0 / 3.0)


def test_compare_response_signatures_identical_signatures() -> None:
    signature = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 5.0},
        ]
    )

    comparison = compare_response_signatures(signature, signature)

    assert comparison["common_target_count"] == 2.0
    assert comparison["pair_agreement_fraction"] == 1.0


def test_compare_response_signatures_tie_changes() -> None:
    baseline = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 5.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 5.0},
        ]
    )
    variant = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 4.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 6.0},
        ]
    )

    comparison = compare_response_signatures(baseline, variant)

    assert comparison["pair_tie_changed_fraction"] > 0.0


def test_stable_response_order_core_deterministic_case() -> None:
    first = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 4.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 6.0},
        ]
    )
    second = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 5.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 6.0},
        ]
    )

    core = stable_response_order_core([first, second])

    assert float(core["stable_pair_fraction"]) == 1.0
    assert np.asarray(core["stable_order_sign_matrix"])[0, 2] == -1


def test_consensus_response_order_matrix_deterministic_case() -> None:
    first = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 4.0},
        ]
    )
    second = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 6.0},
        ]
    )

    consensus = consensus_response_order_matrix([first, second])

    assert consensus[0, 1] == -1
    assert consensus[1, 0] == 1


def test_response_order_cycle_count_cyclic_and_acyclic() -> None:
    cyclic = np.asarray([[0, -1, 1], [1, 0, -1], [-1, 1, 0]], dtype=int)
    acyclic = np.asarray([[0, -1, -1], [1, 0, -1], [1, 1, 0]], dtype=int)

    assert response_order_cycle_count(cyclic) == 1
    assert response_order_cycle_count(acyclic) == 0

