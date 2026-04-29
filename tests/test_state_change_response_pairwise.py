from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
    pairwise_response_dissimilarity_matrix,
    pairwise_response_order_inversion_rate,
    response_pair_comparison_constraints,
    unordered_target_pairs,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _profile() -> EchoResponseProfile:
    delays = np.asarray(
        [
            [1, 2, 3],
            [1, 2, 4],
            [2, 2, 3],
            [1, -1, 3],
        ],
        dtype=int,
    )
    return EchoResponseProfile(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        protocol_labels=["a", "b", "c"],
        delay_rank_matrix=delays,
        reachable_matrix=delays >= 0,
    )


def test_unordered_target_pairs_deterministic_case() -> None:
    pairs = unordered_target_pairs(np.asarray([10, 20, 30]))

    assert np.array_equal(pairs, np.asarray([[0, 1], [0, 2], [1, 2]]))


def test_pairwise_response_dissimilarity_matrix_separation_fraction() -> None:
    protocol = PairwiseResponseComparisonProtocol("sep", "separation_fraction")
    matrix = pairwise_response_dissimilarity_matrix(_profile(), protocol)

    assert np.isclose(matrix[0, 1], 1.0 / 3.0)
    assert np.allclose(matrix, matrix.T, equal_nan=True)


def test_pairwise_response_dissimilarity_matrix_rank_gap_mean() -> None:
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    matrix = pairwise_response_dissimilarity_matrix(_profile(), protocol)

    assert matrix[0, 1] > 0.0
    assert matrix[0, 0] == 0.0


def test_missing_policies() -> None:
    profile = _profile()
    common = PairwiseResponseComparisonProtocol("common", "rank_gap_mean")
    require_all = PairwiseResponseComparisonProtocol(
        "all",
        "rank_gap_mean",
        missing_policy="require_all_reachable",
    )
    penalize = PairwiseResponseComparisonProtocol(
        "penalty",
        "combined_gap_and_mismatch",
        missing_policy="penalize_mismatch",
    )

    common_matrix = pairwise_response_dissimilarity_matrix(profile, common)
    all_matrix = pairwise_response_dissimilarity_matrix(profile, require_all)
    penalty_matrix = pairwise_response_dissimilarity_matrix(profile, penalize)

    assert np.isfinite(common_matrix[0, 3])
    assert np.isnan(all_matrix[0, 3])
    assert np.isfinite(penalty_matrix[0, 3])


def test_protocol_validation() -> None:
    with pytest.raises(ValueError):
        PairwiseResponseComparisonProtocol("bad", "unknown")
    with pytest.raises(ValueError):
        PairwiseResponseComparisonProtocol(
            "bad",
            "rank_gap_mean",
            min_common_protocols=0,
        )


def test_response_pair_comparison_constraints_deterministic_profile() -> None:
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    dissimilarity = pairwise_response_dissimilarity(_profile(), protocol)

    constraints = response_pair_comparison_constraints(dissimilarity, 5, seed=0)

    assert constraints.shape[1] == 4
    assert constraints.shape[0] > 0


def test_pairwise_response_order_inversion_rate_identical_zero() -> None:
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    dissimilarity = pairwise_response_dissimilarity(_profile(), protocol)

    assert pairwise_response_order_inversion_rate(dissimilarity, dissimilarity) == 0.0
