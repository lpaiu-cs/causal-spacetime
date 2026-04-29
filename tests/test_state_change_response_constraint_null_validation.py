from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_null_validation import (
    evaluate_constraint_pool_against_nulls,
    null_separation_by_type,
)
from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _profile() -> EchoResponseProfile:
    return EchoResponseProfile(
        target_event_ids=np.arange(5, dtype=int),
        protocol_labels=["a", "b", "c"],
        delay_rank_matrix=np.asarray(
            [
                [1, 1, 2],
                [2, 2, 3],
                [3, 3, 4],
                [5, 5, 6],
                [8, 8, 9],
            ],
            dtype=int,
        ),
        reachable_matrix=np.ones((5, 3), dtype=bool),
    )


def test_evaluate_constraint_pool_against_nulls_returns_expected_keys() -> None:
    profile = _profile()
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    pool = build_constraint_pool_from_dissimilarity(
        pairwise_response_dissimilarity(profile, protocol),
        max_constraints=30,
        seed=0,
    )

    report = evaluate_constraint_pool_against_nulls(
        profile,
        protocol,
        pool,
        null_repetitions=2,
        seed=0,
    )

    assert "null_z_score" in report
    assert report["null_repetition_count"] == 2.0


def test_null_separation_by_type_returns_one_row_per_null_type() -> None:
    profile = _profile()
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    pool = build_constraint_pool_from_dissimilarity(
        pairwise_response_dissimilarity(profile, protocol),
        max_constraints=30,
        seed=0,
    )

    rows = null_separation_by_type(
        profile,
        protocol,
        pool,
        null_repetitions=2,
        seed=0,
    )

    assert {row["null_type"] for row in rows} == {
        "shuffle_delays",
        "shuffle_reachability",
        "permute_profiles",
        "random_same_marginals",
    }
