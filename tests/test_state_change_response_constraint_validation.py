from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
    bootstrap_constraint_stability,
    bootstrap_profile_protocol_columns,
    heldout_protocol_constraint_validation,
    split_profile_protocol_columns,
    validation_gate_pass_fail,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _profile() -> EchoResponseProfile:
    return EchoResponseProfile(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        protocol_labels=["a", "b", "c", "d"],
        delay_rank_matrix=np.asarray(
            [
                [1, 1, 2, 2],
                [2, 2, 3, 3],
                [4, 4, 5, 5],
                [7, 7, 8, 8],
            ],
            dtype=int,
        ),
        reachable_matrix=np.ones((4, 4), dtype=bool),
    )


def test_split_profile_protocol_columns() -> None:
    train, test = split_profile_protocol_columns(_profile(), 0.5, seed=0)

    assert np.array_equal(train.target_event_ids, test.target_event_ids)
    assert len(train.protocol_labels) == 2
    assert len(test.protocol_labels) == 2


def test_heldout_protocol_constraint_validation_returns_expected_fields() -> None:
    report = heldout_protocol_constraint_validation(
        _profile(),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
        max_constraints=20,
        seed=0,
    )

    assert report["constraint_count"] > 0
    assert "train_protocol_count" in report
    assert "test_valid_pair_fraction" in report


def test_validation_gate_pass_fail_pass_and_fail_cases() -> None:
    gate = ConstraintValidationGate(
        "test",
        min_constraint_count=2,
        min_agreement_fraction=0.5,
    )
    passing = validation_gate_pass_fail(
        {
            "constraint_count": 3.0,
            "evaluable_fraction": 1.0,
            "agreement_fraction": 0.9,
            "inversion_fraction": 0.0,
            "tie_or_unresolved_fraction": 0.1,
        },
        gate,
    )
    failing = validation_gate_pass_fail(
        {
            "constraint_count": 1.0,
            "evaluable_fraction": 0.1,
            "agreement_fraction": 0.2,
            "inversion_fraction": 0.5,
            "tie_or_unresolved_fraction": 0.8,
        },
        gate,
    )

    assert passing["passed"] == 1.0
    assert failing["passed"] == 0.0
    assert "constraint_count" in str(failing["failed_reasons"])


def test_bootstrap_profile_protocol_columns_shape() -> None:
    boot = bootstrap_profile_protocol_columns(_profile(), sample_fraction=0.5, seed=0)

    assert boot.delay_rank_matrix.shape == (4, 2)
    assert boot.reachable_matrix.shape == (4, 2)


def test_bootstrap_constraint_stability_returns_expected_keys() -> None:
    profile = _profile()
    protocol = PairwiseResponseComparisonProtocol("gap", "rank_gap_mean")
    pool = build_constraint_pool_from_dissimilarity(
        pairwise_response_dissimilarity(profile, protocol),
        max_constraints=20,
        seed=0,
    )

    report = bootstrap_constraint_stability(
        profile,
        protocol,
        pool,
        bootstrap_count=3,
        seed=0,
    )

    assert report["bootstrap_count"] == 3.0
    assert "stable_constraint_fraction" in report
