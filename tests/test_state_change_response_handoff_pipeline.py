from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    HandoffDecision,
    HandoffValidationSummary,
    ResponseConstraintHandoffManifest,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
    select_eligible_manifests,
    summarize_handoff_manifests,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _profile() -> EchoResponseProfile:
    return EchoResponseProfile(
        target_event_ids=np.arange(6, dtype=int),
        protocol_labels=["a", "b", "c", "d"],
        delay_rank_matrix=np.asarray(
            [
                [1, 1, 2, 2],
                [2, 2, 3, 3],
                [3, 3, 4, 4],
                [5, 5, 6, 6],
                [8, 8, 9, 9],
                [13, 13, 14, 14],
            ],
            dtype=int,
        ),
        reachable_matrix=np.ones((6, 4), dtype=bool),
    )


def _manifest(eligible: bool) -> ResponseConstraintHandoffManifest:
    return ResponseConstraintHandoffManifest(
        manifest_id="id",
        created_by_milestone="31",
        profile_label="x",
        comparison_protocol_name="gap",
        comparison_method="rank_gap_mean",
        missing_policy="common_reachable",
        min_common_protocols=1,
        min_margin=0.0,
        max_constraints=10,
        constraint_seed=0,
        train_fraction=0.7,
        validation_gate_name="g",
        validation_summary=HandoffValidationSummary(
            0.9,
            0.0,
            1.0,
            0.9,
            0.9,
            2.0,
            0.5,
            1.0,
            1.0,
            10,
        ),
        handoff_decision=HandoffDecision(eligible, [] if eligible else ["x"], []),
        target_event_ids=np.arange(3, dtype=int),
        constraints=np.asarray([[0, 1, 1, 2]], dtype=int),
        margins=np.asarray([0.2], dtype=float),
        train_constraint_indices=np.asarray([0], dtype=int),
        heldout_constraint_indices=np.asarray([], dtype=int),
        null_baseline_labels=["shuffle_delays"],
        forbidden_interpretations=["metric distance"],
    )


def test_build_candidate_handoff_manifest_returns_manifest() -> None:
    manifest = build_candidate_handoff_manifest(
        _profile(),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
        ConstraintValidationGate("gate", min_constraint_count=1),
        max_constraints=20,
        bootstrap_count=3,
        null_repetitions=2,
        constraint_seed=0,
    )

    assert manifest.manifest_id
    assert manifest.constraints.shape[1] == 4
    assert manifest.created_by_milestone == "31"


def test_select_eligible_manifests_filters_correctly() -> None:
    selected = select_eligible_manifests([_manifest(True), _manifest(False)])

    assert len(selected) == 1
    assert selected[0].handoff_decision.eligible


def test_summarize_handoff_manifests_returns_rows() -> None:
    rows = summarize_handoff_manifests([_manifest(True)])

    assert rows
    assert rows[0]["eligible"] == 1.0
    assert "heldout_agreement" in rows[0]
