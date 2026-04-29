from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    HandoffDecision,
    HandoffValidationSummary,
    ResponseConstraintHandoffManifest,
    build_handoff_validation_summary,
    decide_handoff_eligibility,
    forbidden_interpretations_default,
    manifest_digest,
    manifest_to_jsonable,
    read_handoff_manifest,
    split_constraint_indices,
    write_handoff_manifest,
)


def _summary(**overrides: float | int) -> HandoffValidationSummary:
    values: dict[str, float | int] = {
        "heldout_agreement_fraction": 0.9,
        "heldout_inversion_fraction": 0.0,
        "heldout_evaluable_fraction": 1.0,
        "bootstrap_mean_agreement_fraction": 0.9,
        "bootstrap_stable_constraint_fraction": 0.9,
        "null_z_score": 2.0,
        "best_null_agreement_fraction": 0.5,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 1.0,
        "constraint_count": 200,
    }
    values.update(overrides)
    return HandoffValidationSummary(**values)  # type: ignore[arg-type]


def _manifest() -> ResponseConstraintHandoffManifest:
    return ResponseConstraintHandoffManifest(
        manifest_id="",
        created_by_milestone="31",
        profile_label="test",
        comparison_protocol_name="gap",
        comparison_method="rank_gap_mean",
        missing_policy="common_reachable",
        min_common_protocols=1,
        min_margin=0.05,
        max_constraints=10,
        constraint_seed=0,
        train_fraction=0.7,
        validation_gate_name="gate",
        validation_summary=_summary(),
        handoff_decision=HandoffDecision(True, [], []),
        target_event_ids=np.asarray([10, 11, 12], dtype=int),
        constraints=np.asarray([[0, 1, 1, 2]], dtype=int),
        margins=np.asarray([0.3], dtype=float),
        train_constraint_indices=np.asarray([0], dtype=int),
        heldout_constraint_indices=np.asarray([], dtype=int),
        null_baseline_labels=["shuffle_delays"],
        forbidden_interpretations=forbidden_interpretations_default(),
    )


def test_forbidden_interpretations_default_is_nonempty() -> None:
    assert forbidden_interpretations_default()


def test_split_constraint_indices_disjoint_and_exhaustive() -> None:
    train, heldout = split_constraint_indices(10, 0.6, seed=0)

    assert not set(train) & set(heldout)
    assert sorted(np.concatenate([train, heldout]).tolist()) == list(range(10))


def test_manifest_to_jsonable_handles_numpy_arrays() -> None:
    jsonable = manifest_to_jsonable(_manifest())

    assert isinstance(jsonable["target_event_ids"], list)
    assert isinstance(jsonable["constraints"], list)


def test_manifest_digest_stable_for_same_manifest() -> None:
    jsonable = manifest_to_jsonable(_manifest())

    assert manifest_digest(jsonable) == manifest_digest(jsonable)


def test_write_read_handoff_manifest(tmp_path) -> None:  # type: ignore[no-untyped-def]
    manifest = _manifest()
    path = write_handoff_manifest(manifest, tmp_path / "manifest.json")

    loaded = read_handoff_manifest(path)

    assert loaded["profile_label"] == "test"
    assert "forbidden_interpretations" in loaded


def test_decide_handoff_eligibility_pass_case() -> None:
    decision = decide_handoff_eligibility(_summary(), ConstraintValidationGate("g"))

    assert decision.eligible
    assert not decision.failed_reasons


def test_decide_handoff_eligibility_fail_cases() -> None:
    decision = decide_handoff_eligibility(
        _summary(
            constraint_count=1,
            heldout_agreement_fraction=0.1,
            target_coverage_fraction=0.1,
        ),
        ConstraintValidationGate("g"),
    )

    assert not decision.eligible
    assert "too_few_constraints" in decision.failed_reasons
    assert "low_heldout_agreement" in decision.failed_reasons
    assert "low_target_coverage" in decision.failed_reasons


def test_build_handoff_validation_summary_fields() -> None:
    summary = build_handoff_validation_summary(
        {
            "agreement_fraction": 0.8,
            "inversion_fraction": 0.1,
            "evaluable_fraction": 0.9,
        },
        {
            "mean_agreement_fraction": 0.85,
            "stable_constraint_fraction": 0.8,
        },
        {
            "null_z_score": 2.0,
            "best_null_agreement_fraction": 0.6,
        },
        {
            "touched_target_fraction": 1.0,
            "touched_pair_node_fraction": 0.9,
        },
        100,
    )

    assert summary.constraint_count == 100
    assert summary.null_z_score == 2.0
