from __future__ import annotations

from causal_spacetime_lab.state_change_handoff_provenance import (
    HandoffProvenanceMetadata,
    default_bottom_up_handoff_provenance,
    default_hybrid_handoff_provenance,
    default_report_only_handoff_provenance,
    default_top_down_handoff_provenance,
    validate_handoff_provenance,
)


def test_bottom_up_handoff_provenance_validates() -> None:
    provenance = default_bottom_up_handoff_provenance(
        design_source_label="unit",
        design_source_path="tests",
        design_digest="abc",
    )

    assert validate_handoff_provenance(provenance)["valid_provenance"] == 1.0


def test_top_down_handoff_provenance_requires_design_digest() -> None:
    valid = default_top_down_handoff_provenance(
        design_source_label="unit",
        design_source_path="tests/spec.json",
        design_digest="abc",
        template_id="template",
        template_hash="hash",
    )
    invalid = HandoffProvenanceMetadata(
        **{**valid.__dict__, "design_digest": "", "design_source_path": ""}
    )

    assert validate_handoff_provenance(valid)["valid_provenance"] == 1.0
    assert validate_handoff_provenance(invalid)["valid_provenance"] == 0.0


def test_hybrid_handoff_provenance_requires_template_id_and_hash() -> None:
    valid = default_hybrid_handoff_provenance(
        design_source_label="unit",
        design_source_path="tests/spec.json",
        design_digest="abc",
        constraint_selection_policy="profile_instantiated_rank_gap",
        template_id="template",
        template_hash="hash",
    )
    invalid = HandoffProvenanceMetadata(
        **{**valid.__dict__, "template_id": "", "template_hash": ""}
    )

    assert validate_handoff_provenance(valid)["valid_provenance"] == 1.0
    assert validate_handoff_provenance(invalid)["valid_provenance"] == 0.0


def test_provenance_using_evaluation_outputs_is_blocked() -> None:
    provenance = default_hybrid_handoff_provenance(
        design_source_label="unit",
        design_source_path="tests/spec.json",
        design_digest="abc",
        constraint_selection_policy="profile_instantiated_rank_gap",
        template_id="template",
        template_hash="hash",
    )
    fit_used = HandoffProvenanceMetadata(
        **{**provenance.__dict__, "fit_outputs_used": True}
    )
    carry_used = HandoffProvenanceMetadata(
        **{**provenance.__dict__, "carry_forward_outputs_used": True}
    )

    assert validate_handoff_provenance(fit_used)["valid_provenance"] == 0.0
    assert validate_handoff_provenance(carry_used)["valid_provenance"] == 0.0


def test_report_only_control_requires_report_only_true() -> None:
    provenance = default_report_only_handoff_provenance(
        design_source_label="unit",
        design_source_path="tests/spec.json",
        design_digest="abc",
        reason="control",
    )
    invalid = HandoffProvenanceMetadata(
        **{**provenance.__dict__, "report_only": False}
    )

    assert validate_handoff_provenance(provenance)["valid_provenance"] == 1.0
    assert validate_handoff_provenance(invalid)["valid_provenance"] == 0.0
