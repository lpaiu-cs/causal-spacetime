"""Provenance metadata for response-comparison handoff manifests."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

PROVENANCE_TYPES = {
    "bottom_up_profile_derived",
    "top_down_preregistered_template",
    "hybrid_template_instantiated_from_profile",
    "report_only_control",
}
FROZEN_BEFORE_STAGES = {
    "profile_generation",
    "constraint_instantiation",
    "representation_fit",
    "carry_forward_evaluation",
    "stress_test",
}


@dataclass(frozen=True)
class HandoffProvenanceMetadata:
    """Provenance for a frozen response-comparison handoff manifest."""

    provenance_type: str
    design_source_label: str
    design_source_path: str
    design_digest: str
    created_by_milestone: str
    frozen_before_stage: str
    allowed_data_dependencies: list[str]
    forbidden_data_dependencies: list[str]
    constraint_selection_policy: str
    template_id: str
    template_hash: str
    profile_instantiation_required: bool
    fit_outputs_used: bool
    heldout_outputs_used: bool
    carry_forward_outputs_used: bool
    stress_test_outputs_used: bool
    report_only: bool
    reason_if_report_only: str

    def __post_init__(self) -> None:
        if self.provenance_type not in PROVENANCE_TYPES:
            allowed = ", ".join(sorted(PROVENANCE_TYPES))
            raise ValueError(f"provenance_type must be one of: {allowed}")
        if self.frozen_before_stage not in FROZEN_BEFORE_STAGES:
            allowed = ", ".join(sorted(FROZEN_BEFORE_STAGES))
            raise ValueError(f"frozen_before_stage must be one of: {allowed}")


def handoff_provenance_jsonable(
    metadata: HandoffProvenanceMetadata,
) -> dict[str, object]:
    """Return JSON-compatible handoff provenance."""

    return asdict(metadata)


def handoff_provenance_digest(metadata: HandoffProvenanceMetadata) -> str:
    """Return a digest over handoff provenance metadata."""

    payload = json.dumps(
        handoff_provenance_jsonable(metadata),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def default_bottom_up_handoff_provenance(
    *,
    design_source_label: str,
    design_source_path: str = "",
    design_digest: str = "",
) -> HandoffProvenanceMetadata:
    """Return bottom-up profile-derived handoff provenance."""

    return HandoffProvenanceMetadata(
        provenance_type="bottom_up_profile_derived",
        design_source_label=design_source_label,
        design_source_path=design_source_path,
        design_digest=design_digest,
        created_by_milestone="Milestone 41",
        frozen_before_stage="representation_fit",
        allowed_data_dependencies=["generated_response_profile"],
        forbidden_data_dependencies=_forbidden_dependencies(),
        constraint_selection_policy="profile_instantiated_constraints",
        template_id="",
        template_hash="",
        profile_instantiation_required=True,
        fit_outputs_used=False,
        heldout_outputs_used=False,
        carry_forward_outputs_used=False,
        stress_test_outputs_used=False,
        report_only=False,
        reason_if_report_only="",
    )


def default_top_down_handoff_provenance(
    *,
    design_source_label: str,
    design_source_path: str,
    design_digest: str,
    template_id: str,
    template_hash: str,
) -> HandoffProvenanceMetadata:
    """Return top-down preregistered-template handoff provenance."""

    return HandoffProvenanceMetadata(
        provenance_type="top_down_preregistered_template",
        design_source_label=design_source_label,
        design_source_path=design_source_path,
        design_digest=design_digest,
        created_by_milestone="Milestone 41",
        frozen_before_stage="constraint_instantiation",
        allowed_data_dependencies=["preregistered_template"],
        forbidden_data_dependencies=_forbidden_dependencies(),
        constraint_selection_policy="declared_constraint_schema_only",
        template_id=template_id,
        template_hash=template_hash,
        profile_instantiation_required=False,
        fit_outputs_used=False,
        heldout_outputs_used=False,
        carry_forward_outputs_used=False,
        stress_test_outputs_used=False,
        report_only=False,
        reason_if_report_only="",
    )


def default_hybrid_handoff_provenance(
    *,
    design_source_label: str,
    design_source_path: str,
    design_digest: str,
    constraint_selection_policy: str,
    template_id: str,
    template_hash: str,
) -> HandoffProvenanceMetadata:
    """Return hybrid template-instantiated handoff provenance."""

    return HandoffProvenanceMetadata(
        provenance_type="hybrid_template_instantiated_from_profile",
        design_source_label=design_source_label,
        design_source_path=design_source_path,
        design_digest=design_digest,
        created_by_milestone="Milestone 41",
        frozen_before_stage="constraint_instantiation",
        allowed_data_dependencies=[
            "preregistered_template",
            "generated_response_profile",
        ],
        forbidden_data_dependencies=_forbidden_dependencies(),
        constraint_selection_policy=constraint_selection_policy,
        template_id=template_id,
        template_hash=template_hash,
        profile_instantiation_required=True,
        fit_outputs_used=False,
        heldout_outputs_used=False,
        carry_forward_outputs_used=False,
        stress_test_outputs_used=False,
        report_only=False,
        reason_if_report_only="",
    )


def default_report_only_handoff_provenance(
    *,
    design_source_label: str,
    design_source_path: str = "",
    design_digest: str = "",
    reason: str = "report-only control",
) -> HandoffProvenanceMetadata:
    """Return report-only control handoff provenance."""

    return HandoffProvenanceMetadata(
        provenance_type="report_only_control",
        design_source_label=design_source_label,
        design_source_path=design_source_path,
        design_digest=design_digest,
        created_by_milestone="Milestone 41",
        frozen_before_stage="constraint_instantiation",
        allowed_data_dependencies=["accounting_control"],
        forbidden_data_dependencies=_forbidden_dependencies(),
        constraint_selection_policy="report_only_no_constraints",
        template_id="report_only_control",
        template_hash=design_digest,
        profile_instantiation_required=False,
        fit_outputs_used=False,
        heldout_outputs_used=False,
        carry_forward_outputs_used=False,
        stress_test_outputs_used=False,
        report_only=True,
        reason_if_report_only=reason,
    )


def validate_handoff_provenance(
    metadata: HandoffProvenanceMetadata,
) -> dict[str, float | str]:
    """Validate handoff provenance without interpreting it as evidence."""

    failed: list[str] = []
    if metadata.provenance_type == "top_down_preregistered_template":
        if not (metadata.design_source_path or metadata.design_digest):
            failed.append("missing_top_down_design_source")
    if metadata.provenance_type == "hybrid_template_instantiated_from_profile":
        if not metadata.template_id or not metadata.template_hash:
            failed.append("missing_hybrid_template")
        if not metadata.design_digest:
            failed.append("missing_hybrid_design_digest")
    if metadata.fit_outputs_used:
        failed.append("fit_outputs_used")
    if metadata.heldout_outputs_used:
        failed.append("heldout_outputs_used_for_selection")
    if metadata.carry_forward_outputs_used:
        failed.append("carry_forward_outputs_used")
    if metadata.stress_test_outputs_used:
        failed.append("stress_test_outputs_used")
    if metadata.provenance_type == "report_only_control" and not metadata.report_only:
        failed.append("report_only_control_not_report_only")
    if metadata.provenance_type != "report_only_control" and metadata.report_only:
        failed.append("production_manifest_marked_report_only")
    return {
        "provenance_type": metadata.provenance_type,
        "valid_provenance": float(not failed),
        "failed_reasons": ";".join(failed),
        "report_only": float(metadata.report_only),
    }


def _forbidden_dependencies() -> list[str]:
    return [
        "representation_fit_outputs",
        "heldout_fit_results",
        "carry_forward_decisions",
        "stress_test_results",
    ]
