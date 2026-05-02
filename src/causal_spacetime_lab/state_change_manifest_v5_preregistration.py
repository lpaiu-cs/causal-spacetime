"""Planned-only v5 protocol preregistration utilities."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_v5_design import (
    V5ProtocolFamilyDesign,
    default_v5_protocol_family_designs,
    v5_protocol_family_design_table,
)


@dataclass(frozen=True)
class V5ProtocolPreregistrationSpec:
    """Planned-only v5 protocol preregistration specification."""

    spec_id: str
    created_by_milestone: str
    based_on_outputs: list[str]
    fixed_criteria_source: str
    planned_families: list[V5ProtocolFamilyDesign]
    required_metrics: list[str]
    required_output_files: list[str]
    forbidden_actions: list[str]
    forbidden_interpretations: list[str]
    execution_allowed_in_current_milestone: bool
    remediation_iteration_risk_audit_required: bool


def forbidden_v5_interpretations() -> list[str]:
    """Return forbidden interpretations for planned-only v5 design."""

    return [
        "v5 design proves geometry",
        "v5 family will pass",
        "thresholds changed after v4 failure",
        "blocked v4 families can be stress-tested",
        "planned v5 manifests are current results",
        "top-down provenance proves theory",
        "more remediation guarantees success",
    ]


def _forbidden_actions() -> list[str]:
    return [
        "threshold retuning",
        "v5 production manifest generation",
        "representation fitting",
        "carry-forward override",
        "stress tests",
        "rerun v4 fits",
    ]


def _required_output_files() -> list[str]:
    return [
        "outputs/manifests_v5/*.json",
        "outputs/data/v5_protocol_manifest_generation.csv",
        "outputs/data/v5_protocol_manifest_metadata_audit.csv",
        "outputs/data/v5_protocol_manifest_family_fit_summary.csv",
        "outputs/data/v5_protocol_manifest_family_null_taxonomy.csv",
        "outputs/data/v5_protocol_manifest_family_stricter_criteria.csv",
        "outputs/data/v5_protocol_manifest_family_failed_accounting.csv",
        "outputs/data/v5_protocol_manifest_family_coverage_metrics.csv",
        "outputs/data/v5_protocol_manifest_family_restart_stability.csv",
        "outputs/data/v5_protocol_manifest_family_latent_order_stability.csv",
        "outputs/data/v5_protocol_no_retuning_audit.csv",
        "outputs/data/v5_protocol_cross_family_robustness_metrics.csv",
        "outputs/data/v5_protocol_diagnostic_completeness_check.csv",
    ]


def build_v5_protocol_preregistration_spec(
    designs: list[V5ProtocolFamilyDesign],
    *,
    based_on_outputs: list[str] | None = None,
) -> V5ProtocolPreregistrationSpec:
    """Build planned-only v5 preregistration from M46 audit outputs."""

    if not designs:
        designs = default_v5_protocol_family_designs()
    required_metrics = [
        item.metric_name for item in default_diagnostic_metric_requirements()
    ]
    spec = V5ProtocolPreregistrationSpec(
        spec_id="pending",
        created_by_milestone="Milestone 46",
        based_on_outputs=based_on_outputs
        or [
            "outputs/data/v4_protocol_cross_family_robustness_decisions.csv",
            "outputs/data/v4_blocked_root_cause_audit.csv",
            "outputs/data/v4_criterion_margin_summary.csv",
            "outputs/data/v5_remediation_iteration_risk_audit.csv",
        ],
        fixed_criteria_source="Milestone 34 default_cross_family_robustness",
        planned_families=designs,
        required_metrics=required_metrics,
        required_output_files=_required_output_files(),
        forbidden_actions=_forbidden_actions(),
        forbidden_interpretations=forbidden_v5_interpretations(),
        execution_allowed_in_current_milestone=False,
        remediation_iteration_risk_audit_required=True,
    )
    jsonable = v5_protocol_preregistration_to_jsonable(spec)
    digest = v5_protocol_preregistration_digest({**jsonable, "spec_id": "pending"})
    return V5ProtocolPreregistrationSpec(
        spec_id=digest,
        created_by_milestone=spec.created_by_milestone,
        based_on_outputs=spec.based_on_outputs,
        fixed_criteria_source=spec.fixed_criteria_source,
        planned_families=spec.planned_families,
        required_metrics=spec.required_metrics,
        required_output_files=spec.required_output_files,
        forbidden_actions=spec.forbidden_actions,
        forbidden_interpretations=spec.forbidden_interpretations,
        execution_allowed_in_current_milestone=False,
        remediation_iteration_risk_audit_required=True,
    )


def v5_protocol_preregistration_to_jsonable(
    spec: V5ProtocolPreregistrationSpec,
) -> dict[str, object]:
    """Convert a v5 preregistration spec to canonical JSON-compatible data."""

    return {
        "spec_id": spec.spec_id,
        "created_by_milestone": spec.created_by_milestone,
        "based_on_outputs": list(spec.based_on_outputs),
        "fixed_criteria_source": spec.fixed_criteria_source,
        "planned_families": v5_protocol_family_design_table(spec.planned_families),
        "required_metrics": list(spec.required_metrics),
        "required_output_files": list(spec.required_output_files),
        "forbidden_actions": list(spec.forbidden_actions),
        "forbidden_interpretations": list(spec.forbidden_interpretations),
        "execution_allowed_in_current_milestone": bool(
            spec.execution_allowed_in_current_milestone
        ),
        "remediation_iteration_risk_audit_required": bool(
            spec.remediation_iteration_risk_audit_required
        ),
    }


def v5_protocol_preregistration_digest(jsonable: dict[str, object]) -> str:
    """Return SHA256 digest over canonical v5 preregistration JSON."""

    payload = json.dumps(jsonable, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_v5_protocol_preregistration_spec(
    spec: V5ProtocolPreregistrationSpec,
    output_path: Path,
) -> Path:
    """Write planned-only v5 preregistration JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(
            v5_protocol_preregistration_to_jsonable(spec),
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return output_path
