"""Planned-only v4 protocol preregistration utilities."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_v4_design import (
    V4ProtocolFamilyDesign,
    default_v4_protocol_family_designs,
    v4_protocol_family_design_table,
)


@dataclass(frozen=True)
class V4ProtocolPreregistrationSpec:
    """Planned-only v4 protocol preregistration specification."""

    spec_id: str
    created_by_milestone: str
    based_on_outputs: list[str]
    fixed_criteria_source: str
    planned_families: list[V4ProtocolFamilyDesign]
    required_metrics: list[str]
    required_output_files: list[str]
    forbidden_actions: list[str]
    forbidden_interpretations: list[str]
    execution_allowed_in_current_milestone: bool


def forbidden_v4_interpretations() -> list[str]:
    """Return forbidden interpretations for planned-only v4 design."""

    return [
        "v4 design proves geometry",
        "v4 family will pass",
        "thresholds changed after v3 failure",
        "blocked v3 families can be stress-tested",
        "planned v4 manifests are current results",
        "top-down provenance proves theory",
    ]


def _forbidden_actions() -> list[str]:
    return [
        "threshold retuning",
        "v4 production manifest generation",
        "representation fitting",
        "carry-forward override",
        "stress tests",
    ]


def _required_output_files() -> list[str]:
    return [
        "outputs/manifests_v4/*.json",
        "outputs/data/v4_protocol_manifest_generation.csv",
        "outputs/data/v4_protocol_manifest_metadata_audit.csv",
        "outputs/data/v4_protocol_manifest_family_fit_summary.csv",
        "outputs/data/v4_protocol_manifest_family_null_taxonomy.csv",
        "outputs/data/v4_protocol_manifest_family_stricter_criteria.csv",
        "outputs/data/v4_protocol_manifest_family_failed_accounting.csv",
        "outputs/data/v4_protocol_manifest_family_coverage_metrics.csv",
        "outputs/data/v4_protocol_manifest_family_restart_stability.csv",
        "outputs/data/v4_protocol_manifest_family_latent_order_stability.csv",
        "outputs/data/v4_protocol_no_retuning_audit.csv",
        "outputs/data/v4_protocol_cross_family_robustness_metrics.csv",
        "outputs/data/v4_protocol_diagnostic_completeness_check.csv",
    ]


def build_v4_protocol_preregistration_spec(
    designs: list[V4ProtocolFamilyDesign],
    *,
    based_on_outputs: list[str] | None = None,
) -> V4ProtocolPreregistrationSpec:
    """Build planned-only v4 preregistration from M43 audit outputs."""

    if not designs:
        designs = default_v4_protocol_family_designs()
    required_metrics = [
        item.metric_name for item in default_diagnostic_metric_requirements()
    ]
    spec = V4ProtocolPreregistrationSpec(
        spec_id="pending",
        created_by_milestone="Milestone 43",
        based_on_outputs=based_on_outputs
        or [
            "outputs/data/v3_protocol_cross_family_robustness_decisions.csv",
            "outputs/data/v3_protocol_blocked_root_cause_audit.csv",
            "outputs/data/v3_protocol_criterion_margin_summary.csv",
        ],
        fixed_criteria_source="Milestone 34 default_cross_family_robustness",
        planned_families=designs,
        required_metrics=required_metrics,
        required_output_files=_required_output_files(),
        forbidden_actions=_forbidden_actions(),
        forbidden_interpretations=forbidden_v4_interpretations(),
        execution_allowed_in_current_milestone=False,
    )
    jsonable = v4_protocol_preregistration_to_jsonable(spec)
    digest = v4_protocol_preregistration_digest({**jsonable, "spec_id": "pending"})
    return V4ProtocolPreregistrationSpec(
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
    )


def v4_protocol_preregistration_to_jsonable(
    spec: V4ProtocolPreregistrationSpec,
) -> dict[str, object]:
    """Convert a v4 preregistration spec to canonical JSON-compatible data."""

    return {
        "spec_id": spec.spec_id,
        "created_by_milestone": spec.created_by_milestone,
        "based_on_outputs": list(spec.based_on_outputs),
        "fixed_criteria_source": spec.fixed_criteria_source,
        "planned_families": v4_protocol_family_design_table(spec.planned_families),
        "required_metrics": list(spec.required_metrics),
        "required_output_files": list(spec.required_output_files),
        "forbidden_actions": list(spec.forbidden_actions),
        "forbidden_interpretations": list(spec.forbidden_interpretations),
        "execution_allowed_in_current_milestone": bool(
            spec.execution_allowed_in_current_milestone
        ),
    }


def v4_protocol_preregistration_digest(jsonable: dict[str, object]) -> str:
    """Return SHA256 digest over canonical v4 preregistration JSON."""

    payload = json.dumps(jsonable, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_v4_protocol_preregistration_spec(
    spec: V4ProtocolPreregistrationSpec,
    output_path: Path,
) -> Path:
    """Write planned-only v4 preregistration JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(
            v4_protocol_preregistration_to_jsonable(spec),
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return output_path

