"""Preregistration JSON utilities for planned v3 manifest families."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_v3_design import (
    V3ManifestFamilyDesign,
    default_v3_manifest_family_designs,
    v3_manifest_family_design_table,
)


@dataclass(frozen=True)
class V3PreregistrationSpec:
    """Planned-only v3 preregistration specification."""

    spec_id: str
    created_by_milestone: str
    based_on_outputs: list[str]
    fixed_criteria_source: str
    planned_families: list[V3ManifestFamilyDesign]
    required_metrics: list[str]
    required_output_files: list[str]
    forbidden_actions: list[str]
    forbidden_interpretations: list[str]
    execution_allowed_in_current_milestone: bool


def forbidden_v3_interpretations() -> list[str]:
    """Return forbidden interpretations for v3 preregistration."""

    return [
        "v3 design proves geometry",
        "v3 family will pass",
        "thresholds changed after v2 failure",
        "blocked v2 families can be stress-tested",
        "planned v3 manifests are current results",
    ]


def _forbidden_actions() -> list[str]:
    return [
        "threshold retuning",
        "v3 production manifest generation",
        "representation fitting",
        "carry-forward override",
        "stress tests",
    ]


def _required_output_files() -> list[str]:
    return [
        "outputs/manifests_v3/*.json",
        "outputs/data/v3_manifest_generation.csv",
        "outputs/data/v3_manifest_family_fit_summary.csv",
        "outputs/data/v3_manifest_family_null_taxonomy.csv",
        "outputs/data/v3_manifest_family_stricter_criteria.csv",
        "outputs/data/v3_manifest_family_failed_accounting.csv",
        "outputs/data/v3_manifest_family_coverage_metrics.csv",
        "outputs/data/v3_manifest_family_restart_stability.csv",
        "outputs/data/v3_manifest_family_latent_order_stability.csv",
        "outputs/data/v3_no_retuning_audit.csv",
        "outputs/data/v3_cross_family_robustness_metrics.csv",
        "outputs/data/v3_diagnostic_completeness_check.csv",
    ]


def build_v3_preregistration_spec(
    designs: list[V3ManifestFamilyDesign],
    *,
    based_on_outputs: list[str] | None = None,
) -> V3PreregistrationSpec:
    """Build a planned-only v3 preregistration spec."""

    if not designs:
        designs = default_v3_manifest_family_designs()
    required_metrics = [
        item.metric_name for item in default_diagnostic_metric_requirements()
    ]
    spec = V3PreregistrationSpec(
        spec_id="pending",
        created_by_milestone="Milestone 39",
        based_on_outputs=based_on_outputs
        or [
            "outputs/data/v2_cross_family_robustness_decisions.csv",
            "outputs/data/v2_cross_family_robustness_decision_metrics.csv",
            "outputs/data/v2_blocked_root_cause_audit.csv",
        ],
        fixed_criteria_source="Milestone 34 default_cross_family_robustness",
        planned_families=designs,
        required_metrics=required_metrics,
        required_output_files=_required_output_files(),
        forbidden_actions=_forbidden_actions(),
        forbidden_interpretations=forbidden_v3_interpretations(),
        execution_allowed_in_current_milestone=False,
    )
    jsonable = v3_preregistration_to_jsonable(spec)
    digest = v3_preregistration_digest({**jsonable, "spec_id": "pending"})
    return V3PreregistrationSpec(
        spec_id=digest,
        created_by_milestone=spec.created_by_milestone,
        based_on_outputs=spec.based_on_outputs,
        fixed_criteria_source=spec.fixed_criteria_source,
        planned_families=spec.planned_families,
        required_metrics=spec.required_metrics,
        required_output_files=spec.required_output_files,
        forbidden_actions=spec.forbidden_actions,
        forbidden_interpretations=spec.forbidden_interpretations,
        execution_allowed_in_current_milestone=(
            spec.execution_allowed_in_current_milestone
        ),
    )


def v3_preregistration_to_jsonable(
    spec: V3PreregistrationSpec,
) -> dict[str, object]:
    """Convert a v3 preregistration spec to canonical JSON-compatible data."""

    return {
        "spec_id": spec.spec_id,
        "created_by_milestone": spec.created_by_milestone,
        "based_on_outputs": list(spec.based_on_outputs),
        "fixed_criteria_source": spec.fixed_criteria_source,
        "planned_families": v3_manifest_family_design_table(spec.planned_families),
        "required_metrics": list(spec.required_metrics),
        "required_output_files": list(spec.required_output_files),
        "forbidden_actions": list(spec.forbidden_actions),
        "forbidden_interpretations": list(spec.forbidden_interpretations),
        "execution_allowed_in_current_milestone": bool(
            spec.execution_allowed_in_current_milestone
        ),
    }


def v3_preregistration_digest(jsonable: dict[str, object]) -> str:
    """Return SHA256 digest over canonical JSON."""

    payload = json.dumps(jsonable, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_v3_preregistration_spec(
    spec: V3PreregistrationSpec,
    output_path: Path,
) -> Path:
    """Write planned-only v3 preregistration JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(v3_preregistration_to_jsonable(spec), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return output_path
