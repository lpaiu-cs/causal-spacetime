"""Preregistered upstream remediation plan data structures."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    DiagnosticMetricRequirement,
    default_diagnostic_metric_requirements,
)

ACTION_TYPES = {
    "reporting_completeness",
    "metric_collection",
    "profile_generation_design",
    "manifest_generation_design",
    "null_taxonomy_design",
    "stability_diagnostic_design",
    "no_retuning_guard",
}


@dataclass(frozen=True)
class RemediationAction:
    """One planned upstream remediation action."""

    action_name: str
    target_root_cause: str
    action_type: str
    description: str
    expected_effect: str
    requires_new_preregistration: bool
    allowed_in_current_milestone: bool

    def __post_init__(self) -> None:
        if self.action_type not in ACTION_TYPES:
            allowed = ", ".join(sorted(ACTION_TYPES))
            raise ValueError(f"action_type must be one of: {allowed}")


@dataclass(frozen=True)
class RemediationPlan:
    """A report-only plan for future diagnostic-complete manifest generation."""

    plan_id: str
    created_by_milestone: str
    input_failure_report_files: list[str]
    actions: list[RemediationAction]
    diagnostic_requirements: list[DiagnosticMetricRequirement]
    new_manifest_family_specs: list[dict[str, float | str]]
    forbidden_interpretations: list[str]
    execution_allowed_in_current_milestone: bool


def forbidden_remediation_interpretations() -> list[str]:
    """Return forbidden interpretations for remediation planning."""

    return [
        "remediation recovers geometry",
        "missing metric counts as pass",
        "thresholds changed after failure",
        "blocked family can be stress-tested anyway",
        "new manifest design proves metric structure",
    ]


def _action(
    action_name: str,
    target_root_cause: str,
    action_type: str,
    description: str,
    expected_effect: str,
) -> RemediationAction:
    return RemediationAction(
        action_name=action_name,
        target_root_cause=target_root_cause,
        action_type=action_type,
        description=description,
        expected_effect=expected_effect,
        requires_new_preregistration=True,
        allowed_in_current_milestone=False,
    )


def _action_catalog() -> dict[str, RemediationAction]:
    actions = [
        _action(
            "add_protocol_columns",
            "heldout_failure",
            "profile_generation_design",
            "Predeclare response profiles with more protocol columns.",
            "May improve held-out agreement in a future manifest family.",
        ),
        _action(
            "improve_profile_diversity",
            "heldout_failure",
            "profile_generation_design",
            "Predeclare more diverse protocol labels and emissions.",
            "May reduce family dependence on a narrow protocol slice.",
        ),
        _action(
            "increase_rank_resolution",
            "generalization_gap_failure",
            "profile_generation_design",
            "Predeclare rank-resolution enrichment before manifest creation.",
            "May reduce tie-heavy profiles and weak held-out transfer.",
        ),
        _action(
            "reduce_overfit_constraints",
            "generalization_gap_failure",
            "manifest_generation_design",
            "Predeclare constraint sampling balance and margin reporting.",
            "May reduce train-specific constraint pool effects.",
        ),
        _action(
            "strengthen_null_taxonomy_reporting",
            "null_separation_failure",
            "null_taxonomy_design",
            "Report destructive, symmetry-control, and marginal-preserving nulls.",
            "Improves null-taxonomy completeness in future bundles.",
        ),
        _action(
            "improve_target_pair_coverage",
            "coverage_failure",
            "metric_collection",
            "Predeclare target and pair-node coverage production.",
            "Makes coverage weakness directly measurable.",
        ),
        _action(
            "add_restart_stability_outputs",
            "restart_instability",
            "stability_diagnostic_design",
            "Predeclare family-level restart-stability outputs.",
            "Makes optimizer sensitivity visible in future reports.",
        ),
        _action(
            "add_latent_order_stability_outputs",
            "latent_order_instability",
            "stability_diagnostic_design",
            "Predeclare latent-order stability outputs.",
            "Makes fitted latent-order variability visible in future reports.",
        ),
        _action(
            "add_missing_metric_collection",
            "missing_metric",
            "metric_collection",
            "Predeclare production of every required diagnostic metric.",
            "Separates unavailable diagnostics from measured failures.",
        ),
        _action(
            "add_output_bundle_completeness_check",
            "missing_output",
            "reporting_completeness",
            "Predeclare output-bundle completeness checks.",
            "Prevents absent CSV inputs from being silently ignored.",
        ),
        _action(
            "preserve_no_retuning_guard",
            "no_retuning",
            "no_retuning_guard",
            "Keep fixed criteria tied to their exported threshold source.",
            "Prevents post-result threshold changes in future execution.",
        ),
    ]
    return {action.action_name: action for action in actions}


def _add_action(
    selected: dict[str, RemediationAction],
    catalog: dict[str, RemediationAction],
    action_name: str,
) -> None:
    selected[action_name] = catalog[action_name]


def build_remediation_actions_from_failure_summary(
    failure_summary_rows: list[dict[str, str | float]],
    completeness_rows: list[dict[str, str | float]],
) -> list[RemediationAction]:
    """Map failure summaries and missing metrics to planned remediation actions."""

    catalog = _action_catalog()
    selected: dict[str, RemediationAction] = {}
    root_to_actions = {
        "heldout_failure": ["add_protocol_columns", "improve_profile_diversity"],
        "generalization_gap_failure": [
            "increase_rank_resolution",
            "reduce_overfit_constraints",
        ],
        "null_separation_failure": ["strengthen_null_taxonomy_reporting"],
        "coverage_failure": ["improve_target_pair_coverage"],
        "restart_instability": ["add_restart_stability_outputs"],
        "latent_order_instability": ["add_latent_order_stability_outputs"],
        "missing_metric": ["add_missing_metric_collection"],
        "missing_output": ["add_output_bundle_completeness_check"],
    }
    for row in failure_summary_rows:
        root = str(row.get("root_cause_category", ""))
        for action_name in root_to_actions.get(root, []):
            _add_action(selected, catalog, action_name)
    if not failure_summary_rows and not completeness_rows:
        _add_action(selected, catalog, "add_output_bundle_completeness_check")
    for row in completeness_rows:
        missing_metrics = str(row.get("missing_metrics", ""))
        if missing_metrics:
            _add_action(selected, catalog, "add_missing_metric_collection")
        if "target_coverage_fraction" in missing_metrics:
            _add_action(selected, catalog, "improve_target_pair_coverage")
        if "pair_node_coverage_fraction" in missing_metrics:
            _add_action(selected, catalog, "improve_target_pair_coverage")
        if "restart_std" in missing_metrics:
            _add_action(selected, catalog, "add_restart_stability_outputs")
        if "latent_order_disagreement" in missing_metrics:
            _add_action(selected, catalog, "add_latent_order_stability_outputs")
    _add_action(selected, catalog, "preserve_no_retuning_guard")
    return [selected[name] for name in sorted(selected)]


def default_new_manifest_family_specs_v2() -> list[dict[str, float | str]]:
    """Return planned v2 manifest-family specifications."""

    common = {
        "profile_column_policy": "diagnostic_complete_predeclared_columns",
        "target_inclusion_policy": "predeclared_reachability_and_identity_filter",
        "pair_node_coverage_policy": "report_target_and_pair_node_coverage",
        "null_taxonomy_policy": "destructive_symmetry_and_marginal_nulls",
        "restart_stability_required": 1.0,
        "latent_order_stability_required": 1.0,
        "execution_status": "planned_only",
    }
    return [
        {
            "family_name": "rank_gap_more_protocol_columns_v2",
            "family_kind": "structured",
            **common,
        },
        {
            "family_name": "rank_gap_coverage_enriched_v2",
            "family_kind": "structured",
            **common,
        },
        {
            "family_name": "rank_gap_rank_resolution_enriched_v2",
            "family_kind": "structured",
            **common,
        },
        {
            "family_name": "combined_diagnostic_complete_v2",
            "family_kind": "structured",
            **common,
        },
        {
            "family_name": "failed_controls_v2",
            "family_kind": "failed_control",
            **common,
        },
    ]


def build_preregistered_remediation_plan(
    failure_summary_rows: list[dict[str, str | float]],
    completeness_rows: list[dict[str, str | float]],
) -> RemediationPlan:
    """Build the Milestone 36 report-only remediation plan."""

    return RemediationPlan(
        plan_id="remediation_plan_m36",
        created_by_milestone="Milestone 36",
        input_failure_report_files=[
            "outputs/data/carry_forward_failure_summary.csv",
            "outputs/data/cross_family_diagnostic_completeness_audit.csv",
        ],
        actions=build_remediation_actions_from_failure_summary(
            failure_summary_rows,
            completeness_rows,
        ),
        diagnostic_requirements=default_diagnostic_metric_requirements(),
        new_manifest_family_specs=default_new_manifest_family_specs_v2(),
        forbidden_interpretations=forbidden_remediation_interpretations(),
        execution_allowed_in_current_milestone=False,
    )


def remediation_plan_to_jsonable(plan: RemediationPlan) -> dict[str, object]:
    """Convert a remediation plan to canonical JSON-compatible data."""

    return asdict(plan)


def remediation_plan_digest(plan_jsonable: dict[str, object]) -> str:
    """Return a SHA256 digest for canonical remediation-plan JSON."""

    payload = json.dumps(plan_jsonable, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_remediation_plan(
    plan: RemediationPlan,
    output_path: Path,
) -> Path:
    """Write a remediation plan JSON artifact."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = remediation_plan_to_jsonable(plan)
    payload["plan_digest"] = remediation_plan_digest(payload)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return output_path
