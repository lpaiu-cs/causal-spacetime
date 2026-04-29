"""Future manifest-run specification derived from a remediation plan."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    RemediationPlan,
)

RUN_STATUSES = {
    "planned_only",
    "requires_new_preregistration",
    "ready_for_future_execution",
}


@dataclass(frozen=True)
class FutureManifestRunSpec:
    """A non-executable future-run specification."""

    run_name: str
    run_status: str
    allowed_to_execute_now: bool
    prerequisite_milestone: str
    planned_manifest_families: list[str]
    required_output_files: list[str]
    required_metrics: list[str]
    fixed_threshold_source: str
    forbidden_actions: list[str]

    def __post_init__(self) -> None:
        if self.run_status not in RUN_STATUSES:
            allowed = ", ".join(sorted(RUN_STATUSES))
            raise ValueError(f"run_status must be one of: {allowed}")


def build_future_manifest_run_spec_from_plan(
    plan: RemediationPlan,
) -> FutureManifestRunSpec:
    """Build a future-run spec from a preregistered remediation plan."""

    required_files = sorted(
        {
            requirement.source_output_file
            for requirement in plan.diagnostic_requirements
            if requirement.source_output_file
        }
    )
    return FutureManifestRunSpec(
        run_name=f"{plan.plan_id}_future_manifest_run",
        run_status="requires_new_preregistration",
        allowed_to_execute_now=False,
        prerequisite_milestone="Milestone 37 or later",
        planned_manifest_families=[
            str(spec["family_name"]) for spec in plan.new_manifest_family_specs
        ],
        required_output_files=required_files,
        required_metrics=[
            requirement.metric_name
            for requirement in plan.diagnostic_requirements
        ],
        fixed_threshold_source="outputs/data/cross_family_robustness_criteria_table.csv",
        forbidden_actions=[
            "threshold retuning",
            "stress tests",
            "new representation fits",
            "production handoff manifest export",
        ],
    )


def future_run_spec_to_jsonable(spec: FutureManifestRunSpec) -> dict[str, object]:
    """Convert a future-run spec to JSON-compatible data."""

    return asdict(spec)
