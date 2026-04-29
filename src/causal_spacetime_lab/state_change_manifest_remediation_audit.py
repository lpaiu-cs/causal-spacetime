"""Audits for remediation plans and future-run specifications."""

from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    FutureManifestRunSpec,
)
from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    RemediationPlan,
)


def remediation_plan_execution_audit(
    plan: RemediationPlan,
) -> dict[str, float | str]:
    """Audit that a remediation plan is report-only in the current milestone."""

    execution_changing_actions_ok = all(
        action.requires_new_preregistration
        for action in plan.actions
        if action.action_type != "no_retuning_guard"
    )
    planned_only_specs = all(
        spec.get("execution_status") == "planned_only"
        for spec in plan.new_manifest_family_specs
    )
    execution_not_allowed = not plan.execution_allowed_in_current_milestone
    forbidden_nonempty = bool(plan.forbidden_interpretations)
    passed = all(
        [
            execution_not_allowed,
            execution_changing_actions_ok,
            forbidden_nonempty,
            planned_only_specs,
        ]
    )
    return {
        "plan_id": plan.plan_id,
        "execution_allowed_in_current_milestone": float(
            plan.execution_allowed_in_current_milestone
        ),
        "execution_not_allowed": float(execution_not_allowed),
        "execution_changing_actions_require_preregistration": float(
            execution_changing_actions_ok
        ),
        "forbidden_interpretations_nonempty": float(forbidden_nonempty),
        "new_manifest_family_specs_planned_only": float(planned_only_specs),
        "passed": float(passed),
    }


def future_run_spec_audit(
    spec: FutureManifestRunSpec,
) -> dict[str, float | str]:
    """Audit that a future-run spec is non-executable now."""

    execution_not_allowed = not spec.allowed_to_execute_now
    required_metrics_nonempty = bool(spec.required_metrics)
    forbidden_text = ";".join(spec.forbidden_actions)
    forbidden_actions_present = (
        "threshold retuning" in forbidden_text
        and "stress tests" in forbidden_text
    )
    fixed_threshold_source_declared = bool(spec.fixed_threshold_source)
    passed = all(
        [
            execution_not_allowed,
            required_metrics_nonempty,
            forbidden_actions_present,
            fixed_threshold_source_declared,
        ]
    )
    return {
        "run_name": spec.run_name,
        "allowed_to_execute_now": float(spec.allowed_to_execute_now),
        "execution_not_allowed": float(execution_not_allowed),
        "required_metrics_nonempty": float(required_metrics_nonempty),
        "forbidden_actions_include_required_stops": float(
            forbidden_actions_present
        ),
        "fixed_threshold_source_declared": float(
            fixed_threshold_source_declared
        ),
        "passed": float(passed),
    }
