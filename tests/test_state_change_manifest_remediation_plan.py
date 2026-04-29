from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
    build_remediation_actions_from_failure_summary,
    default_new_manifest_family_specs_v2,
    remediation_plan_digest,
    remediation_plan_to_jsonable,
)


def test_build_remediation_actions_maps_heldout_failure() -> None:
    actions = build_remediation_actions_from_failure_summary(
        [{"root_cause_category": "heldout_failure"}],
        [],
    )

    assert any(action.action_name == "add_protocol_columns" for action in actions)


def test_build_remediation_actions_maps_missing_metric() -> None:
    actions = build_remediation_actions_from_failure_summary(
        [{"root_cause_category": "missing_metric"}],
        [{"missing_metrics": "target_coverage_fraction"}],
    )

    assert any(
        action.action_name == "add_missing_metric_collection"
        for action in actions
    )


def test_default_new_manifest_family_specs_v2_all_planned_only() -> None:
    specs = default_new_manifest_family_specs_v2()

    assert specs
    assert all(spec["execution_status"] == "planned_only" for spec in specs)


def test_build_preregistered_remediation_plan_execution_disallowed() -> None:
    plan = build_preregistered_remediation_plan([], [])

    assert not plan.execution_allowed_in_current_milestone
    assert plan.diagnostic_requirements


def test_remediation_plan_digest_stable() -> None:
    plan = build_preregistered_remediation_plan([], [])
    payload = remediation_plan_to_jsonable(plan)

    assert remediation_plan_digest(payload) == remediation_plan_digest(payload)
