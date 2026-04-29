from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    build_future_manifest_run_spec_from_plan,
)
from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
)


def test_build_future_manifest_run_spec_from_plan_disallows_execution_now() -> None:
    plan = build_preregistered_remediation_plan([], [])
    spec = build_future_manifest_run_spec_from_plan(plan)

    assert not spec.allowed_to_execute_now
    assert spec.run_status == "requires_new_preregistration"
    assert spec.required_metrics
