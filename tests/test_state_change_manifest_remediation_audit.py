from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    build_future_manifest_run_spec_from_plan,
)
from causal_spacetime_lab.state_change_manifest_remediation_audit import (
    future_run_spec_audit,
    remediation_plan_execution_audit,
)
from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
)


def test_remediation_plan_execution_audit_passes_for_report_only_plan() -> None:
    plan = build_preregistered_remediation_plan([], [])
    row = remediation_plan_execution_audit(plan)

    assert row["passed"] == 1.0


def test_future_run_spec_audit_passes_for_non_executable_spec() -> None:
    plan = build_preregistered_remediation_plan([], [])
    spec = build_future_manifest_run_spec_from_plan(plan)
    row = future_run_spec_audit(spec)

    assert row["passed"] == 1.0
