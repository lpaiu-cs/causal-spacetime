from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    CarryForwardFamilyRecord,
    CarryForwardRegistry,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)
from causal_spacetime_lab.state_change_manifest_v2_stress_handoff import (
    build_v2_stress_test_plan,
    evaluate_v2_stress_test_stop_condition,
)


def test_evaluate_v2_stress_test_stop_condition_stop_case() -> None:
    decisions = [
        FamilyRobustnessDecision(
            family_name="family",
            family_kind="structured",
            decision="blocked",
            passed=False,
            failed_reasons=[],
            warning_reasons=[],
            key_metrics={},
        )
    ]

    report = evaluate_v2_stress_test_stop_condition(decisions)

    assert not report.stress_tests_allowed


def test_build_v2_stress_test_plan_blocks_blocked_families() -> None:
    registry = CarryForwardRegistry(
        registry_id="r",
        created_by_milestone="Milestone 38",
        criteria_name="default",
        records=[
            CarryForwardFamilyRecord(
                family_name="family",
                family_kind="structured",
                decision="blocked",
                manifest_ids=[],
                eligible_manifest_count=0,
                ineligible_manifest_count=0,
                failed_reasons=[],
                warning_reasons=[],
                key_metrics={},
                allowed_future_use="report_only",
                forbidden_interpretations=[],
            )
        ],
        stop_rules=[],
        forbidden_interpretations=[],
    )

    entries = build_v2_stress_test_plan(registry)

    assert not entries[0].allowed
