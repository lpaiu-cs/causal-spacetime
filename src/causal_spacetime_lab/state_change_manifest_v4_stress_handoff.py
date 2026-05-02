"""V4 protocol future stress-test handoff utilities."""

from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    CarryForwardRegistry,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)
from causal_spacetime_lab.state_change_manifest_stop_condition import (
    StopConditionReport,
    evaluate_stress_test_stop_condition,
)
from causal_spacetime_lab.state_change_manifest_stress_handoff import (
    StressTestPlanEntry,
    build_stress_test_plan,
    stress_test_plan_table,
)


def evaluate_v4_protocol_stress_test_stop_condition(
    decisions: list[FamilyRobustnessDecision],
) -> StopConditionReport:
    """Evaluate the v4 protocol future stress-test stop condition."""

    return evaluate_stress_test_stop_condition(
        [
            {
                "family_name": decision.family_name,
                "decision": decision.decision,
            }
            for decision in decisions
        ]
    )


def build_v4_protocol_stress_test_plan(
    registry: CarryForwardRegistry,
) -> list[StressTestPlanEntry]:
    """Build a future stress-test handoff plan from the v4 protocol registry."""

    return build_stress_test_plan(registry)


def v4_protocol_stress_test_plan_rows(
    entries: list[StressTestPlanEntry],
) -> list[dict[str, float | str]]:
    """Convert v4 protocol stress-test handoff entries to CSV rows."""

    return stress_test_plan_table(entries)
