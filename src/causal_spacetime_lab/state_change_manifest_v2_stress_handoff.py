"""V2 future stress-test handoff utilities."""

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


def evaluate_v2_stress_test_stop_condition(
    decisions: list[FamilyRobustnessDecision],
) -> StopConditionReport:
    """Evaluate the v2 stress-test stop condition from family decisions."""

    return evaluate_stress_test_stop_condition(
        [
            {
                "family_name": decision.family_name,
                "decision": decision.decision,
            }
            for decision in decisions
        ]
    )


def build_v2_stress_test_plan(
    registry: CarryForwardRegistry,
) -> list[StressTestPlanEntry]:
    """Build a v2 future stress-test plan from a carry-forward registry."""

    return build_stress_test_plan(registry)


def v2_stress_test_plan_rows(
    entries: list[StressTestPlanEntry],
) -> list[dict[str, float | str]]:
    """Convert v2 stress-test handoff entries to CSV rows."""

    return stress_test_plan_table(entries)
