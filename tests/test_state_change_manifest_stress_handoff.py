from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    build_carry_forward_registry,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)
from causal_spacetime_lab.state_change_manifest_stress_handoff import (
    build_stress_test_plan,
    default_allowed_future_stress_tests,
)


def test_default_allowed_future_stress_tests_by_decision() -> None:
    assert "harder_nulls" in default_allowed_future_stress_tests("carry_forward")
    assert default_allowed_future_stress_tests("blocked") == []


def test_stress_test_plan_blocks_blocked_families() -> None:
    registry = build_carry_forward_registry(
        [
            FamilyRobustnessDecision(
                family_name="blocked_family",
                family_kind="structured",
                decision="blocked",
                passed=False,
                failed_reasons=["high_heldout_violation"],
                warning_reasons=[],
                key_metrics={},
            )
        ]
    )

    plan = build_stress_test_plan(registry)

    assert not plan[0].allowed
