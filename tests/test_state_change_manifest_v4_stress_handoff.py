from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)
from causal_spacetime_lab.state_change_manifest_v4_registry import (
    build_v4_protocol_carry_forward_registry,
)
from causal_spacetime_lab.state_change_manifest_v4_stress_handoff import (
    build_v4_protocol_stress_test_plan,
    evaluate_v4_protocol_stress_test_stop_condition,
)


def test_evaluate_v4_protocol_stress_test_stop_condition_stop_case() -> None:
    report = evaluate_v4_protocol_stress_test_stop_condition(
        [
            FamilyRobustnessDecision(
                "family",
                "structured",
                "blocked",
                False,
                ["high_heldout_violation"],
                [],
                {},
            )
        ]
    )

    assert not report.stress_tests_allowed


def test_build_v4_protocol_stress_test_plan_blocks_ineligible_families(
    m44_manifest_dir,
) -> None:
    registry = build_v4_protocol_carry_forward_registry(
        [
            FamilyRobustnessDecision(
                "rank_gap_earliest_full_stability_v4",
                "structured",
                "blocked",
                False,
                ["high_heldout_violation"],
                [],
                {},
            )
        ],
        m44_manifest_dir,
    )

    entries = build_v4_protocol_stress_test_plan(registry)

    assert entries
    assert not entries[0].allowed
