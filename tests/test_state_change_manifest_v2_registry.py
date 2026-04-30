from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)
from causal_spacetime_lab.state_change_manifest_v2_registry import (
    build_v2_carry_forward_registry,
    v2_registry_summary_rows,
)


def test_build_v2_carry_forward_registry_serializes(tmp_path: Path) -> None:
    decision = FamilyRobustnessDecision(
        family_name="family",
        family_kind="structured",
        decision="blocked",
        passed=False,
        failed_reasons=["low_manifest_count"],
        warning_reasons=[],
        key_metrics={"manifest_count": 1.0},
    )

    registry = build_v2_carry_forward_registry([decision], tmp_path)
    rows = v2_registry_summary_rows(registry)

    assert registry.created_by_milestone == "Milestone 38"
    assert rows[0]["decision"] == "blocked"
