from __future__ import annotations

import json
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_registry import (
    build_v3_protocol_carry_forward_registry,
    write_v3_protocol_carry_forward_registry,
)


def test_build_v3_protocol_carry_forward_registry_serializes(
    tmp_path: Path,
    m41_manifest_dir: Path,
) -> None:
    registry = build_v3_protocol_carry_forward_registry(
        [
            FamilyRobustnessDecision(
                family_name="rank_gap_earliest_full_reference_v3",
                family_kind="structured",
                decision="blocked",
                passed=False,
                failed_reasons=["high_heldout_violation"],
                warning_reasons=[],
                key_metrics={"manifest_count": 2.0},
            )
        ],
        m41_manifest_dir,
    )
    path = write_v3_protocol_carry_forward_registry(registry, tmp_path)

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["created_by_milestone"] == "Milestone 42"

