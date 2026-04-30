from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_design import (
    default_v3_manifest_family_designs,
)


def test_default_v3_structured_families_have_replicates() -> None:
    designs = default_v3_manifest_family_designs()

    assert any(design.family_name == "rank_gap_multi_manifest_v3" for design in designs)
    assert all(
        design.planned_manifest_count >= 3
        for design in designs
        if design.family_kind == "structured"
    )
    assert all(design.execution_status == "planned_only" for design in designs)

