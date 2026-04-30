from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
    v4_protocol_family_design_table,
)


def test_default_v4_protocol_family_designs_returns_planned_only() -> None:
    designs = default_v4_protocol_family_designs()

    assert {design.execution_status for design in designs} == {"planned_only"}
    assert any(design.family_name == "failed_controls_v4" for design in designs)


def test_structured_v4_families_have_planned_count_at_least_three() -> None:
    designs = default_v4_protocol_family_designs()

    assert all(
        design.planned_manifest_count >= 3
        for design in designs
        if design.family_kind == "structured"
    )


def test_v4_protocol_family_design_table_serializes_root_causes() -> None:
    rows = v4_protocol_family_design_table(default_v4_protocol_family_designs())

    assert isinstance(rows[0]["linked_v3_root_causes"], str)
