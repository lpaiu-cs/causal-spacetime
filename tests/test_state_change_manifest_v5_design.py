from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
)


def test_default_v5_protocol_family_designs_are_planned_only() -> None:
    designs = default_v5_protocol_family_designs()

    assert designs
    assert all(design.execution_status == "planned_only" for design in designs)
    assert all(
        design.planned_manifest_count >= 3
        for design in designs
        if design.family_kind == "structured"
    )
