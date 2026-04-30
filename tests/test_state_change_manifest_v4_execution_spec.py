from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    load_v4_protocol_execution_specs,
    mark_v4_protocol_execution_specs_executed,
)


def test_load_v4_protocol_execution_specs_reads_m43_preregistration(
    m43_v4_prereg_path,
) -> None:
    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)
    structured = [spec for spec in specs if spec.family_kind == "structured"]
    controls = [spec for spec in specs if spec.family_kind != "structured"]

    assert specs
    assert structured
    assert controls
    assert all(spec.planned_manifest_count >= 3 for spec in structured)
    assert all(spec.execution_status == "planned_only" for spec in specs)


def test_mark_v4_protocol_execution_specs_executed(m43_v4_prereg_path) -> None:
    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)
    executed = mark_v4_protocol_execution_specs_executed(specs)

    assert all(spec.execution_status == "planned_only" for spec in specs)
    assert all(spec.execution_status == "executed" for spec in executed)
