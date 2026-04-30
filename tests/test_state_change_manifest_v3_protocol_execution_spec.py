from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    load_v3_protocol_execution_specs,
    mark_v3_protocol_execution_specs_executed,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    parameter_complete_for_protocol,
)


def test_load_v3_protocol_execution_specs_reads_m40_preregistration(
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    structured = [spec for spec in specs if spec.family_kind == "structured"]

    assert specs
    assert structured
    assert all(spec.admissible_for_pairwise_dissimilarity for spec in structured)
    assert all(
        parameter_complete_for_protocol(spec.measurement_protocol)
        for spec in structured
    )
    assert all(spec.planned_manifest_count >= 3 for spec in structured)


def test_mark_v3_protocol_execution_specs_executed(m40_prereg_path) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    executed = mark_v3_protocol_execution_specs_executed(specs)

    assert all(spec.execution_status == "executed" for spec in executed)
    assert all(spec.execution_status == "planned_only" for spec in specs)
