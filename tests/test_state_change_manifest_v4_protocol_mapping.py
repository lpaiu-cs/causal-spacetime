from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    load_v4_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v4_protocol_mapping import (
    comparison_protocol_for_v4_family,
    measurement_protocol_for_v4_family,
    v4_protocol_id,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    parameter_complete_for_protocol,
)


def test_measurement_protocol_for_v4_family_is_parameter_complete(
    m43_v4_prereg_path,
) -> None:
    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)

    for spec in specs:
        protocol = measurement_protocol_for_v4_family(spec)
        params = protocol.protocol_parameters or {}
        comparison = comparison_protocol_for_v4_family(spec, protocol)

        assert parameter_complete_for_protocol(protocol)
        assert float(params["margin_value"]) == spec.margin_value
        assert comparison.missing_policy == protocol.missing_policy


def test_v4_protocol_id_is_stable_for_identical_parameters(
    m43_v4_prereg_path,
) -> None:
    spec = load_v4_protocol_execution_specs(m43_v4_prereg_path)[0]
    protocol_a = measurement_protocol_for_v4_family(spec)
    protocol_b = measurement_protocol_for_v4_family(spec)

    assert v4_protocol_id(protocol_a) == v4_protocol_id(protocol_b)
