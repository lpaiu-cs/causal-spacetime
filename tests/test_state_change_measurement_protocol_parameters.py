from __future__ import annotations

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    default_earliest_full_protocol,
    default_gated_full_protocol,
    measurement_protocol_hash,
    measurement_protocol_id,
    parameter_complete_for_protocol,
)


def test_protocol_parameters_are_included_in_hash() -> None:
    base = default_gated_full_protocol(gate_delay_rank=2)
    changed = default_gated_full_protocol(gate_delay_rank=3)

    assert measurement_protocol_hash(base) != measurement_protocol_hash(changed)
    assert measurement_protocol_id(base) != measurement_protocol_id(changed)


def test_changing_margin_value_changes_protocol_id() -> None:
    base = default_earliest_full_protocol(margin_value=0.05)
    changed = default_earliest_full_protocol(margin_value=0.10)

    assert measurement_protocol_id(base) != measurement_protocol_id(changed)


def test_reference_chain_ids_do_not_change_protocol_id() -> None:
    protocol = default_earliest_full_protocol()
    reference_chain_ids_a = ["r1", "r2"]
    reference_chain_ids_b = ["x1", "x2", "x3"]

    assert reference_chain_ids_a != reference_chain_ids_b
    assert measurement_protocol_id(protocol) == measurement_protocol_id(protocol)


def test_parameter_incomplete_protocol_is_not_complete() -> None:
    incomplete = MeasurementProtocolSpec(
        echo_rule="gated_return",
        emission_rule="fixed_position",
        gate_rule="fixed_min_delay",
        reference_subsampling_rule="none",
        spectrum_type="full_transitive",
        normalization_rule="per_protocol_range",
        missing_policy="common_reachable",
        tie_policy="tie_as_unresolved",
        margin_policy="fixed_margin",
        protocol_parameters={"margin_value": 0.05},
    )

    assert not parameter_complete_for_protocol(incomplete)

