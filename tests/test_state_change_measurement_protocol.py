from __future__ import annotations

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    default_earliest_full_protocol,
    measurement_protocol_hash,
)


def test_identical_measurement_protocol_specs_have_same_hash() -> None:
    assert measurement_protocol_hash(
        default_earliest_full_protocol()
    ) == measurement_protocol_hash(default_earliest_full_protocol())


def test_reference_chain_ids_do_not_change_protocol_hash() -> None:
    protocol = default_earliest_full_protocol()
    reference_chain_ids_a = ["r1", "r2"]
    reference_chain_ids_b = ["x", "y", "z"]

    assert reference_chain_ids_a != reference_chain_ids_b
    assert measurement_protocol_hash(protocol) == measurement_protocol_hash(protocol)


def test_changing_measurement_rules_changes_protocol_hash() -> None:
    base = default_earliest_full_protocol()
    changed_specs = [
        MeasurementProtocolSpec(
            **{
                **base.__dict__,
                "echo_rule": "gated_return",
                "gate_rule": "fixed_min_delay",
            }
        ),
        MeasurementProtocolSpec(**{**base.__dict__, "gate_rule": "declared_gate"}),
        MeasurementProtocolSpec(
            **{
                **base.__dict__,
                "spectrum_type": "retained_reference",
                "reference_subsampling_rule": "retained_reference_set",
            }
        ),
        MeasurementProtocolSpec(
            **{**base.__dict__, "missing_policy": "penalize_mismatch"}
        ),
    ]

    base_hash = measurement_protocol_hash(base)
    assert all(measurement_protocol_hash(item) != base_hash for item in changed_specs)
