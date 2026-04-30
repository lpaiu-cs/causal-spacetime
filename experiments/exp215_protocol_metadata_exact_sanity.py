"""Exact sanity checks for measurement-protocol metadata."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    default_earliest_full_protocol,
    measurement_protocol_id,
    same_measurement_protocol,
)


def run_experiment() -> list[dict[str, float | str]]:
    base = default_earliest_full_protocol()
    same = default_earliest_full_protocol()
    echo_changed = MeasurementProtocolSpec(
        **{**base.__dict__, "echo_rule": "gated_return", "gate_rule": "fixed_min_delay"}
    )
    gate_changed = MeasurementProtocolSpec(
        **{**base.__dict__, "gate_rule": "declared_gate"}
    )
    spectrum_changed = MeasurementProtocolSpec(
        **{
            **base.__dict__,
            "spectrum_type": "retained_reference",
            "reference_subsampling_rule": "retained_reference_set",
        }
    )
    missing_changed = MeasurementProtocolSpec(
        **{**base.__dict__, "missing_policy": "penalize_mismatch"}
    )
    base_id = measurement_protocol_id(base)
    reference_id_a = measurement_protocol_id(base)
    reference_id_b = measurement_protocol_id(base)
    return [
        {
            "check": "identical_specs_same_hash",
            "passed": float(same_measurement_protocol(base, same)),
        },
        {
            "check": "reference_ids_do_not_change_protocol_id",
            "passed": float(reference_id_a == reference_id_b == base_id),
        },
        {
            "check": "changing_echo_rule_changes_id",
            "passed": float(measurement_protocol_id(echo_changed) != base_id),
        },
        {
            "check": "changing_gate_rule_changes_id",
            "passed": float(measurement_protocol_id(gate_changed) != base_id),
        },
        {
            "check": "changing_spectrum_type_changes_id",
            "passed": float(measurement_protocol_id(spectrum_changed) != base_id),
        },
        {
            "check": "changing_missing_policy_changes_id",
            "passed": float(measurement_protocol_id(missing_changed) != base_id),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/protocol_metadata_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote protocol metadata exact sanity: {path}")


if __name__ == "__main__":
    main()
