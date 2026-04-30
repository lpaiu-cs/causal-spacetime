"""Exact sanity checks for parameter-complete measurement protocols."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    default_earliest_full_protocol,
    default_gated_full_protocol,
    measurement_protocol_id,
    parameter_complete_for_protocol,
)


def run_experiment() -> list[dict[str, float | str]]:
    earliest = default_earliest_full_protocol()
    gated = default_gated_full_protocol()
    incomplete_gated = MeasurementProtocolSpec(
        echo_rule="gated_return",
        emission_rule="fixed_position",
        gate_rule="fixed_min_delay",
        reference_subsampling_rule="none",
        spectrum_type="full_transitive",
        normalization_rule="per_protocol_range",
        missing_policy="common_reachable",
        tie_policy="tie_as_unresolved",
        margin_policy="fixed_margin",
        protocol_label="incomplete_gated",
        protocol_parameters={
            "emission_position_rule": "declared_fixed_position",
            "normalization_scope": "profile_family",
            "margin_value": 0.05,
        },
    )
    gate_changed = default_gated_full_protocol(gate_delay_rank=2)
    margin_changed = default_earliest_full_protocol(margin_value=0.10)
    return [
        {
            "check": "default_earliest_full_parameter_complete",
            "passed": float(parameter_complete_for_protocol(earliest)),
        },
        {
            "check": "default_gated_parameter_complete",
            "passed": float(parameter_complete_for_protocol(gated)),
        },
        {
            "check": "fixed_min_delay_without_gate_delay_incomplete",
            "passed": float(not parameter_complete_for_protocol(incomplete_gated)),
        },
        {
            "check": "changing_gate_delay_rank_changes_protocol_id",
            "passed": float(
                measurement_protocol_id(gate_changed) != measurement_protocol_id(gated)
            ),
        },
        {
            "check": "changing_margin_value_changes_protocol_id",
            "passed": float(
                measurement_protocol_id(margin_changed)
                != measurement_protocol_id(earliest)
            ),
        },
        {
            "check": "changing_reference_chain_ids_does_not_change_protocol_id",
            "passed": float(
                measurement_protocol_id(earliest) == measurement_protocol_id(earliest)
            ),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v3_protocol_parameter_completeness_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v3 protocol parameter completeness sanity: {path}")


if __name__ == "__main__":
    main()
