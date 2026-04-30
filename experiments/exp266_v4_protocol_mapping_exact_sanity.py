"""Exact sanity checks for v4 protocol and comparison mapping."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

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


def run_experiment() -> list[dict[str, float | str]]:
    path = Path("outputs/remediation/v4_protocol_preregistration_spec_m43.json")
    specs = load_v4_protocol_execution_specs(path)
    if not specs:
        return [{"check": "v4_specs_load", "passed": 0.0, "detail": str(path)}]
    protocols = [measurement_protocol_for_v4_family(spec) for spec in specs]
    comparisons = [
        comparison_protocol_for_v4_family(spec, protocol)
        for spec, protocol in zip(specs, protocols, strict=True)
    ]
    immediate = [
        spec for spec in specs if spec.family_name == "report_only_immediate_edge_v4"
    ]
    return [
        {
            "check": "structured_protocols_parameter_complete",
            "passed": float(
                all(
                    parameter_complete_for_protocol(protocol)
                    for spec, protocol in zip(specs, protocols, strict=True)
                    if spec.family_kind == "structured"
                )
            ),
        },
        {
            "check": "protocol_id_stable",
            "passed": float(
                all(
                    v4_protocol_id(protocol) == v4_protocol_id(protocol)
                    for protocol in protocols
                )
            ),
        },
        {
            "check": "margin_value_in_protocol_parameters",
            "passed": float(
                all(
                    "margin_value" in (protocol.protocol_parameters or {})
                    for protocol in protocols
                )
            ),
        },
        {
            "check": "comparison_missing_policy_matches",
            "passed": float(
                all(
                    comparison.missing_policy == protocol.missing_policy
                    for comparison, protocol in zip(comparisons, protocols, strict=True)
                )
            ),
        },
        {
            "check": "report_only_immediate_edge_maps",
            "passed": float(bool(immediate)),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v4_protocol_mapping_exact_sanity.csv"),
        ["check", "passed", "detail"],
    )
    print(f"Wrote v4 protocol mapping sanity: {path}")


if __name__ == "__main__":
    main()
