"""Exact sanity checks for M41 v3 protocol execution specs."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    load_v3_protocol_execution_specs,
    mark_v3_protocol_execution_specs_executed,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    parameter_complete_for_protocol,
)


def run_experiment(
    output_dir: Path = Path("outputs"),
) -> list[dict[str, float | str]]:
    path = output_dir / "remediation" / "v3_protocol_patched_preregistration_m40.json"
    specs = load_v3_protocol_execution_specs(path)
    structured = [spec for spec in specs if spec.family_kind == "structured"]
    executed = mark_v3_protocol_execution_specs_executed(specs)
    return [
        {"check": "patched_preregistration_exists", "passed": float(path.exists())},
        {
            "check": "structured_patched_families_load",
            "passed": float(bool(structured)),
        },
        {
            "check": "structured_families_admissible",
            "passed": float(
                bool(structured)
                and all(
                    spec.admissible_for_pairwise_dissimilarity for spec in structured
                )
            ),
        },
        {
            "check": "structured_families_parameter_complete",
            "passed": float(
                bool(structured)
                and all(
                    parameter_complete_for_protocol(spec.measurement_protocol)
                    for spec in structured
                )
            ),
        },
        {
            "check": "planned_manifest_count_at_least_3",
            "passed": float(
                bool(structured)
                and all(spec.planned_manifest_count >= 3 for spec in structured)
            ),
        },
        {
            "check": "handoff_provenance_type_exists",
            "passed": float(all(bool(spec.handoff_provenance_type) for spec in specs)),
        },
        {
            "check": "mark_executed_works",
            "passed": float(
                bool(executed)
                and all(spec.execution_status == "executed" for spec in executed)
            ),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v3_protocol_execution_spec_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v3 protocol execution spec sanity: {path}")


if __name__ == "__main__":
    main()
