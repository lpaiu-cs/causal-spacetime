"""Export planned-only v5 protocol preregistration JSON."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_blocked_v5_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v5_preregistration import (
    build_v5_protocol_preregistration_spec,
    write_v5_protocol_preregistration_spec,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V5 preregistration export.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    spec = build_v5_protocol_preregistration_spec(
        default_v5_protocol_family_designs()
    )
    json_path = write_v5_protocol_preregistration_spec(
        spec,
        config.output_dir / "remediation" / "v5_protocol_preregistration_spec_m46.json",
    )
    return [
        {
            "spec_id": spec.spec_id,
            "created_by_milestone": spec.created_by_milestone,
            "planned_family_count": float(len(spec.planned_families)),
            "execution_allowed_in_current_milestone": float(
                spec.execution_allowed_in_current_milestone
            ),
            "remediation_iteration_risk_audit_required": float(
                spec.remediation_iteration_risk_audit_required
            ),
            "json_path": str(json_path),
        }
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v5_protocol_preregistration_export.csv"),
        ["spec_id", "json_path"],
    )
    print(f"Wrote v5 preregistration export: {path}")
    print(rows[0]["json_path"])


if __name__ == "__main__":
    main()
