"""Export planned-only v4 protocol preregistration JSON."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_blocked_v4_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v4_preregistration import (
    build_v4_protocol_preregistration_spec,
    write_v4_protocol_preregistration_spec,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Export v4 preregistration.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[Path, list[dict[str, float | str]]]:
    spec = build_v4_protocol_preregistration_spec(
        default_v4_protocol_family_designs(),
    )
    path = write_v4_protocol_preregistration_spec(
        spec,
        config.output_dir / "remediation" / "v4_protocol_preregistration_spec_m43.json",
    )
    rows = [
        {
            "spec_id": spec.spec_id,
            "created_by_milestone": spec.created_by_milestone,
            "planned_family_count": float(len(spec.planned_families)),
            "execution_allowed_in_current_milestone": float(
                spec.execution_allowed_in_current_milestone
            ),
            "output_path": str(path),
        }
    ]
    return path, rows


def main() -> None:
    config = parse_args()
    spec_path, rows = run_experiment(config)
    summary_path = write_csv(
        rows,
        data_path(config.output_dir, "v4_protocol_preregistration_export.csv"),
        ["spec_id", "output_path"],
    )
    print(f"Wrote v4 preregistration export: {summary_path}")
    print(f"Wrote v4 preregistration JSON: {spec_path}")


if __name__ == "__main__":
    main()
