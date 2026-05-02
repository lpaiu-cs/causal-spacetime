"""Audit v5 design proposals for remediation-iteration risk."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_blocked_v5_experiment_helpers import (
    data_path,
    read_csv_rows,
    to_float_rows,
)

from causal_spacetime_lab.state_change_manifest_remediation_iteration_risk import (
    audit_v5_remediation_iteration_risk,
    remediation_iteration_risk_table,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V5 remediation risk audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    records = audit_v5_remediation_iteration_risk(
        to_float_rows(
            read_csv_rows(
                data_path(config.output_dir, "v4_blocked_root_cause_summary.csv")
            )
        ),
        to_float_rows(
            read_csv_rows(
                data_path(config.output_dir, "v3_to_v4_metric_delta_summary.csv")
            )
        ),
        to_float_rows(
            read_csv_rows(data_path(config.output_dir, "v5_protocol_family_design.csv"))
        ),
    )
    return remediation_iteration_risk_table(records)


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v5_remediation_iteration_risk_audit.csv"),
        ["proposed_action", "risk_category", "risk_level"],
    )
    print(f"Wrote v5 remediation iteration risk audit: {path}")


if __name__ == "__main__":
    main()
