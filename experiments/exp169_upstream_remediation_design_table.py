"""Map carry-forward failure root causes to report-only remediation proposals."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import (
    data_path,
    failure_record_from_row,
    read_csv_rows,
)
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_remediation_design import (
    default_remediation_proposals,
    propose_remediations_for_failure_records,
    remediation_proposal_table,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for remediation design table."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Upstream remediation design table.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Generate remediation proposal rows for observed failures."""

    failure_rows = read_csv_rows(
        data_path(config.output_dir, "carry_forward_failure_decomposition.csv")
    )
    records = [failure_record_from_row(row) for row in failure_rows]
    proposals = (
        propose_remediations_for_failure_records(records)
        if records
        else default_remediation_proposals()
    )
    rows = remediation_proposal_table(proposals)
    for row in rows:
        row["milestone35_action"] = "report_only_design"
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write remediation design table."""

    return write_csv(
        rows,
        output_dir / "data" / "upstream_remediation_design_table.csv",
        [
            "proposal_name",
            "target_root_cause",
            "description",
            "expected_effect",
            "requires_new_preregistration",
            "allowed_in_current_milestone",
            "milestone35_action",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote upstream remediation design table: {path}")


if __name__ == "__main__":
    main()
