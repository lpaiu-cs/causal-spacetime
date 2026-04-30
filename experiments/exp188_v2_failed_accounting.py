"""Report v2 failed, ineligible, and report-only controls."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    data_path,
    figure_path,
    save_metric_bar,
)

from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_failed_accounting_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 failed accounting."""

    manifest_dir: Path = Path("outputs/manifests_v2")
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 failed accounting.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(manifest_dir=args.manifest_dir, output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run v2 failed accounting."""

    return compute_v2_failed_accounting_rows(config.manifest_dir)


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 failed accounting rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_manifest_family_failed_accounting.csv"),
        ["row_type", "family_name", "family_kind", "reason", "count"],
    )


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 failed accounting figure."""

    reason_rows = [row for row in rows if row.get("row_type") == "failure_reason"]
    return [
        save_metric_bar(
            reason_rows,
            label_key="reason",
            value_key="count",
            path=figure_path(
                output_dir,
                "v2_manifest_family_failed_reason_counts.png",
            ),
            ylabel="failed reason count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote v2 failed accounting: {output_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
