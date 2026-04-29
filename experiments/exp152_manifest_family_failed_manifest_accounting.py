"""Failed and ineligible manifest accounting for family reports."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import (
    assignments_for_datasets,
    ensure_failed_control_manifest,
    load_or_create_family_datasets,
    save_bar_figure,
    write_csv,
)

from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    failed_manifest_accounting_summary,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for failed-manifest accounting."""

    manifest_dir: Path = Path("outputs/manifests")
    generate_failed_controls: bool = True
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Manifest-family failed accounting.")
    parser.add_argument("--manifest-dir", type=Path, default=Path("outputs/manifests"))
    parser.add_argument(
        "--generate-failed-controls",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        manifest_dir=args.manifest_dir,
        generate_failed_controls=args.generate_failed_controls,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run failed/ineligible manifest accounting."""

    datasets = load_or_create_family_datasets(config.manifest_dir, config.output_dir)
    if config.generate_failed_controls and not any(
        not dataset.eligible for dataset in datasets
    ):
        datasets.append(ensure_failed_control_manifest(config.output_dir))
    assignments = assignments_for_datasets(datasets)
    rows = failed_manifest_accounting_summary(assignments)
    rows.append(
        {
            "row_type": "eligibility_summary",
            "family_name": "",
            "family_kind": "",
            "reason": "eligible_count",
            "count": float(sum(assignment.eligible for assignment in assignments)),
        }
    )
    rows.append(
        {
            "row_type": "eligibility_summary",
            "family_name": "",
            "family_kind": "",
            "reason": "ineligible_count",
            "count": float(sum(not assignment.eligible for assignment in assignments)),
        }
    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write failed-manifest accounting CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "manifest_family_failed_manifest_accounting.csv",
        ["row_type", "family_name", "family_kind", "reason", "count"],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save failed reason count figure."""

    reason_rows = [row for row in rows if row["row_type"] == "failure_reason"]
    labels = [str(row["reason"]) for row in reason_rows]
    values = [float(row["count"]) for row in reason_rows]
    return [
        save_bar_figure(
            labels,
            values,
            output_dir / "figures" / "manifest_family_failed_reason_counts.png",
            ylabel="failed manifest count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote manifest family failed-manifest accounting: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
