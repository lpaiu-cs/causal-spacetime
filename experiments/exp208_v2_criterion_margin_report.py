"""Report fixed-criterion margins for v2 families."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows
from v2_manifest_experiment_helpers import save_metric_bar

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    decompose_v2_blocking_by_family,
)
from causal_spacetime_lab.state_change_manifest_v2_criterion_margins import (
    criterion_margins_from_blocking_records,
    summarize_criterion_margins,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for criterion-margin report."""

    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V2 criterion-margin report.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _metric_rows(output_dir: Path) -> list[dict[str, float | str]]:
    rows = read_csv_rows(
        data_path(output_dir, "v2_cross_family_robustness_decision_metrics.csv")
    )
    if not rows:
        rows = read_csv_rows(
            data_path(output_dir, "v2_cross_family_robustness_metrics.csv")
        )
    return [{key: value for key, value in row.items()} for row in rows]


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    records = decompose_v2_blocking_by_family(
        _metric_rows(config.output_dir),
        default_cross_family_robustness_criteria(),
    )
    margins = criterion_margins_from_blocking_records(records)
    return [asdict(margin) for margin in margins], summarize_criterion_margins(margins)


def write_outputs(
    rows: list[dict[str, float | str]],
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    return (
        write_csv(
            rows,
            data_path(output_dir, "v2_criterion_margin_report.csv"),
            [
                "family_name",
                "criterion_name",
                "observed_value",
                "threshold_value",
                "margin",
                "normalized_margin",
                "passed",
                "blocking_type",
            ],
        ),
        write_csv(
            summary,
            data_path(output_dir, "v2_criterion_margin_summary.csv"),
            [
                "family_name",
                "worst_margin",
                "worst_normalized_margin",
                "worst_criterion",
                "failed_criterion_count",
                "structural_failure_count",
                "measured_failure_count",
                "diagnostic_failure_count",
            ],
        ),
    )


def save_figures(
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    rows = [
        {
            "family_name": row["family_name"],
            "worst_normalized_margin": row["worst_normalized_margin"],
        }
        for row in summary
    ]
    return save_metric_bar(
        rows,
        label_key="family_name",
        value_key="worst_normalized_margin",
        path=output_dir / "figures" / "v2_worst_criterion_margins.png",
        ylabel="worst normalized margin",
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = write_outputs(rows, summary, config.output_dir)
    figure = save_figures(summary, config.output_dir)
    print("Wrote v2 criterion margins: " + ", ".join(str(p) for p in paths))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
