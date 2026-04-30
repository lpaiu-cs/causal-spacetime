"""Audit actual v2 blocked decisions by blocking type."""

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
    summarize_v2_blocking_records,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 blocked root-cause audit."""

    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V2 blocked root-cause audit.")
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
    summary = summarize_v2_blocking_records(records)
    return [asdict(record) for record in records], summary


def write_outputs(
    rows: list[dict[str, float | str]],
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    return (
        write_csv(
            rows,
            data_path(output_dir, "v2_blocked_root_cause_audit.csv"),
            [
                "family_name",
                "family_kind",
                "criterion_name",
                "observed_value",
                "threshold_value",
                "criterion_direction",
                "blocking_type",
                "status",
                "explanation",
            ],
        ),
        write_csv(
            summary,
            data_path(output_dir, "v2_blocked_root_cause_summary.csv"),
            ["family_name", "family_kind", "blocking_type", "status", "count"],
        ),
    )


def save_figures(
    summary: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    counts: dict[str, float] = {}
    for row in summary:
        if row["status"] != "fail":
            continue
        key = str(row["blocking_type"])
        counts[key] = counts.get(key, 0.0) + float(row["count"])
    rows = [
        {"blocking_type": key, "count": value}
        for key, value in sorted(counts.items())
    ]
    return save_metric_bar(
        rows,
        label_key="blocking_type",
        value_key="count",
        path=output_dir / "figures" / "v2_blocking_type_counts.png",
        ylabel="failed criterion count",
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = write_outputs(rows, summary, config.output_dir)
    figure = save_figures(summary, config.output_dir)
    print("Wrote v2 blocked root-cause audit: " + ", ".join(str(p) for p in paths))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
