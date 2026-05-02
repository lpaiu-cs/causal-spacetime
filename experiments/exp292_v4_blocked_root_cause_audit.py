"""Audit actual M45 v4 protocol blocked decisions by root cause."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v4_blocked_v5_experiment_helpers import (
    data_path,
    missing_input_row,
    read_csv_rows,
    to_float_rows,
)

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_blocking_analysis import (
    blocking_record_to_row,
    decompose_v4_blocking_by_family,
    summarize_v4_blocking_records,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 root-cause audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    metrics_path = data_path(
        config.output_dir,
        "v4_protocol_cross_family_robustness_decision_metrics.csv",
    )
    decisions_path = data_path(
        config.output_dir,
        "v4_protocol_cross_family_robustness_decisions.csv",
    )
    preconditions_path = data_path(
        config.output_dir,
        "v4_protocol_precondition_audit.csv",
    )
    if not metrics_path.exists() or not decisions_path.exists():
        missing = [
            missing_input_row(str(path))
            for path in [metrics_path, decisions_path]
            if not path.exists()
        ]
        return missing, missing
    records = decompose_v4_blocking_by_family(
        to_float_rows(read_csv_rows(metrics_path)),
        to_float_rows(read_csv_rows(decisions_path)),
        to_float_rows(read_csv_rows(preconditions_path)),
        default_cross_family_robustness_criteria(),
    )
    return (
        [blocking_record_to_row(record) for record in records],
        summarize_v4_blocking_records(records),
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = [
        write_csv(
            rows,
            data_path(config.output_dir, "v4_blocked_root_cause_audit.csv"),
            ["family_name", "criterion_name", "status"],
        ),
        write_csv(
            summary,
            data_path(config.output_dir, "v4_blocked_root_cause_summary.csv"),
            ["family_name", "blocking_type", "root_cause_category", "status", "count"],
        ),
    ]
    counts: dict[str, float] = {}
    for row in summary:
        if row.get("status") == "fail":
            root = str(row.get("root_cause_category", ""))
            counts[root] = counts.get(root, 0.0) + float(row.get("count", 0.0))
    figure = save_metric_bar(
        [
            {"root_cause_category": root, "count": count}
            for root, count in sorted(counts.items())
        ],
        label_key="root_cause_category",
        value_key="count",
        path=config.output_dir / "figures" / "v4_blocking_root_cause_counts.png",
        ylabel="failure count",
    )
    print("Wrote v4 root-cause audit: " + ", ".join(map(str, paths)))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
