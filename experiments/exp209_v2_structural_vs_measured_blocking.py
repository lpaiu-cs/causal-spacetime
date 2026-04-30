"""Report-only v2 structural versus measured blocking counterfactuals."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    decompose_v2_blocking_by_family,
)
from causal_spacetime_lab.state_change_manifest_v2_counterfactuals import (
    structural_count_counterfactual_report,
    would_remain_blocked_ignoring_structural_count,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for report-only structural counterfactuals."""

    output_dir: Path = Path("outputs")
    hypothetical_manifest_count: int = 3


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V2 structural counterfactual.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--hypothetical-manifest-count", type=int, default=3)
    args = parser.parse_args()
    return ExperimentConfig(
        output_dir=args.output_dir,
        hypothetical_manifest_count=args.hypothetical_manifest_count,
    )


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
    metrics = _metric_rows(config.output_dir)
    criteria = default_cross_family_robustness_criteria()
    records = decompose_v2_blocking_by_family(metrics, criteria)
    return (
        structural_count_counterfactual_report(
            metrics,
            criteria,
            hypothetical_manifest_count=config.hypothetical_manifest_count,
        ),
        would_remain_blocked_ignoring_structural_count(records),
    )


def write_outputs(
    count_rows: list[dict[str, float | str]],
    remaining_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    return (
        write_csv(
            count_rows,
            data_path(output_dir, "v2_structural_count_counterfactual.csv"),
            [
                "family_name",
                "family_kind",
                "original_manifest_count",
                "hypothetical_manifest_count",
                "manifest_count_would_pass",
                "remaining_measured_blocker_count",
                "remaining_measured_blockers",
                "analysis_label",
            ],
        ),
        write_csv(
            remaining_rows,
            data_path(
                output_dir,
                "v2_remaining_measured_blockers_after_structural_ignore.csv",
            ),
            [
                "family_name",
                "family_kind",
                "ignored_structural_count_failure",
                "measured_blocker_count",
                "measured_blockers",
                "would_remain_blocked",
                "analysis_label",
            ],
        ),
    )


def main() -> None:
    config = parse_args()
    count_rows, remaining_rows = run_experiment(config)
    paths = write_outputs(count_rows, remaining_rows, config.output_dir)
    print("Wrote v2 structural counterfactuals: " + ", ".join(str(p) for p in paths))


if __name__ == "__main__":
    main()
