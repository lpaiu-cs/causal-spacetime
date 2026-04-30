"""Count v2 blocked, provisional, failed, and report-only families."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, decisions_from_csv
from v2_manifest_experiment_helpers import save_metric_bar

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 accounting."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 failed/provisional accounting.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build v2 accounting rows."""

    decisions = decisions_from_csv(
        data_path(config.output_dir, "v2_cross_family_robustness_decisions.csv")
    )
    rows: list[dict[str, float | str]] = []
    decision_counts = Counter(decision.decision for decision in decisions)
    for decision, count in sorted(decision_counts.items()):
        rows.append(
            {
                "row_type": "decision_count",
                "decision": decision,
                "reason": "",
                "count": float(count),
            }
        )
    reason_counts: Counter[str] = Counter(
        reason
        for decision in decisions
        for reason in decision.failed_reasons
    )
    for reason, count in sorted(reason_counts.items()):
        rows.append(
            {
                "row_type": "failed_reason_count",
                "decision": "",
                "reason": reason,
                "count": float(count),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write v2 accounting rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_failed_provisional_accounting.csv"),
        ["row_type", "decision", "reason", "count"],
    )


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save v2 accounting figure."""

    return [
        save_metric_bar(
            [row for row in rows if row["row_type"] == "decision_count"],
            label_key="decision",
            value_key="count",
            path=output_dir / "figures" / "v2_failed_provisional_counts.png",
            ylabel="family count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote v2 failed/provisional accounting: {output_path}")
    for figure in figure_paths:
        print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
