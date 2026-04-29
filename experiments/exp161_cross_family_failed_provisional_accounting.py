"""Audit blocked, provisional, failed, and ineligible family accounting."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_bar_figure

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for failed/provisional accounting."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Cross-family accounting audit.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _read_decisions(output_dir: Path) -> list[dict[str, str]]:
    path = output_dir / "data" / "cross_family_robustness_decisions.csv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Count decision types and failed reasons."""

    decisions = _read_decisions(config.output_dir)
    decision_counts = Counter(row.get("decision", "") for row in decisions)
    rows: list[dict[str, float | str]] = [
        {
            "row_type": "decision_count",
            "name": decision,
            "count": float(count),
        }
        for decision, count in sorted(decision_counts.items())
    ]
    reason_counts: Counter[str] = Counter()
    for row in decisions:
        for reason in row.get("failed_reasons", "").split(";"):
            if reason:
                reason_counts[reason] += 1
    rows.extend(
        {
            "row_type": "failed_reason",
            "name": reason,
            "count": float(count),
        }
        for reason, count in sorted(reason_counts.items())
    )
    omitted = 0.0 if decisions else 1.0
    rows.append(
        {
            "row_type": "omission_check",
            "name": "blocked_provisional_failed_visible",
            "count": 1.0 - omitted,
        }
    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write accounting audit CSV."""

    path = output_dir / "data" / "cross_family_failed_provisional_accounting.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["row_type", "name", "count"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save decision count figure."""

    decision_rows = [row for row in rows if row["row_type"] == "decision_count"]
    return [
        save_bar_figure(
            [str(row["name"]) for row in decision_rows],
            [float(row["count"]) for row in decision_rows],
            output_dir / "figures" / "cross_family_failed_provisional_counts.png",
            ylabel="family count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    figures = save_figures(rows, config.output_dir)
    print(f"Wrote cross-family failed/provisional accounting: {path}")
    for figure in figures:
        print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
