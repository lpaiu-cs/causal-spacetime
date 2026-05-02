"""Count v4 protocol carry-forward decision outcomes and failed reasons."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v4_protocol_carry_forward_experiment_helpers import data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 protocol accounting.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    rows = read_csv_rows(
        data_path(
            config.output_dir, "v4_protocol_cross_family_robustness_decisions.csv"
        )
    )
    counts: Counter[tuple[str, str]] = Counter()
    for row in rows:
        counts[("decision", row.get("decision", ""))] += 1
        for reason in row.get("failed_reasons", "").split(";"):
            if reason:
                counts[("failed_reason", reason)] += 1
    return [
        {"row_type": row_type, "label": label, "count": float(count)}
        for (row_type, label), count in sorted(counts.items())
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v4_protocol_failed_provisional_accounting.csv"),
        ["row_type", "label", "count"],
    )
    figure = save_metric_bar(
        [row for row in rows if row["row_type"] == "decision"],
        label_key="label",
        value_key="count",
        path=config.output_dir
        / "figures"
        / "v4_protocol_failed_provisional_counts.png",
        ylabel="family count",
    )
    print(f"Wrote v4 protocol accounting: {path}")
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
