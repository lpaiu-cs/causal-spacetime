"""Report fixed-criteria margins for v4 protocol families."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v4_blocked_v5_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_v4_blocking_analysis import (
    V4ProtocolBlockingRecord,
)
from causal_spacetime_lab.state_change_manifest_v4_criterion_margins import (
    criterion_margin_to_row,
    criterion_margins_from_v4_blocking_records,
    summarize_v4_criterion_margins,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 criterion margins.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def _records_from_rows(rows: list[dict[str, str]]) -> list[V4ProtocolBlockingRecord]:
    records: list[V4ProtocolBlockingRecord] = []
    for row in rows:
        if not row.get("family_name"):
            continue
        records.append(
            V4ProtocolBlockingRecord(
                family_name=row["family_name"],
                family_kind=row["family_kind"],
                criterion_name=row["criterion_name"],
                observed_value=float(row["observed_value"]),
                threshold_value=float(row["threshold_value"]),
                criterion_direction=row["criterion_direction"],
                blocking_type=row["blocking_type"],
                status=row["status"],
                root_cause_category=row["root_cause_category"],
                explanation=row["explanation"],
            )
        )
    return records


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    audit_rows = read_csv_rows(
        data_path(config.output_dir, "v4_blocked_root_cause_audit.csv")
    )
    margins = criterion_margins_from_v4_blocking_records(
        _records_from_rows(audit_rows)
    )
    return (
        [criterion_margin_to_row(margin) for margin in margins],
        summarize_v4_criterion_margins(margins),
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = [
        write_csv(
            rows,
            data_path(config.output_dir, "v4_criterion_margin_report.csv"),
            ["family_name", "criterion_name", "margin"],
        ),
        write_csv(
            summary,
            data_path(config.output_dir, "v4_criterion_margin_summary.csv"),
            ["family_name", "worst_criterion", "worst_margin"],
        ),
    ]
    figure = save_metric_bar(
        summary,
        label_key="family_name",
        value_key="worst_margin",
        path=config.output_dir / "figures" / "v4_worst_criterion_margins.png",
        ylabel="worst signed margin",
    )
    print("Wrote v4 criterion margins: " + ", ".join(map(str, paths)))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
