"""Aggregate manifest-family comparison outputs into a report card."""

from __future__ import annotations

import csv
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_family_comparison import (
    family_report_card,
    summarize_nulls_by_taxonomy,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Build a family-level report card from available outputs."""

    data_dir = output_dir / "data"
    family_summary = _read_rows(data_dir / "manifest_family_fit_summary.csv")
    null_rows = _read_rows(data_dir / "manifest_family_null_taxonomy.csv")
    failure_rows = _read_rows(
        data_dir / "manifest_family_failed_manifest_accounting.csv"
    )
    null_taxonomy = summarize_nulls_by_taxonomy(null_rows)
    return family_report_card(family_summary, null_taxonomy, failure_rows)


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write family report-card CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "manifest_family_report_card.csv",
        [
            "family_name",
            "family_kind",
            "fitted_count",
            "no_fit_count",
            "best_heldout_violation",
            "mean_heldout_violation",
            "null_comparison_summary",
            "failure_reason_summary",
            "interpretation_warning",
        ],
    )


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote manifest family report card: {output_path}")


if __name__ == "__main__":
    main()
