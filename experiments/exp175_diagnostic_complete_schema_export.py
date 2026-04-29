"""Export diagnostic-complete metric requirements."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    diagnostic_metric_requirement_table,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment() -> tuple[
    list[dict[str, float | str]],
    list[dict[str, float | str]],
]:
    """Return schema rows and grouped summary rows."""

    rows = diagnostic_metric_requirement_table()
    counts = Counter(
        (str(row["criterion_type"]), str(row["missing_status"]))
        for row in rows
    )
    summary = [
        {
            "criterion_type": criterion_type,
            "missing_status": missing_status,
            "metric_count": float(count),
        }
        for (criterion_type, missing_status), count in sorted(counts.items())
    ]
    return rows, summary


def write_outputs(
    rows: list[dict[str, float | str]],
    summary: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> tuple[Path, Path]:
    """Write diagnostic-complete schema CSVs."""

    schema_path = write_csv(
        rows,
        output_dir / "data" / "diagnostic_complete_schema.csv",
        [
            "metric_name",
            "criterion_type",
            "required_for_carry_forward",
            "source_output_file",
            "source_milestone",
            "missing_status",
            "description",
        ],
    )
    summary_path = write_csv(
        summary,
        output_dir / "data" / "diagnostic_complete_schema_summary.csv",
        ["criterion_type", "missing_status", "metric_count"],
    )
    return schema_path, summary_path


def main() -> None:
    rows, summary = run_experiment()
    paths = write_outputs(rows, summary)
    print(f"Wrote diagnostic-complete schema: {', '.join(map(str, paths))}")


if __name__ == "__main__":
    main()
