"""Shared helpers for carry-forward failure decomposition experiments."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from manifest_family_experiment_helpers import save_bar_figure, write_csv

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
    failure_record_to_dict,
)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows, returning an empty list when the file is absent."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def data_path(output_dir: Path, filename: str) -> Path:
    """Return a file path under the output data directory."""

    return output_dir / "data" / filename


def write_failure_records(
    records: list[CriterionFailureRecord],
    path: Path,
) -> Path:
    """Write criterion failure records to CSV."""

    rows = [failure_record_to_dict(record) for record in records]
    return write_csv(
        rows,
        path,
        [
            "family_name",
            "family_kind",
            "criterion_name",
            "criterion_type",
            "observed_value",
            "threshold_value",
            "status",
            "root_cause_category",
            "explanation",
        ],
    )


def failure_record_from_row(row: dict[str, str]) -> CriterionFailureRecord:
    """Reconstruct a criterion failure record from a CSV row."""

    return CriterionFailureRecord(
        family_name=row.get("family_name", ""),
        family_kind=row.get("family_kind", ""),
        criterion_name=row.get("criterion_name", ""),
        criterion_type=row.get("criterion_type", ""),
        observed_value=float(row.get("observed_value", "nan")),
        threshold_value=float(row.get("threshold_value", "nan")),
        status=row.get("status", ""),
        root_cause_category=row.get("root_cause_category", ""),
        explanation=row.get("explanation", ""),
    )


def save_count_figure(
    rows: list[dict[str, float | str]],
    key: str,
    path: Path,
    *,
    ylabel: str,
) -> Path:
    """Save a count bar chart for a categorical row field."""

    counts: Counter[str] = Counter(str(row.get(key, "")) for row in rows)
    labels = sorted(label for label in counts if label)
    return save_bar_figure(
        labels,
        [float(counts[label]) for label in labels],
        path,
        ylabel=ylabel,
    )
