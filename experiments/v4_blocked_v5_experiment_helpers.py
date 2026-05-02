"""Shared helpers for Milestone 46 v4 blocked/v5 preregistration experiments."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def data_path(output_dir: Path, filename: str) -> Path:
    """Return an M46 output data path."""

    return output_dir / "data" / filename


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows, returning empty rows when absent."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def to_float_rows(rows: list[dict[str, str]]) -> list[dict[str, float | str]]:
    """Convert numeric-looking CSV values to floats."""

    converted: list[dict[str, float | str]] = []
    for row in rows:
        new_row: dict[str, float | str] = {}
        for key, value in row.items():
            try:
                new_row[key] = float(value)
            except (TypeError, ValueError):
                new_row[key] = value
        converted.append(new_row)
    return converted


def missing_input_row(filename: str) -> dict[str, float | str]:
    """Return an explicit missing-input row."""

    return {
        "row_type": "missing_input",
        "input_file": filename,
        "status": "missing_input",
        "count": 0.0,
    }


def read_json(path: Path) -> dict[str, object]:
    """Read JSON if present."""

    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))
