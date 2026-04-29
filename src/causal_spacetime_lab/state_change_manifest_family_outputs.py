"""Load Milestone 33 manifest-family output bundles."""

from __future__ import annotations

import csv
from pathlib import Path

FAMILY_OUTPUT_FILES = {
    "fit_summary": "manifest_family_fit_summary.csv",
    "null_taxonomy": "manifest_family_null_taxonomy.csv",
    "stricter_criteria": "manifest_family_stricter_criteria.csv",
    "failed_accounting": "manifest_family_failed_manifest_accounting.csv",
    "no_retuning_audit": "manifest_family_no_retuning_audit.csv",
    "report_card": "manifest_family_report_card.csv",
}


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    """Load CSV rows, returning an empty list when the file is absent."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def load_family_output_bundle(output_dir: Path) -> dict[str, list[dict[str, str]]]:
    """Load expected Milestone 33 family-output CSVs."""

    data_dir = output_dir / "data"
    return {
        key: load_csv_rows(data_dir / filename)
        for key, filename in FAMILY_OUTPUT_FILES.items()
    }


def missing_family_output_inputs(
    bundle: dict[str, list[dict[str, str]]],
) -> list[str]:
    """Return bundle keys with no loaded rows."""

    return [key for key in FAMILY_OUTPUT_FILES if not bundle.get(key)]


def family_names_from_bundle(bundle: dict[str, list[dict[str, str]]]) -> list[str]:
    """Return sorted family names found in loaded family outputs."""

    names: set[str] = set()
    for key in ["fit_summary", "stricter_criteria", "report_card"]:
        for row in bundle.get(key, []):
            family_name = row.get("family_name", "")
            if family_name:
                names.add(family_name)
    return sorted(names)
