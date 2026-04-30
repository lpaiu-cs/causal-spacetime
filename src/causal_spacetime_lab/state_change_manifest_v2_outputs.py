"""Load diagnostic-complete v2 output bundles for carry-forward evaluation."""

from __future__ import annotations

import csv
from pathlib import Path

V2_OUTPUT_FILES = {
    "metrics": "v2_cross_family_robustness_metrics.csv",
    "completeness": "v2_diagnostic_completeness_check.csv",
    "bundle_report": "v2_diagnostic_complete_bundle_report.csv",
    "generation": "v2_manifest_generation.csv",
    "fit_summary": "v2_manifest_family_fit_summary.csv",
    "null_taxonomy": "v2_manifest_family_null_taxonomy.csv",
    "stricter_criteria": "v2_manifest_family_stricter_criteria.csv",
    "failed_accounting": "v2_manifest_family_failed_accounting.csv",
    "coverage": "v2_manifest_family_coverage_metrics.csv",
    "restart_stability": "v2_manifest_family_restart_stability.csv",
    "latent_order_stability": "v2_manifest_family_latent_order_stability.csv",
    "no_retuning_audit": "v2_no_retuning_audit.csv",
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def load_v2_output_bundle(output_dir: Path) -> dict[str, list[dict[str, str]]]:
    """Load all required v2 output CSV files from ``outputs/data``."""

    data_dir = output_dir / "data"
    return {
        key: _read_csv(data_dir / filename)
        for key, filename in V2_OUTPUT_FILES.items()
    }


def missing_v2_bundle_inputs(
    bundle: dict[str, list[dict[str, str]]],
) -> list[str]:
    """Return v2 bundle input names that are missing or empty."""

    return [
        key
        for key in V2_OUTPUT_FILES
        if not bundle.get(key)
    ]


def v2_family_names_from_bundle(
    bundle: dict[str, list[dict[str, str]]],
) -> list[str]:
    """Return sorted v2 family names found in bundle rows."""

    names: set[str] = set()
    for rows in bundle.values():
        for row in rows:
            name = row.get("family_name", "")
            if name and not name.startswith("__"):
                names.add(name)
    return sorted(names)


def v2_bundle_input_report(
    bundle: dict[str, list[dict[str, str]]],
) -> list[dict[str, float | str]]:
    """Return one presence row per required v2 bundle input."""

    rows: list[dict[str, float | str]] = []
    for key, filename in V2_OUTPUT_FILES.items():
        row_count = len(bundle.get(key, []))
        rows.append(
            {
                "input_name": key,
                "filename": filename,
                "row_count": float(row_count),
                "present": float(row_count > 0),
            }
        )
    return rows
