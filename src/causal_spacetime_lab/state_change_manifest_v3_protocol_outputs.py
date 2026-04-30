"""Load M41 protocol-invariant v3 output bundles for M42 evaluation."""

from __future__ import annotations

import csv
from pathlib import Path

V3_PROTOCOL_OUTPUT_FILES = {
    "metrics": "v3_protocol_cross_family_robustness_metrics.csv",
    "completeness": "v3_protocol_diagnostic_completeness_check.csv",
    "bundle_report": "v3_protocol_diagnostic_complete_bundle_report.csv",
    "generation": "v3_protocol_manifest_generation.csv",
    "metadata_audit": "v3_protocol_manifest_metadata_audit.csv",
    "fit_summary": "v3_protocol_manifest_family_fit_summary.csv",
    "null_taxonomy": "v3_protocol_manifest_family_null_taxonomy.csv",
    "stricter_criteria": "v3_protocol_manifest_family_stricter_criteria.csv",
    "failed_accounting": "v3_protocol_manifest_family_failed_accounting.csv",
    "coverage": "v3_protocol_manifest_family_coverage_metrics.csv",
    "restart_stability": "v3_protocol_manifest_family_restart_stability.csv",
    "latent_order_stability": "v3_protocol_manifest_family_latent_order_stability.csv",
    "no_retuning_audit": "v3_protocol_no_retuning_audit.csv",
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def load_v3_protocol_output_bundle(
    output_dir: Path,
) -> dict[str, list[dict[str, str]]]:
    """Load required M41 v3 protocol CSV outputs from ``outputs/data``."""

    data_dir = output_dir / "data"
    return {
        key: _read_csv(data_dir / filename)
        for key, filename in V3_PROTOCOL_OUTPUT_FILES.items()
    }


def missing_v3_protocol_bundle_inputs(
    bundle: dict[str, list[dict[str, str]]],
) -> list[str]:
    """Return missing or empty v3 protocol bundle inputs."""

    return [key for key in V3_PROTOCOL_OUTPUT_FILES if not bundle.get(key)]


def v3_protocol_family_names_from_bundle(
    bundle: dict[str, list[dict[str, str]]],
) -> list[str]:
    """Return sorted family names present in the v3 protocol bundle."""

    names: set[str] = set()
    for rows in bundle.values():
        for row in rows:
            name = row.get("family_name", "")
            if name and not name.startswith("__"):
                names.add(name)
    return sorted(names)


def v3_protocol_bundle_input_report(
    bundle: dict[str, list[dict[str, str]]],
) -> list[dict[str, float | str]]:
    """Return one present/missing row per required v3 protocol input."""

    rows: list[dict[str, float | str]] = []
    for key, filename in V3_PROTOCOL_OUTPUT_FILES.items():
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

