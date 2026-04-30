"""Manifest-level drilldown utilities for blocked v3 protocol decisions."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class V3ProtocolManifestDrilldownRow:
    """Manifest-level diagnostic row for v3 protocol blocked-decision audits."""

    manifest_id: str
    family_name: str
    family_kind: str
    handoff_provenance_type: str
    measurement_protocol_id: str
    profile_invariance_status: str
    eligible: bool
    heldout_violation_rate: float
    train_violation_rate: float
    generalization_gap: float
    target_coverage_fraction: float
    pair_node_coverage_fraction: float
    restart_std: float
    latent_order_disagreement: float
    failed_reasons: list[str]
    warning_reasons: list[str]


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _mean(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.mean(finite)) if finite else float("nan")


def _family_kind_for(family: str) -> str:
    if "failed" in family:
        return "failed_control"
    if "report_only" in family:
        return "report_only"
    return "structured"


def load_v3_protocol_manifest_drilldown_rows(
    output_dir: Path,
) -> list[V3ProtocolManifestDrilldownRow]:
    """Load manifest-level diagnostic rows from M41 v3 protocol outputs."""

    data_dir = output_dir / "data"
    fit_rows = _read_csv(data_dir / "v3_protocol_manifest_family_fit_comparison.csv")
    coverage_rows = _read_csv(
        data_dir / "v3_protocol_manifest_family_coverage_metrics.csv"
    )
    restart_rows = _read_csv(
        data_dir / "v3_protocol_manifest_family_restart_stability.csv"
    )
    latent_rows = _read_csv(
        data_dir / "v3_protocol_manifest_family_latent_order_stability.csv"
    )
    metadata_rows = _read_csv(data_dir / "v3_protocol_manifest_metadata_audit.csv")
    failed_rows = _read_csv(
        data_dir / "v3_protocol_manifest_family_failed_accounting.csv"
    )

    fit_by_manifest: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in fit_rows:
        fit_by_manifest[row.get("manifest_id", "")].append(row)
    coverage_by_manifest = {row.get("manifest_id", ""): row for row in coverage_rows}
    restart_by_manifest = {row.get("manifest_id", ""): row for row in restart_rows}
    latent_by_manifest = {row.get("manifest_id", ""): row for row in latent_rows}
    metadata_by_manifest = {row.get("manifest_id", ""): row for row in metadata_rows}
    failed_by_manifest = {row.get("manifest_id", ""): row for row in failed_rows}

    manifest_ids = sorted(
        set(fit_by_manifest)
        | set(coverage_by_manifest)
        | set(restart_by_manifest)
        | set(latent_by_manifest)
        | set(metadata_by_manifest)
        | set(failed_by_manifest)
    )
    rows: list[V3ProtocolManifestDrilldownRow] = []
    for manifest_id in manifest_ids:
        fits = fit_by_manifest.get(manifest_id, [])
        fit_family = fits[0] if fits else {}
        coverage = coverage_by_manifest.get(manifest_id, {})
        restart = restart_by_manifest.get(manifest_id, {})
        latent = latent_by_manifest.get(manifest_id, {})
        metadata = metadata_by_manifest.get(manifest_id, {})
        failed = failed_by_manifest.get(manifest_id, {})
        family_name = str(
            fit_family.get("family_name")
            or coverage.get("family_name")
            or metadata.get("family_name")
            or failed.get("family_name")
            or ""
        )
        family_kind = str(
            fit_family.get("family_kind")
            or coverage.get("family_kind")
            or failed.get("family_kind")
            or _family_kind_for(family_name)
        )
        train = _mean([_to_float(row.get("train_violation_rate")) for row in fits])
        heldout = _mean([_to_float(row.get("heldout_violation_rate")) for row in fits])
        rows.append(
            V3ProtocolManifestDrilldownRow(
                manifest_id=manifest_id,
                family_name=family_name,
                family_kind=family_kind,
                handoff_provenance_type=str(
                    failed.get("handoff_provenance_type", "")
                ),
                measurement_protocol_id=str(
                    metadata.get("measurement_protocol_id", "")
                ),
                profile_invariance_status=str(
                    metadata.get("profile_invariance_status", "")
                ),
                eligible=_to_float(failed.get("eligible", "1")) == 1.0,
                heldout_violation_rate=heldout,
                train_violation_rate=train,
                generalization_gap=heldout - train,
                target_coverage_fraction=_to_float(
                    coverage.get("target_coverage_fraction")
                ),
                pair_node_coverage_fraction=_to_float(
                    coverage.get("pair_node_coverage_fraction")
                ),
                restart_std=_to_float(restart.get("restart_std")),
                latent_order_disagreement=_to_float(
                    latent.get("latent_order_disagreement")
                ),
                failed_reasons=[
                    item for item in str(failed.get("failed_reasons", "")).split(";")
                    if item
                ],
                warning_reasons=[],
            )
        )
    return rows


def manifest_drilldown_row_to_csv(
    row: V3ProtocolManifestDrilldownRow,
) -> dict[str, float | str]:
    """Convert a manifest drilldown row to a CSV-safe row."""

    data = asdict(row)
    data["eligible"] = float(row.eligible)
    data["failed_reasons"] = ";".join(row.failed_reasons)
    data["warning_reasons"] = ";".join(row.warning_reasons)
    return data  # type: ignore[return-value]


def summarize_manifest_drilldown_by_family(
    rows: list[V3ProtocolManifestDrilldownRow],
) -> list[dict[str, float | str]]:
    """Summarize manifest-level variance by family."""

    grouped: dict[str, list[V3ProtocolManifestDrilldownRow]] = defaultdict(list)
    for row in rows:
        grouped[row.family_name].append(row)
    summary: list[dict[str, float | str]] = []
    for family_name, family_rows in sorted(grouped.items()):
        heldout = [row.heldout_violation_rate for row in family_rows]
        gaps = [row.generalization_gap for row in family_rows]
        finite_heldout = [value for value in heldout if np.isfinite(value)]
        finite_gaps = [value for value in gaps if np.isfinite(value)]
        worst = max(
            family_rows,
            key=lambda item: (
                item.heldout_violation_rate
                if np.isfinite(item.heldout_violation_rate)
                else -1.0
            ),
        )
        summary.append(
            {
                "family_name": family_name,
                "manifest_count": float(len(family_rows)),
                "mean_heldout": _mean(heldout),
                "std_heldout": (
                    float(np.std(finite_heldout))
                    if finite_heldout
                    else float("nan")
                ),
                "max_heldout": (
                    float(max(finite_heldout))
                    if finite_heldout
                    else float("nan")
                ),
                "mean_generalization_gap": _mean(gaps),
                "std_generalization_gap": (
                    float(np.std(finite_gaps)) if finite_gaps else float("nan")
                ),
                "mean_restart_std": _mean(
                    [row.restart_std for row in family_rows]
                ),
                "mean_latent_order_disagreement": _mean(
                    [row.latent_order_disagreement for row in family_rows]
                ),
                "worst_manifest_id": worst.manifest_id,
            }
        )
    return summary
