"""V2 manifest-family specification loading."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path

EXECUTION_STATUSES = {"planned_only", "executed"}


@dataclass(frozen=True)
class V2ManifestFamilySpec:
    """Specification for one diagnostic-complete v2 manifest family."""

    family_name: str
    family_kind: str
    profile_column_policy: str
    target_inclusion_policy: str
    pair_node_coverage_policy: str
    null_taxonomy_policy: str
    restart_stability_required: bool
    latent_order_stability_required: bool
    execution_status: str

    def __post_init__(self) -> None:
        if self.execution_status not in EXECUTION_STATUSES:
            allowed = ", ".join(sorted(EXECUTION_STATUSES))
            raise ValueError(f"execution_status must be one of: {allowed}")


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"1", "1.0", "true", "yes"}


def _spec_from_row(row: dict[str, object]) -> V2ManifestFamilySpec:
    return V2ManifestFamilySpec(
        family_name=str(row.get("family_name", "")),
        family_kind=str(row.get("family_kind", "")),
        profile_column_policy=str(row.get("profile_column_policy", "")),
        target_inclusion_policy=str(row.get("target_inclusion_policy", "")),
        pair_node_coverage_policy=str(row.get("pair_node_coverage_policy", "")),
        null_taxonomy_policy=str(row.get("null_taxonomy_policy", "")),
        restart_stability_required=_as_bool(
            row.get("restart_stability_required", False)
        ),
        latent_order_stability_required=_as_bool(
            row.get("latent_order_stability_required", False)
        ),
        execution_status=str(row.get("execution_status", "planned_only")),
    )


def load_v2_family_specs_from_design_csv(
    path: Path,
) -> list[V2ManifestFamilySpec]:
    """Load v2 family specs from ``new_manifest_family_design_v2.csv``."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return [_spec_from_row(row) for row in csv.DictReader(csv_file)]


def load_v2_family_specs_from_remediation_plan(
    path: Path,
) -> list[V2ManifestFamilySpec]:
    """Load v2 family specs from the Milestone 36 remediation plan JSON."""

    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    specs = payload.get("new_manifest_family_specs", [])
    if not isinstance(specs, list):
        return []
    return [_spec_from_row(dict(row)) for row in specs]


def mark_v2_family_specs_executed(
    specs: list[V2ManifestFamilySpec],
) -> list[V2ManifestFamilySpec]:
    """Return specs marked as executed for the v2 generation run."""

    return [replace(spec, execution_status="executed") for spec in specs]


def v2_family_spec_table(
    specs: list[V2ManifestFamilySpec],
) -> list[dict[str, float | str]]:
    """Return v2 family specs as flat rows."""

    rows: list[dict[str, float | str]] = []
    for spec in specs:
        row = asdict(spec)
        row["restart_stability_required"] = float(spec.restart_stability_required)
        row["latent_order_stability_required"] = float(
            spec.latent_order_stability_required
        )
        rows.append(row)
    return rows
