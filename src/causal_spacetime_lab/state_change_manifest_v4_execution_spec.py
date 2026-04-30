"""Load v4 protocol execution specs from the M43 preregistration."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path

EXECUTION_STATUSES = {"planned_only", "executed"}


@dataclass(frozen=True)
class V4ProtocolExecutionFamilySpec:
    """Execution spec for one preregistered v4 protocol family."""

    family_name: str
    family_kind: str
    planned_manifest_count: int
    measurement_protocol_family: str
    handoff_provenance_type: str
    profile_resolution_policy: str
    target_inclusion_policy: str
    null_taxonomy_policy: str
    stability_policy: str
    comparison_method: str
    margin_policy: str
    margin_value: float
    linked_v3_root_causes: list[str]
    execution_status: str

    def __post_init__(self) -> None:
        if self.execution_status not in EXECUTION_STATUSES:
            allowed = ", ".join(sorted(EXECUTION_STATUSES))
            raise ValueError(f"execution_status must be one of: {allowed}")
        if self.family_kind == "structured" and self.planned_manifest_count < 3:
            raise ValueError("structured v4 specs require at least 3 manifests")
        if self.family_kind != "structured" and self.planned_manifest_count < 1:
            raise ValueError("v4 control specs require at least 1 manifest")


def load_v4_protocol_preregistration(path: Path) -> dict[str, object]:
    """Load the M43 v4 preregistration JSON if present."""

    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return dict(payload) if isinstance(payload, dict) else {}


def load_v4_protocol_execution_specs(
    path: Path,
) -> list[V4ProtocolExecutionFamilySpec]:
    """Read v4 execution specs from the M43 preregistration."""

    payload = load_v4_protocol_preregistration(path)
    if not payload:
        return []
    criteria_source = str(payload.get("fixed_criteria_source", ""))
    if "M34" not in criteria_source and "Milestone 34" not in criteria_source:
        raise ValueError("v4 preregistration must retain fixed M34 criteria")
    families = payload.get("planned_families", [])
    if not isinstance(families, list):
        return []
    specs: list[V4ProtocolExecutionFamilySpec] = []
    for item in families:
        if not isinstance(item, dict):
            continue
        execution_status = str(item.get("execution_status", "planned_only"))
        if execution_status != "planned_only":
            raise ValueError("v4 specs must be planned_only before M44 execution")
        specs.append(_spec_from_row(item))
    return specs


def mark_v4_protocol_execution_specs_executed(
    specs: list[V4ProtocolExecutionFamilySpec],
) -> list[V4ProtocolExecutionFamilySpec]:
    """Return copies of v4 execution specs marked executed."""

    return [replace(spec, execution_status="executed") for spec in specs]


def v4_protocol_execution_spec_table(
    specs: list[V4ProtocolExecutionFamilySpec],
) -> list[dict[str, float | str]]:
    """Return v4 execution specs as CSV-safe rows."""

    rows: list[dict[str, float | str]] = []
    for spec in specs:
        rows.append(
            {
                "family_name": spec.family_name,
                "family_kind": spec.family_kind,
                "planned_manifest_count": float(spec.planned_manifest_count),
                "measurement_protocol_family": spec.measurement_protocol_family,
                "handoff_provenance_type": spec.handoff_provenance_type,
                "profile_resolution_policy": spec.profile_resolution_policy,
                "target_inclusion_policy": spec.target_inclusion_policy,
                "null_taxonomy_policy": spec.null_taxonomy_policy,
                "stability_policy": spec.stability_policy,
                "comparison_method": spec.comparison_method,
                "margin_policy": spec.margin_policy,
                "margin_value": float(spec.margin_value),
                "linked_v3_root_causes": ";".join(spec.linked_v3_root_causes),
                "execution_status": spec.execution_status,
            }
        )
    return rows


def _spec_from_row(row: dict[str, object]) -> V4ProtocolExecutionFamilySpec:
    linked = row.get("linked_v3_root_causes", [])
    if isinstance(linked, str):
        linked_causes = [item for item in linked.split(";") if item]
    elif isinstance(linked, list):
        linked_causes = [str(item) for item in linked]
    else:
        linked_causes = []
    return V4ProtocolExecutionFamilySpec(
        family_name=str(row.get("family_name", "")),
        family_kind=str(row.get("family_kind", "")),
        planned_manifest_count=int(float(row.get("planned_manifest_count", 1))),
        measurement_protocol_family=str(row.get("measurement_protocol_family", "")),
        handoff_provenance_type=str(row.get("handoff_provenance_type", "")),
        profile_resolution_policy=str(row.get("profile_resolution_policy", "")),
        target_inclusion_policy=str(row.get("target_inclusion_policy", "")),
        null_taxonomy_policy=str(row.get("null_taxonomy_policy", "")),
        stability_policy=str(row.get("stability_policy", "")),
        comparison_method=str(row.get("comparison_method", "")),
        margin_policy=str(row.get("margin_policy", "fixed_margin")),
        margin_value=float(row.get("margin_value", 0.05)),
        linked_v3_root_causes=linked_causes,
        execution_status=str(row.get("execution_status", "planned_only")),
    )
