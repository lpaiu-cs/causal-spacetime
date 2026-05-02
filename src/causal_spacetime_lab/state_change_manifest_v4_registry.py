"""V4 protocol carry-forward registry utilities."""

from __future__ import annotations

import re
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    CarryForwardRegistry,
    build_carry_forward_registry,
    write_carry_forward_registry,
)
from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_datasets
from causal_spacetime_lab.state_change_manifest_family import (
    ManifestFamilyAssignment,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)


def _family_from_label(label: str) -> str:
    return re.sub(r"_m\d+$", "", label)


def _v4_assignments_from_manifest_dir(
    manifest_dir: Path,
) -> list[ManifestFamilyAssignment]:
    assignments: list[ManifestFamilyAssignment] = []
    for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True):
        family_name = _family_from_label(
            str(dataset.manifest_json.get("profile_label", dataset.manifest_id))
        )
        if "failed" in family_name:
            family_kind = "failed_control"
        elif "report_only" in family_name:
            family_kind = "report_only"
        else:
            family_kind = "structured"
        assignments.append(
            ManifestFamilyAssignment(
                manifest_id=dataset.manifest_id,
                family_name=family_name,
                family_kind=family_kind,
                eligible=dataset.eligible,
                failed_reasons=list(dataset.failed_reasons),
            )
        )
    return assignments


def build_v4_protocol_carry_forward_registry(
    decisions: list[FamilyRobustnessDecision],
    manifest_dir: Path,
    criteria_name: str = "fixed_m34_v4_protocol",
) -> CarryForwardRegistry:
    """Build a v4 protocol carry-forward registry from family decisions."""

    registry = build_carry_forward_registry(
        decisions,
        manifest_assignments=_v4_assignments_from_manifest_dir(manifest_dir),
        criteria_name=criteria_name,
    )
    registry.created_by_milestone = "Milestone 45"
    return registry


def write_v4_protocol_carry_forward_registry(
    registry: CarryForwardRegistry,
    output_dir: Path,
) -> Path:
    """Write the v4 protocol carry-forward registry JSON."""

    return write_carry_forward_registry(
        registry,
        output_dir
        / "carry_forward_v4_protocol"
        / "carry_forward_registry_v4_protocol.json",
    )


def v4_protocol_registry_summary_rows(
    registry: CarryForwardRegistry,
) -> list[dict[str, float | str]]:
    """Return CSV-safe summary rows for a v4 protocol registry."""

    rows: list[dict[str, float | str]] = []
    for record in registry.records:
        rows.append(
            {
                "registry_id": registry.registry_id,
                "criteria_name": registry.criteria_name,
                "family_name": record.family_name,
                "family_kind": record.family_kind,
                "decision": record.decision,
                "manifest_count": float(len(record.manifest_ids)),
                "eligible_manifest_count": float(record.eligible_manifest_count),
                "ineligible_manifest_count": float(record.ineligible_manifest_count),
                "failed_reasons": ";".join(record.failed_reasons),
                "warning_reasons": ";".join(record.warning_reasons),
                "allowed_future_use": record.allowed_future_use,
            }
        )
    return rows
