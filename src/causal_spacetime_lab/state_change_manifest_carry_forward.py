"""Carry-forward registry for manifest-family robustness decisions."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_manifest_family import (
    ManifestFamilyAssignment,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)


@dataclass
class CarryForwardFamilyRecord:
    """One family decision exported for future stress-test planning."""

    family_name: str
    family_kind: str
    decision: str
    manifest_ids: list[str]
    eligible_manifest_count: int
    ineligible_manifest_count: int
    failed_reasons: list[str]
    warning_reasons: list[str]
    key_metrics: dict[str, float]
    allowed_future_use: str
    forbidden_interpretations: list[str]


@dataclass
class CarryForwardRegistry:
    """Frozen family-level carry-forward registry."""

    registry_id: str
    created_by_milestone: str
    criteria_name: str
    records: list[CarryForwardFamilyRecord]
    stop_rules: list[str]
    forbidden_interpretations: list[str]


def forbidden_carry_forward_interpretations() -> list[str]:
    """Return forbidden interpretations for carry-forward decisions."""

    return [
        "carry-forward proves geometry",
        "family-level fit recovers space",
        "latent coordinates are physical coordinates",
        "selected dimension is physical dimension",
        "null separation proves metric structure",
    ]


def _allowed_future_use(decision: str) -> str:
    if decision == "carry_forward":
        return "eligible_for_future_stress_tests"
    if decision == "provisional":
        return "limited_stress_tests_with_caveats"
    return "report_only"


def _assignment_group(
    assignments: list[ManifestFamilyAssignment] | None,
    family_name: str,
) -> list[ManifestFamilyAssignment]:
    if assignments is None:
        return []
    return [
        assignment
        for assignment in assignments
        if assignment.family_name == family_name
    ]


def _clean_float(value: float) -> float | str:
    return float(value) if np.isfinite(value) else "nan"


def build_carry_forward_registry(
    decisions: list[FamilyRobustnessDecision],
    manifest_assignments: list[ManifestFamilyAssignment] | None = None,
    criteria_name: str = "default",
) -> CarryForwardRegistry:
    """Build a frozen carry-forward registry from robustness decisions."""

    forbidden = forbidden_carry_forward_interpretations()
    records: list[CarryForwardFamilyRecord] = []
    for decision in decisions:
        assignments = _assignment_group(manifest_assignments, decision.family_name)
        assignment_reasons = sorted(
            {
                reason
                for assignment in assignments
                for reason in assignment.failed_reasons
            }
        )
        failed_reasons = sorted(set(decision.failed_reasons + assignment_reasons))
        records.append(
            CarryForwardFamilyRecord(
                family_name=decision.family_name,
                family_kind=decision.family_kind,
                decision=decision.decision,
                manifest_ids=[assignment.manifest_id for assignment in assignments],
                eligible_manifest_count=sum(
                    assignment.eligible for assignment in assignments
                ),
                ineligible_manifest_count=sum(
                    not assignment.eligible for assignment in assignments
                ),
                failed_reasons=failed_reasons,
                warning_reasons=list(decision.warning_reasons),
                key_metrics=dict(decision.key_metrics),
                allowed_future_use=_allowed_future_use(decision.decision),
                forbidden_interpretations=forbidden,
            )
        )
    registry = CarryForwardRegistry(
        registry_id="",
        created_by_milestone="Milestone 34",
        criteria_name=criteria_name,
        records=records,
        stop_rules=[
            "if_no_carry_forward_family_then_future_stress_tests_stop_or_report_only",
            "provisional_families_require_explicit_caveats",
            "blocked_families_are_report_only",
        ],
        forbidden_interpretations=forbidden,
    )
    jsonable = registry_to_jsonable(registry)
    registry.registry_id = registry_digest(jsonable)
    return registry


def registry_to_jsonable(registry: CarryForwardRegistry) -> dict[str, object]:
    """Convert registry dataclasses and float values to JSON-safe objects."""

    def convert(value: object) -> object:
        if isinstance(value, float):
            return _clean_float(value)
        if isinstance(value, dict):
            return {str(key): convert(item) for key, item in value.items()}
        if isinstance(value, list):
            return [convert(item) for item in value]
        return value

    return convert(asdict(registry))  # type: ignore[return-value]


def registry_digest(registry_jsonable: dict[str, object]) -> str:
    """Return SHA256 digest for canonical registry JSON."""

    payload = dict(registry_jsonable)
    payload["registry_id"] = ""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_carry_forward_registry(
    registry: CarryForwardRegistry,
    output_path: Path,
) -> Path:
    """Write carry-forward registry JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = registry_to_jsonable(registry)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return output_path


def registry_from_jsonable(payload: dict[str, object]) -> CarryForwardRegistry:
    """Reconstruct a carry-forward registry from JSON data."""

    records = [
        CarryForwardFamilyRecord(
            family_name=str(record["family_name"]),
            family_kind=str(record["family_kind"]),
            decision=str(record["decision"]),
            manifest_ids=list(record.get("manifest_ids", [])),
            eligible_manifest_count=int(record.get("eligible_manifest_count", 0)),
            ineligible_manifest_count=int(record.get("ineligible_manifest_count", 0)),
            failed_reasons=list(record.get("failed_reasons", [])),
            warning_reasons=list(record.get("warning_reasons", [])),
            key_metrics={
                str(key): float(value)
                for key, value in dict(record.get("key_metrics", {})).items()
                if str(value) != "nan"
            },
            allowed_future_use=str(record.get("allowed_future_use", "")),
            forbidden_interpretations=list(
                record.get("forbidden_interpretations", [])
            ),
        )
        for record in payload.get("records", [])
        if isinstance(record, dict)
    ]
    return CarryForwardRegistry(
        registry_id=str(payload.get("registry_id", "")),
        created_by_milestone=str(payload.get("created_by_milestone", "")),
        criteria_name=str(payload.get("criteria_name", "")),
        records=records,
        stop_rules=list(payload.get("stop_rules", [])),
        forbidden_interpretations=list(payload.get("forbidden_interpretations", [])),
    )
