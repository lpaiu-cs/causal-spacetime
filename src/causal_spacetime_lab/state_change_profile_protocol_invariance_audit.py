"""Audit response-profile protocol invariance metadata."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    measurement_protocol_id,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    ResponseProfileMetadata,
    profile_metadata_from_protocols,
)


@dataclass(frozen=True)
class ProfileProtocolInvarianceAuditRow:
    """One protocol-invariance audit row."""

    family_id: str
    manifest_id: str
    measurement_protocol_ids_used: list[str]
    reference_set_id: str
    number_of_protocols_mixed: int
    number_of_references: int
    profile_invariance_status: str
    admissible_for_pairwise_dissimilarity: bool
    reason_if_blocked: str
    source: str


def audit_profile_metadata(
    metadata: ResponseProfileMetadata,
    *,
    manifest_id: str = "",
    source: str = "",
) -> ProfileProtocolInvarianceAuditRow:
    """Audit explicit response-profile metadata."""

    if metadata.profile_invariance_status == "protocol_invariant":
        protocol_ids = [metadata.measurement_protocol_id]
        mixed_count = 1
    elif metadata.profile_invariance_status == "protocol_mixed":
        protocol_ids = [
            item
            for item in metadata.measurement_protocol_hash.split(";")
            if item
        ]
        mixed_count = max(2, len(protocol_ids))
    else:
        protocol_ids = []
        mixed_count = 0
    return ProfileProtocolInvarianceAuditRow(
        family_id=metadata.profile_family_id,
        manifest_id=manifest_id,
        measurement_protocol_ids_used=protocol_ids,
        reference_set_id=metadata.reference_set_id,
        number_of_protocols_mixed=mixed_count,
        number_of_references=len(metadata.reference_chain_ids),
        profile_invariance_status=metadata.profile_invariance_status,
        admissible_for_pairwise_dissimilarity=(
            metadata.admissible_for_pairwise_dissimilarity
        ),
        reason_if_blocked=metadata.reason_if_blocked,
        source=source,
    )


def audit_manifest_protocol_metadata(
    manifest_json: dict[str, object],
) -> ProfileProtocolInvarianceAuditRow:
    """Audit measurement-protocol metadata embedded in a manifest JSON."""

    manifest_id = str(manifest_json.get("manifest_id", ""))
    family_id = str(
        manifest_json.get("profile_family_id")
        or manifest_json.get("profile_label")
        or manifest_id
    )
    metadata_json = manifest_json.get("response_profile_metadata")
    if not isinstance(metadata_json, dict):
        return ProfileProtocolInvarianceAuditRow(
            family_id=family_id,
            manifest_id=manifest_id,
            measurement_protocol_ids_used=[],
            reference_set_id="",
            number_of_protocols_mixed=0,
            number_of_references=0,
            profile_invariance_status="underspecified",
            admissible_for_pairwise_dissimilarity=False,
            reason_if_blocked="missing measurement protocol metadata",
            source="manifest",
        )
    protocols_json = metadata_json.get("measurement_protocols")
    reference_ids = _string_list(metadata_json.get("reference_chain_ids"))
    reference_set_id = str(metadata_json.get("reference_set_id", ""))
    if not isinstance(protocols_json, list):
        protocol_json = metadata_json.get("measurement_protocol")
        protocols_json = [protocol_json] if isinstance(protocol_json, dict) else []
    protocols = [
        _protocol_from_json(item)
        for item in protocols_json
        if isinstance(item, dict)
    ]
    metadata = profile_metadata_from_protocols(
        profile_family_id=family_id,
        protocols=protocols,
        reference_chain_ids=reference_ids,
        reference_set_id=reference_set_id,
        exploratory_mixed_context=bool(
            metadata_json.get("exploratory_mixed_context", False)
        ),
    )
    return audit_profile_metadata(metadata, manifest_id=manifest_id, source="manifest")


def audit_manifest_directory_protocol_invariance(
    manifest_dir: Path,
    *,
    source: str,
) -> list[ProfileProtocolInvarianceAuditRow]:
    """Audit all manifest JSON files in a directory."""

    rows: list[ProfileProtocolInvarianceAuditRow] = []
    if not manifest_dir.exists():
        return rows
    for path in sorted(manifest_dir.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        row = audit_manifest_protocol_metadata(payload)
        rows.append(
            ProfileProtocolInvarianceAuditRow(
                family_id=row.family_id,
                manifest_id=row.manifest_id,
                measurement_protocol_ids_used=row.measurement_protocol_ids_used,
                reference_set_id=row.reference_set_id,
                number_of_protocols_mixed=row.number_of_protocols_mixed,
                number_of_references=row.number_of_references,
                profile_invariance_status=row.profile_invariance_status,
                admissible_for_pairwise_dissimilarity=(
                    row.admissible_for_pairwise_dissimilarity
                ),
                reason_if_blocked=row.reason_if_blocked,
                source=source,
            )
        )
    return rows


def audit_v3_design_protocol_invariance(
    design_rows: list[dict[str, str | float]],
) -> list[ProfileProtocolInvarianceAuditRow]:
    """Audit planned v3 design rows before protocol metadata patching."""

    rows: list[ProfileProtocolInvarianceAuditRow] = []
    for design in design_rows:
        family_id = str(design.get("family_name", ""))
        rows.append(
            ProfileProtocolInvarianceAuditRow(
                family_id=family_id,
                manifest_id="",
                measurement_protocol_ids_used=[],
                reference_set_id="",
                number_of_protocols_mixed=0,
                number_of_references=0,
                profile_invariance_status="underspecified",
                admissible_for_pairwise_dissimilarity=False,
                reason_if_blocked="missing measurement protocol metadata",
                source="v3_design",
            )
        )
    return rows


def audit_row_to_dict(
    row: ProfileProtocolInvarianceAuditRow,
) -> dict[str, float | str]:
    """Convert an audit row to CSV data."""

    return {
        "family_id": row.family_id,
        "manifest_id": row.manifest_id,
        "measurement_protocol_ids_used": ";".join(row.measurement_protocol_ids_used),
        "reference_set_id": row.reference_set_id,
        "number_of_protocols_mixed": float(row.number_of_protocols_mixed),
        "number_of_references": float(row.number_of_references),
        "profile_invariance_status": row.profile_invariance_status,
        "admissible_for_pairwise_dissimilarity": float(
            row.admissible_for_pairwise_dissimilarity
        ),
        "reason_if_blocked": row.reason_if_blocked,
        "source": row.source,
    }


def _string_list(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(item) for item in value]
    return []


def _protocol_from_json(payload: dict[str, Any]) -> MeasurementProtocolSpec:
    return MeasurementProtocolSpec(
        echo_rule=str(payload["echo_rule"]),
        emission_rule=str(payload["emission_rule"]),
        gate_rule=str(payload["gate_rule"]),
        reference_subsampling_rule=str(payload["reference_subsampling_rule"]),
        spectrum_type=str(payload["spectrum_type"]),
        normalization_rule=str(payload["normalization_rule"]),
        missing_policy=str(payload["missing_policy"]),
        tie_policy=str(payload["tie_policy"]),
        margin_policy=str(payload["margin_policy"]),
        protocol_label=str(payload.get("protocol_label", "")),
    )


def protocol_ids_for_specs(protocols: list[MeasurementProtocolSpec]) -> list[str]:
    """Return stable identifiers for measurement-protocol specs."""

    return [measurement_protocol_id(protocol) for protocol in protocols]
