"""Planned-only v3 protocol-invariant family patches."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    default_earliest_full_protocol,
    default_earliest_retained_protocol,
    default_gated_full_protocol,
    default_immediate_edge_protocol,
    measurement_protocol_hash,
    measurement_protocol_id,
)

EXECUTION_STATUSES = {"planned_only"}


@dataclass(frozen=True)
class V3ProtocolInvariantFamilyPatch:
    """Planned-only protocol-invariant v3 family patch."""

    original_family_name: str
    patched_family_name: str
    family_kind: str
    measurement_protocol_id: str
    measurement_protocol_hash: str
    echo_rule: str
    emission_rule: str
    gate_rule: str
    reference_subsampling_rule: str
    spectrum_type: str
    normalization_rule: str
    missing_policy: str
    tie_policy: str
    margin_policy: str
    planned_manifest_count: int
    profile_family_semantics: str
    admissible_for_pairwise_dissimilarity: bool
    execution_status: str

    def __post_init__(self) -> None:
        if self.execution_status not in EXECUTION_STATUSES:
            allowed = ", ".join(sorted(EXECUTION_STATUSES))
            raise ValueError(f"execution_status must be one of: {allowed}")


def default_v3_protocol_invariant_family_patches() -> list[
    V3ProtocolInvariantFamilyPatch
]:
    """Return planned-only protocol-invariant v3 family patches."""

    return [
        _patch(
            original_family_name="rank_gap_multi_manifest_v3",
            patched_family_name="rank_gap_earliest_full_reference_v3",
            family_kind="structured",
            protocol=default_earliest_full_protocol(),
            planned_manifest_count=5,
            admissible=True,
        ),
        _patch(
            original_family_name="rank_gap_more_protocol_columns_v3",
            patched_family_name="rank_gap_earliest_retained_reference_v3",
            family_kind="structured",
            protocol=default_earliest_retained_protocol(),
            planned_manifest_count=5,
            admissible=True,
        ),
        _patch(
            original_family_name="rank_gap_rank_resolution_enriched_v3",
            patched_family_name="rank_gap_gated_full_reference_v3",
            family_kind="structured",
            protocol=default_gated_full_protocol(),
            planned_manifest_count=5,
            admissible=True,
        ),
        _patch(
            original_family_name="rank_gap_coverage_enriched_v3",
            patched_family_name="rank_gap_coverage_enriched_full_reference_v3",
            family_kind="structured",
            protocol=default_earliest_full_protocol(),
            planned_manifest_count=5,
            admissible=True,
            semantics=(
                "protocol-invariant coverage-enriched response profile family"
            ),
        ),
        _patch(
            original_family_name="combined_diagnostic_complete_v3",
            patched_family_name="combined_earliest_full_reference_v3",
            family_kind="structured",
            protocol=default_earliest_full_protocol(
                missing_policy="penalize_mismatch"
            ),
            planned_manifest_count=5,
            admissible=True,
        ),
        _patch(
            original_family_name="new_report_only_protocol_family",
            patched_family_name="rank_gap_immediate_edge_report_only_v3",
            family_kind="report_only",
            protocol=default_immediate_edge_protocol(),
            planned_manifest_count=5,
            admissible=False,
            semantics="report-only immediate-edge protocol family",
        ),
        _patch(
            original_family_name="failed_controls_v3",
            patched_family_name="failed_controls_v3",
            family_kind="failed_control",
            protocol=default_earliest_full_protocol(),
            planned_manifest_count=3,
            admissible=False,
            semantics="report-only failed controls with explicit metadata",
        ),
    ]


def v3_protocol_patch_table(
    patches: list[V3ProtocolInvariantFamilyPatch],
) -> list[dict[str, float | str]]:
    """Convert v3 protocol patches to rows."""

    return [asdict(patch) for patch in patches]


def write_v3_protocol_patch_json(
    patches: list[V3ProtocolInvariantFamilyPatch],
    output_path: Path,
) -> Path:
    """Write planned-only v3 protocol patch JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = v3_protocol_patch_table(patches)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return output_path


def read_v3_protocol_patch_json(path: Path) -> list[dict[str, object]]:
    """Read planned-only v3 protocol patch JSON rows."""

    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        return []
    return [dict(item) for item in payload if isinstance(item, dict)]


def _patch(
    *,
    original_family_name: str,
    patched_family_name: str,
    family_kind: str,
    protocol: MeasurementProtocolSpec,
    planned_manifest_count: int,
    admissible: bool,
    semantics: str = "protocol-invariant multi-reference response profile",
) -> V3ProtocolInvariantFamilyPatch:
    return V3ProtocolInvariantFamilyPatch(
        original_family_name=original_family_name,
        patched_family_name=patched_family_name,
        family_kind=family_kind,
        measurement_protocol_id=measurement_protocol_id(protocol),
        measurement_protocol_hash=measurement_protocol_hash(protocol),
        echo_rule=protocol.echo_rule,
        emission_rule=protocol.emission_rule,
        gate_rule=protocol.gate_rule,
        reference_subsampling_rule=protocol.reference_subsampling_rule,
        spectrum_type=protocol.spectrum_type,
        normalization_rule=protocol.normalization_rule,
        missing_policy=protocol.missing_policy,
        tie_policy=protocol.tie_policy,
        margin_policy=protocol.margin_policy,
        planned_manifest_count=planned_manifest_count,
        profile_family_semantics=semantics,
        admissible_for_pairwise_dissimilarity=admissible,
        execution_status="planned_only",
    )
