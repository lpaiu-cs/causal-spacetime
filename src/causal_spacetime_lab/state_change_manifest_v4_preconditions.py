"""Provenance-aware preconditions for v4 protocol carry-forward evaluation."""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_handoff_provenance import (
    HandoffProvenanceMetadata,
    validate_handoff_provenance,
)


@dataclass(frozen=True)
class V4ProtocolPreconditionReport:
    """Family-level metadata/provenance precondition report."""

    family_name: str
    family_kind: str
    manifest_count: int
    diagnostic_complete: bool
    all_manifests_have_measurement_protocol: bool
    all_manifests_have_profile_metadata: bool
    all_manifests_have_handoff_provenance: bool
    all_structured_protocol_invariant: bool
    all_structured_parameter_complete: bool
    all_structured_admissible_for_pairwise_dissimilarity: bool
    all_structured_valid_provenance: bool
    report_only_controls_ineligible: bool
    failed_controls_ineligible: bool
    preconditions_passed: bool
    failed_preconditions: list[str]
    warning_preconditions: list[str]


def evaluate_v4_protocol_preconditions(
    manifest_dir: Path,
    bundle: dict[str, list[dict[str, str]]],
) -> list[V4ProtocolPreconditionReport]:
    """Evaluate v4 admissibility metadata preconditions by family."""

    completeness = _diagnostic_complete_by_family(bundle)
    family_kinds = _family_kind_by_family(bundle)
    manifests = _manifest_jsons_by_family(manifest_dir)
    families = sorted(set(family_kinds) | set(completeness) | set(manifests))
    reports: list[V4ProtocolPreconditionReport] = []
    for family_name in families:
        family_kind = family_kinds.get(family_name, _infer_family_kind(family_name))
        family_manifests = manifests.get(family_name, [])
        diagnostic_complete = completeness.get(family_name, False)
        failed: list[str] = []
        warnings: list[str] = []
        if not family_manifests:
            failed.append("missing_manifest_files")
        if not diagnostic_complete:
            failed.append("diagnostic_incomplete")

        has_measurement = _all_have(family_manifests, "measurement_protocol")
        has_profile = _all_have(family_manifests, "profile_metadata")
        has_provenance = _all_have(family_manifests, "handoff_provenance")
        if not has_measurement:
            failed.append("missing_measurement_protocol")
        if not has_profile:
            failed.append("missing_profile_metadata")
        if not has_provenance:
            failed.append("missing_handoff_provenance")

        structured = family_kind == "structured"
        protocol_invariant = True
        parameter_complete = True
        admissible = True
        valid_provenance = True
        report_only_ineligible = True
        failed_control_ineligible = True
        if structured:
            protocol_invariant = all(
                item.get("profile_invariance_status") == "protocol_invariant"
                for item in family_manifests
            )
            parameter_complete = all(
                _truthy(item.get("profile_metadata", {}).get("measurement_protocol"))
                for item in family_manifests
            )
            admissible = all(
                bool(item.get("admissible_for_pairwise_dissimilarity"))
                for item in family_manifests
            )
            valid_provenance = all(
                _manifest_provenance_valid(item) for item in family_manifests
            )
            if not protocol_invariant:
                failed.append("protocol_mixed")
            if not parameter_complete:
                failed.append("parameter_incomplete")
            if not admissible:
                failed.append("structured_not_admissible")
            if not valid_provenance:
                failed.append("invalid_provenance")
        else:
            controls_ineligible = all(
                not bool(item.get("handoff_decision", {}).get("eligible"))
                for item in family_manifests
            )
            if family_kind == "failed_control":
                failed_control_ineligible = controls_ineligible
                if not failed_control_ineligible:
                    failed.append("failed_control_marked_eligible")
                warnings.append("failed_controls_remain_ineligible")
            else:
                report_only_ineligible = controls_ineligible
                if not report_only_ineligible:
                    failed.append("report_only_marked_eligible")
                warnings.append("report_only_controls_remain_ineligible")
        reports.append(
            V4ProtocolPreconditionReport(
                family_name=family_name,
                family_kind=family_kind,
                manifest_count=len(family_manifests),
                diagnostic_complete=diagnostic_complete,
                all_manifests_have_measurement_protocol=has_measurement,
                all_manifests_have_profile_metadata=has_profile,
                all_manifests_have_handoff_provenance=has_provenance,
                all_structured_protocol_invariant=protocol_invariant,
                all_structured_parameter_complete=parameter_complete,
                all_structured_admissible_for_pairwise_dissimilarity=admissible,
                all_structured_valid_provenance=valid_provenance,
                report_only_controls_ineligible=report_only_ineligible,
                failed_controls_ineligible=failed_control_ineligible,
                preconditions_passed=not failed,
                failed_preconditions=sorted(set(failed)),
                warning_preconditions=sorted(set(warnings)),
            )
        )
    return reports


def v4_protocol_precondition_report_to_row(
    report: V4ProtocolPreconditionReport,
) -> dict[str, float | str]:
    """Convert one v4 precondition report to a CSV row."""

    row = asdict(report)
    for key, value in list(row.items()):
        if isinstance(value, bool):
            row[key] = float(value)
        elif isinstance(value, int):
            row[key] = float(value)
        elif isinstance(value, list):
            row[key] = ";".join(str(item) for item in value)
    return row  # type: ignore[return-value]


def v4_protocol_precondition_summary(
    rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Summarize v4 precondition pass/fail counts by family kind."""

    counts: Counter[tuple[str, str]] = Counter()
    for row in rows:
        status = (
            "passed"
            if float(row.get("preconditions_passed", 0.0)) == 1.0
            else "failed"
        )
        counts[(str(row.get("family_kind", "")), status)] += 1
    return [
        {"family_kind": kind, "status": status, "count": float(count)}
        for (kind, status), count in sorted(counts.items())
    ]


def _diagnostic_complete_by_family(
    bundle: dict[str, list[dict[str, str]]],
) -> dict[str, bool]:
    status: dict[str, bool] = {}
    for row in bundle.get("completeness", []):
        family = row.get("family_name", "")
        if family:
            status[family] = float(row.get("diagnostic_complete", "0") or 0) == 1.0
    return status


def _family_kind_by_family(bundle: dict[str, list[dict[str, str]]]) -> dict[str, str]:
    kinds: dict[str, str] = {}
    for key in ("metrics", "generation", "completeness"):
        for row in bundle.get(key, []):
            family = row.get("family_name", "")
            kind = row.get("family_kind", "")
            if family and kind:
                kinds[family] = kind
    return kinds


def _manifest_jsons_by_family(manifest_dir: Path) -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for path in sorted(manifest_dir.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        family = _family_from_profile_label(str(payload.get("profile_label", "")))
        if family:
            grouped.setdefault(family, []).append(payload)
    return grouped


def _family_from_profile_label(label: str) -> str:
    return re.sub(r"_m\d+$", "", label)


def _infer_family_kind(family_name: str) -> str:
    if "failed" in family_name:
        return "failed_control"
    if "report_only" in family_name:
        return "report_only"
    return "structured"


def _all_have(manifests: list[dict[str, object]], key: str) -> bool:
    return bool(manifests) and all(bool(item.get(key)) for item in manifests)


def _truthy(value: object) -> bool:
    return bool(value)


def _manifest_provenance_valid(manifest: dict[str, object]) -> bool:
    provenance = manifest.get("handoff_provenance")
    if not isinstance(provenance, dict):
        return False
    try:
        metadata = HandoffProvenanceMetadata(
            provenance_type=str(provenance.get("provenance_type", "")),
            design_source_label=str(provenance.get("design_source_label", "")),
            design_source_path=str(provenance.get("design_source_path", "")),
            design_digest=str(provenance.get("design_digest", "")),
            created_by_milestone=str(provenance.get("created_by_milestone", "")),
            frozen_before_stage=str(provenance.get("frozen_before_stage", "")),
            allowed_data_dependencies=list(
                provenance.get("allowed_data_dependencies", [])
            ),
            forbidden_data_dependencies=list(
                provenance.get("forbidden_data_dependencies", [])
            ),
            constraint_selection_policy=str(
                provenance.get("constraint_selection_policy", "")
            ),
            template_id=str(provenance.get("template_id", "")),
            template_hash=str(provenance.get("template_hash", "")),
            profile_instantiation_required=bool(
                provenance.get("profile_instantiation_required", False)
            ),
            fit_outputs_used=bool(provenance.get("fit_outputs_used", False)),
            heldout_outputs_used=bool(provenance.get("heldout_outputs_used", False)),
            carry_forward_outputs_used=bool(
                provenance.get("carry_forward_outputs_used", False)
            ),
            stress_test_outputs_used=bool(
                provenance.get("stress_test_outputs_used", False)
            ),
            report_only=bool(provenance.get("report_only", False)),
            reason_if_report_only=str(provenance.get("reason_if_report_only", "")),
        )
    except ValueError:
        return False
    return float(validate_handoff_provenance(metadata)["valid_provenance"]) == 1.0
