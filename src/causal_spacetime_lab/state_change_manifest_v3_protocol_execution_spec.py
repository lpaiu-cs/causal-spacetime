"""Load protocol-patched v3 execution specs from M40 preregistration."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    measurement_protocol_hash,
    measurement_protocol_id,
    parameter_complete_for_protocol,
)

EXECUTION_STATUSES = {"planned_only", "executed"}


@dataclass(frozen=True)
class V3ProtocolExecutionFamilySpec:
    """Execution spec for one protocol-invariant patched v3 family."""

    family_name: str
    family_kind: str
    measurement_protocol: MeasurementProtocolSpec
    measurement_protocol_id: str
    measurement_protocol_hash: str
    planned_manifest_count: int
    profile_family_semantics: str
    admissible_for_pairwise_dissimilarity: bool
    execution_status: str
    comparison_method: str
    min_margin: float
    handoff_provenance_type: str
    top_down_template_required: bool

    def __post_init__(self) -> None:
        if self.execution_status not in EXECUTION_STATUSES:
            allowed = ", ".join(sorted(EXECUTION_STATUSES))
            raise ValueError(f"execution_status must be one of: {allowed}")
        if (
            self.family_kind == "structured"
            and not self.admissible_for_pairwise_dissimilarity
        ):
            raise ValueError("structured v3 protocol specs must be admissible")
        if self.family_kind == "structured" and self.planned_manifest_count < 3:
            raise ValueError(
                "structured v3 protocol specs require at least 3 manifests"
            )
        if self.family_kind == "structured" and not parameter_complete_for_protocol(
            self.measurement_protocol
        ):
            raise ValueError("structured v3 protocol specs must be parameter-complete")


def load_v3_protocol_patched_preregistration(path: Path) -> dict[str, object]:
    """Load the M40 protocol-patched v3 preregistration JSON."""

    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return dict(payload) if isinstance(payload, dict) else {}


def load_v3_protocol_execution_specs(
    path: Path,
) -> list[V3ProtocolExecutionFamilySpec]:
    """Read execution specs from the M40 patched preregistration."""

    payload = load_v3_protocol_patched_preregistration(path)
    patches = payload.get("patched_families", [])
    if not isinstance(patches, list):
        return []
    specs: list[V3ProtocolExecutionFamilySpec] = []
    for item in patches:
        if not isinstance(item, dict):
            continue
        protocol = _protocol_from_patch(item)
        family_name = str(item.get("patched_family_name", ""))
        family_kind = str(item.get("family_kind", ""))
        comparison_method = (
            "combined_gap_and_mismatch"
            if family_name.startswith("combined_")
            else "rank_gap_mean"
        )
        provenance_type = (
            "report_only_control"
            if family_kind in {"failed_control", "report_only"}
            else "hybrid_template_instantiated_from_profile"
        )
        specs.append(
            V3ProtocolExecutionFamilySpec(
                family_name=family_name,
                family_kind=family_kind,
                measurement_protocol=protocol,
                measurement_protocol_id=measurement_protocol_id(protocol),
                measurement_protocol_hash=measurement_protocol_hash(protocol),
                planned_manifest_count=int(
                    float(item.get("planned_manifest_count", 1))
                ),
                profile_family_semantics=str(item.get("profile_family_semantics", "")),
                admissible_for_pairwise_dissimilarity=bool(
                    item.get("admissible_for_pairwise_dissimilarity", False)
                ),
                execution_status=str(item.get("execution_status", "planned_only")),
                comparison_method=comparison_method,
                min_margin=_margin_value(protocol),
                handoff_provenance_type=provenance_type,
                top_down_template_required=provenance_type
                != "bottom_up_profile_derived",
            )
        )
    return specs


def mark_v3_protocol_execution_specs_executed(
    specs: list[V3ProtocolExecutionFamilySpec],
) -> list[V3ProtocolExecutionFamilySpec]:
    """Return copies of execution specs marked executed."""

    return [replace(spec, execution_status="executed") for spec in specs]


def v3_protocol_execution_spec_table(
    specs: list[V3ProtocolExecutionFamilySpec],
) -> list[dict[str, float | str]]:
    """Return execution specs as CSV rows."""

    rows: list[dict[str, float | str]] = []
    for spec in specs:
        rows.append(
            {
                "family_name": spec.family_name,
                "family_kind": spec.family_kind,
                "measurement_protocol_id": spec.measurement_protocol_id,
                "measurement_protocol_hash": spec.measurement_protocol_hash,
                "planned_manifest_count": float(spec.planned_manifest_count),
                "profile_family_semantics": spec.profile_family_semantics,
                "admissible_for_pairwise_dissimilarity": float(
                    spec.admissible_for_pairwise_dissimilarity
                ),
                "execution_status": spec.execution_status,
                "comparison_method": spec.comparison_method,
                "min_margin": float(spec.min_margin),
                "handoff_provenance_type": spec.handoff_provenance_type,
                "top_down_template_required": float(spec.top_down_template_required),
                "parameter_complete": float(
                    parameter_complete_for_protocol(spec.measurement_protocol)
                ),
            }
        )
    return rows


def _protocol_from_patch(item: dict[str, object]) -> MeasurementProtocolSpec:
    parameters: dict[str, str | int | float | bool] = {
        "emission_position_rule": "declared_fixed_position",
        "normalization_scope": "profile_family",
        "margin_value": 0.05,
    }
    if str(item.get("gate_rule", "")) == "fixed_min_delay":
        parameters["gate_delay_rank"] = 1
    if str(item.get("reference_subsampling_rule", "")) == "retained_reference_set":
        parameters["retained_reference_policy"] = "declared_retained_reference_set"
    if str(item.get("reference_subsampling_rule", "")) == "fixed_stride":
        parameters["reference_stride"] = 1
    return MeasurementProtocolSpec(
        echo_rule=str(item["echo_rule"]),
        emission_rule=str(item["emission_rule"]),
        gate_rule=str(item["gate_rule"]),
        reference_subsampling_rule=str(item["reference_subsampling_rule"]),
        spectrum_type=str(item["spectrum_type"]),
        normalization_rule=str(item["normalization_rule"]),
        missing_policy=str(item["missing_policy"]),
        tie_policy=str(item["tie_policy"]),
        margin_policy=str(item["margin_policy"]),
        protocol_label=str(item.get("patched_family_name", "")),
        protocol_parameters=parameters,
    )


def _margin_value(protocol: MeasurementProtocolSpec) -> float:
    parameters = protocol.protocol_parameters or {}
    return float(parameters.get("margin_value", 0.05))
