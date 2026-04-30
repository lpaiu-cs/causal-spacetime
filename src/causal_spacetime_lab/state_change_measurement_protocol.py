"""Measurement-protocol metadata for response-profile construction."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from typing import TypeAlias

ProtocolParameterValue: TypeAlias = str | int | float | bool

ECHO_RULES = {
    "earliest_return",
    "gated_return",
    "immediate_edge_return",
}
EMISSION_RULES = {
    "fixed_position",
    "central_position",
    "declared_emission_set",
}
GATE_RULES = {
    "none",
    "fixed_min_delay",
    "declared_gate",
}
REFERENCE_SUBSAMPLING_RULES = {
    "none",
    "fixed_stride",
    "retained_reference_set",
}
SPECTRUM_TYPES = {
    "full_transitive",
    "retained_reference",
    "immediate_edge",
}
NORMALIZATION_RULES = {
    "none",
    "per_protocol_range",
    "declared_scale",
}
MISSING_POLICIES = {
    "common_reachable",
    "penalize_mismatch",
    "require_all_reachable",
}
TIE_POLICIES = {
    "tie_as_unresolved",
    "tie_as_equal",
    "declared_tolerance",
}
MARGIN_POLICIES = {
    "fixed_margin",
    "zero_margin",
    "declared_margin_schedule",
}


@dataclass(frozen=True)
class MeasurementProtocolSpec:
    """Fixed measurement rules for a multi-reference response profile."""

    echo_rule: str
    emission_rule: str
    gate_rule: str
    reference_subsampling_rule: str
    spectrum_type: str
    normalization_rule: str
    missing_policy: str
    tie_policy: str
    margin_policy: str
    protocol_label: str = ""
    protocol_parameters: dict[str, ProtocolParameterValue] | None = None

    def __post_init__(self) -> None:
        _validate_allowed("echo_rule", self.echo_rule, ECHO_RULES)
        _validate_allowed("emission_rule", self.emission_rule, EMISSION_RULES)
        _validate_allowed("gate_rule", self.gate_rule, GATE_RULES)
        _validate_allowed(
            "reference_subsampling_rule",
            self.reference_subsampling_rule,
            REFERENCE_SUBSAMPLING_RULES,
        )
        _validate_allowed("spectrum_type", self.spectrum_type, SPECTRUM_TYPES)
        _validate_allowed(
            "normalization_rule",
            self.normalization_rule,
            NORMALIZATION_RULES,
        )
        _validate_allowed("missing_policy", self.missing_policy, MISSING_POLICIES)
        _validate_allowed("tie_policy", self.tie_policy, TIE_POLICIES)
        _validate_allowed("margin_policy", self.margin_policy, MARGIN_POLICIES)


def _validate_allowed(field_name: str, value: str, allowed_values: set[str]) -> None:
    if value not in allowed_values:
        allowed = ", ".join(sorted(allowed_values))
        raise ValueError(f"{field_name} must be one of: {allowed}")


def measurement_protocol_jsonable(
    spec: MeasurementProtocolSpec,
) -> dict[str, object]:
    """Return JSON-compatible measurement-protocol metadata."""

    payload = asdict(spec)
    payload["protocol_parameters"] = normalized_protocol_parameters(
        spec.protocol_parameters
    )
    return payload


def _hash_payload(spec: MeasurementProtocolSpec) -> dict[str, object]:
    payload = measurement_protocol_jsonable(spec)
    payload.pop("protocol_label", None)
    return payload


def measurement_protocol_hash(spec: MeasurementProtocolSpec) -> str:
    """Return a protocol hash independent of reference-chain identifiers."""

    payload = json.dumps(_hash_payload(spec), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def measurement_protocol_id(spec: MeasurementProtocolSpec) -> str:
    """Return a short stable protocol identifier."""

    label = spec.protocol_label.strip()
    digest = measurement_protocol_hash(spec)[:12]
    if not label:
        return f"measurement_protocol_{digest}"
    safe = label.lower().replace(" ", "_").replace("/", "_")
    return f"{safe}_{digest}"


def same_measurement_protocol(
    a: MeasurementProtocolSpec,
    b: MeasurementProtocolSpec,
) -> bool:
    """Return whether two specs use the same measurement protocol."""

    return measurement_protocol_hash(a) == measurement_protocol_hash(b)


def default_earliest_full_protocol(
    *,
    missing_policy: str = "common_reachable",
    tie_policy: str = "tie_as_unresolved",
    margin_policy: str = "fixed_margin",
    margin_value: float = 0.05,
) -> MeasurementProtocolSpec:
    """Return the default earliest-return full-spectrum protocol."""

    return MeasurementProtocolSpec(
        echo_rule="earliest_return",
        emission_rule="fixed_position",
        gate_rule="none",
        reference_subsampling_rule="none",
        spectrum_type="full_transitive",
        normalization_rule="per_protocol_range",
        missing_policy=missing_policy,
        tie_policy=tie_policy,
        margin_policy=margin_policy,
        protocol_label="earliest_full",
        protocol_parameters={
            "emission_position_rule": "declared_fixed_position",
            "normalization_scope": "profile_family",
            "margin_value": float(margin_value),
        },
    )


def default_gated_full_protocol(
    *,
    missing_policy: str = "common_reachable",
    tie_policy: str = "tie_as_unresolved",
    margin_policy: str = "fixed_margin",
    gate_delay_rank: int = 1,
    margin_value: float = 0.05,
) -> MeasurementProtocolSpec:
    """Return the default gated-return full-spectrum protocol."""

    return MeasurementProtocolSpec(
        echo_rule="gated_return",
        emission_rule="fixed_position",
        gate_rule="fixed_min_delay",
        reference_subsampling_rule="none",
        spectrum_type="full_transitive",
        normalization_rule="per_protocol_range",
        missing_policy=missing_policy,
        tie_policy=tie_policy,
        margin_policy=margin_policy,
        protocol_label="gated_full",
        protocol_parameters={
            "emission_position_rule": "declared_fixed_position",
            "gate_delay_rank": int(gate_delay_rank),
            "normalization_scope": "profile_family",
            "margin_value": float(margin_value),
        },
    )


def default_earliest_retained_protocol(
    *,
    missing_policy: str = "common_reachable",
    tie_policy: str = "tie_as_unresolved",
    margin_policy: str = "fixed_margin",
    retained_reference_policy: str = "declared_retained_reference_set",
    margin_value: float = 0.05,
) -> MeasurementProtocolSpec:
    """Return the default earliest-return retained-reference protocol."""

    return MeasurementProtocolSpec(
        echo_rule="earliest_return",
        emission_rule="fixed_position",
        gate_rule="none",
        reference_subsampling_rule="retained_reference_set",
        spectrum_type="retained_reference",
        normalization_rule="per_protocol_range",
        missing_policy=missing_policy,
        tie_policy=tie_policy,
        margin_policy=margin_policy,
        protocol_label="earliest_retained",
        protocol_parameters={
            "emission_position_rule": "declared_fixed_position",
            "retained_reference_policy": retained_reference_policy,
            "normalization_scope": "profile_family",
            "margin_value": float(margin_value),
        },
    )


def default_immediate_edge_protocol(
    *,
    missing_policy: str = "common_reachable",
    tie_policy: str = "tie_as_unresolved",
    margin_policy: str = "fixed_margin",
    margin_value: float = 0.05,
) -> MeasurementProtocolSpec:
    """Return the default immediate-edge response protocol."""

    return MeasurementProtocolSpec(
        echo_rule="immediate_edge_return",
        emission_rule="fixed_position",
        gate_rule="none",
        reference_subsampling_rule="none",
        spectrum_type="immediate_edge",
        normalization_rule="per_protocol_range",
        missing_policy=missing_policy,
        tie_policy=tie_policy,
        margin_policy=margin_policy,
        protocol_label="immediate_edge",
        protocol_parameters={
            "emission_position_rule": "declared_fixed_position",
            "normalization_scope": "profile_family",
            "margin_value": float(margin_value),
        },
    )


def normalized_protocol_parameters(
    parameters: dict[str, ProtocolParameterValue] | None,
) -> dict[str, ProtocolParameterValue]:
    """Return canonical protocol parameters used by protocol hashing."""

    if not parameters:
        return {}
    normalized: dict[str, ProtocolParameterValue] = {}
    for key in sorted(parameters):
        value = parameters[key]
        if isinstance(value, bool):
            normalized[str(key)] = bool(value)
        elif isinstance(value, int) and not isinstance(value, bool):
            normalized[str(key)] = int(value)
        elif isinstance(value, float):
            normalized[str(key)] = float(value)
        else:
            normalized[str(key)] = str(value)
    return normalized


def measurement_protocol_parameter_table(
    spec: MeasurementProtocolSpec,
) -> list[dict[str, str]]:
    """Return protocol parameters as CSV-friendly rows."""

    parameters = normalized_protocol_parameters(spec.protocol_parameters)
    return [
        {
            "measurement_protocol_id": measurement_protocol_id(spec),
            "parameter_name": key,
            "parameter_value": str(value),
        }
        for key, value in parameters.items()
    ]


def parameter_complete_for_protocol(spec: MeasurementProtocolSpec) -> bool:
    """Return whether required protocol-rule parameters are declared."""

    return not missing_required_protocol_parameters(spec)


def missing_required_protocol_parameters(spec: MeasurementProtocolSpec) -> list[str]:
    """Return required parameter names missing from a protocol spec."""

    parameters = normalized_protocol_parameters(spec.protocol_parameters)
    required: list[str] = []
    if spec.gate_rule == "fixed_min_delay":
        required.append("gate_delay_rank")
    if spec.emission_rule == "fixed_position":
        required.extend(["emission_position_rule|emission_position"])
    if spec.reference_subsampling_rule == "fixed_stride":
        required.append("reference_stride")
    if spec.reference_subsampling_rule == "retained_reference_set":
        required.append("retained_reference_policy")
    if spec.normalization_rule == "per_protocol_range":
        required.append("normalization_scope")
    if spec.margin_policy == "fixed_margin":
        required.append("margin_value")

    missing: list[str] = []
    for item in required:
        if "|" in item:
            choices = item.split("|")
            if not any(choice in parameters for choice in choices):
                missing.append(item)
        elif item not in parameters:
            missing.append(item)
    return missing
