"""Map v4 preregistered families to measurement and comparison protocols."""

from __future__ import annotations

from dataclasses import replace

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    V4ProtocolExecutionFamilySpec,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    default_earliest_full_protocol,
    default_earliest_retained_protocol,
    default_gated_full_protocol,
    default_immediate_edge_protocol,
    measurement_protocol_hash,
    measurement_protocol_id,
    parameter_complete_for_protocol,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)


def measurement_protocol_for_v4_family(
    spec: V4ProtocolExecutionFamilySpec,
) -> MeasurementProtocolSpec:
    """Return the parameter-complete measurement protocol for a v4 family."""

    margin = float(spec.margin_value)
    missing_policy = (
        "penalize_mismatch"
        if spec.comparison_method == "combined_gap_and_mismatch"
        else "common_reachable"
    )
    if spec.measurement_protocol_family == "earliest_retained_reference":
        protocol = default_earliest_retained_protocol(
            missing_policy=missing_policy,
            margin_policy=spec.margin_policy,
            margin_value=margin,
        )
    elif spec.measurement_protocol_family == "gated_full_reference":
        protocol = default_gated_full_protocol(
            missing_policy=missing_policy,
            margin_policy=spec.margin_policy,
            gate_delay_rank=2,
            margin_value=margin,
        )
    elif spec.measurement_protocol_family == "immediate_edge_reference":
        protocol = default_immediate_edge_protocol(
            missing_policy=missing_policy,
            margin_policy=spec.margin_policy,
            margin_value=margin,
        )
    else:
        protocol = default_earliest_full_protocol(
            missing_policy=missing_policy,
            margin_policy=spec.margin_policy,
            margin_value=margin,
        )
    parameters = dict(protocol.protocol_parameters or {})
    parameters["v4_profile_resolution_policy"] = spec.profile_resolution_policy
    parameters["v4_target_inclusion_policy"] = spec.target_inclusion_policy
    if spec.family_name == "rank_gap_high_resolution_reference_set_v4":
        parameters["reference_set_policy"] = "high_resolution_diverse_reference_set"
        parameters["reference_chain_diversity_level"] = 2
    protocol = replace(
        protocol,
        protocol_label=spec.family_name,
        protocol_parameters=parameters,
    )
    if not parameter_complete_for_protocol(protocol):
        raise ValueError(f"v4 protocol is parameter-incomplete: {spec.family_name}")
    return protocol


def comparison_protocol_for_v4_family(
    spec: V4ProtocolExecutionFamilySpec,
    measurement_protocol: MeasurementProtocolSpec,
) -> PairwiseResponseComparisonProtocol:
    """Return the pairwise comparison protocol for a v4 family."""

    method = (
        "combined_gap_and_mismatch"
        if spec.family_name.startswith("combined_")
        else "rank_gap_mean"
    )
    penalty = 0.5 if method == "combined_gap_and_mismatch" else 1.0
    protocol = PairwiseResponseComparisonProtocol(
        name=f"{spec.family_name}_{method}",
        method=method,
        missing_policy=measurement_protocol.missing_policy,
        min_common_protocols=1,
        reachability_mismatch_penalty=penalty,
    )
    if protocol.missing_policy != measurement_protocol.missing_policy:
        raise ValueError("comparison missing policy must match measurement protocol")
    return protocol


def v4_protocol_id(protocol: MeasurementProtocolSpec) -> str:
    """Return stable v4 measurement protocol id."""

    return measurement_protocol_id(protocol)


def v4_protocol_hash(protocol: MeasurementProtocolSpec) -> str:
    """Return stable v4 measurement protocol hash."""

    return measurement_protocol_hash(protocol)
