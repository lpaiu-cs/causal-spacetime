"""Metadata and admissibility checks for response profiles."""

from __future__ import annotations

from dataclasses import dataclass

from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    measurement_protocol_hash,
    measurement_protocol_id,
    measurement_protocol_jsonable,
    missing_required_protocol_parameters,
    parameter_complete_for_protocol,
)

PROFILE_INVARIANCE_STATUSES = {
    "protocol_invariant",
    "protocol_mixed",
    "underspecified",
}


@dataclass(frozen=True)
class ResponseProfileMetadata:
    """Measurement metadata for a multi-reference response profile."""

    profile_family_id: str
    profile_family_semantics: str
    measurement_protocol_id: str
    measurement_protocol_hash: str
    measurement_protocol: MeasurementProtocolSpec
    reference_set_id: str
    reference_chain_ids: list[str]
    profile_invariance_status: str
    admissible_for_pairwise_dissimilarity: bool
    exploratory_mixed_context: bool = False
    reason_if_blocked: str = ""

    def __post_init__(self) -> None:
        if self.profile_invariance_status not in PROFILE_INVARIANCE_STATUSES:
            allowed = ", ".join(sorted(PROFILE_INVARIANCE_STATUSES))
            raise ValueError(
                f"profile_invariance_status must be one of: {allowed}"
            )
        if self.profile_invariance_status == "protocol_invariant":
            if not self.admissible_for_pairwise_dissimilarity:
                raise ValueError(
                    "protocol_invariant profiles must be admissible for "
                    "pairwise response-profile dissimilarity"
                )
        if self.profile_invariance_status == "underspecified":
            if self.admissible_for_pairwise_dissimilarity:
                raise ValueError("underspecified profiles must be blocked")
        if (
            self.profile_invariance_status == "protocol_mixed"
            and self.admissible_for_pairwise_dissimilarity
        ):
            raise ValueError(
                "protocol_mixed profiles are report-only and not production "
                "admissible"
            )


def profile_metadata_jsonable(
    metadata: ResponseProfileMetadata,
) -> dict[str, object]:
    """Return JSON-compatible profile metadata."""

    return {
        "profile_family_id": metadata.profile_family_id,
        "profile_family_semantics": metadata.profile_family_semantics,
        "measurement_protocol_id": metadata.measurement_protocol_id,
        "measurement_protocol_hash": metadata.measurement_protocol_hash,
        "measurement_protocol": measurement_protocol_jsonable(
            metadata.measurement_protocol
        ),
        "reference_set_id": metadata.reference_set_id,
        "reference_chain_ids": list(metadata.reference_chain_ids),
        "profile_invariance_status": metadata.profile_invariance_status,
        "admissible_for_pairwise_dissimilarity": bool(
            metadata.admissible_for_pairwise_dissimilarity
        ),
        "exploratory_mixed_context": bool(metadata.exploratory_mixed_context),
        "reason_if_blocked": metadata.reason_if_blocked,
    }


def profile_metadata_from_protocols(
    profile_family_id: str,
    protocols: list[MeasurementProtocolSpec],
    reference_chain_ids: list[str],
    reference_set_id: str,
    exploratory_mixed_context: bool = False,
) -> ResponseProfileMetadata:
    """Build response-profile metadata from per-column protocol declarations."""

    if not protocols:
        fallback = _fallback_protocol()
        return ResponseProfileMetadata(
            profile_family_id=profile_family_id,
            profile_family_semantics="underspecified measurement protocol",
            measurement_protocol_id="underspecified",
            measurement_protocol_hash="",
            measurement_protocol=fallback,
            reference_set_id=reference_set_id,
            reference_chain_ids=list(reference_chain_ids),
            profile_invariance_status="underspecified",
            admissible_for_pairwise_dissimilarity=False,
            exploratory_mixed_context=exploratory_mixed_context,
            reason_if_blocked="missing measurement protocol metadata",
        )

    incomplete = [
        (
            measurement_protocol_id(protocol),
            missing_required_protocol_parameters(protocol),
        )
        for protocol in protocols
        if not parameter_complete_for_protocol(protocol)
    ]
    if incomplete:
        fallback = protocols[0]
        reason = ";".join(
            f"{protocol_id}:{','.join(missing)}"
            for protocol_id, missing in incomplete
        )
        return ResponseProfileMetadata(
            profile_family_id=profile_family_id,
            profile_family_semantics="underspecified measurement protocol parameters",
            measurement_protocol_id=measurement_protocol_id(fallback),
            measurement_protocol_hash=measurement_protocol_hash(fallback),
            measurement_protocol=fallback,
            reference_set_id=reference_set_id,
            reference_chain_ids=list(reference_chain_ids),
            profile_invariance_status="underspecified",
            admissible_for_pairwise_dissimilarity=False,
            exploratory_mixed_context=exploratory_mixed_context,
            reason_if_blocked=f"missing protocol parameters: {reason}",
        )

    if not reference_chain_ids or not reference_set_id:
        fallback = protocols[0]
        return ResponseProfileMetadata(
            profile_family_id=profile_family_id,
            profile_family_semantics="underspecified reference metadata",
            measurement_protocol_id=measurement_protocol_id(fallback),
            measurement_protocol_hash=measurement_protocol_hash(fallback),
            measurement_protocol=fallback,
            reference_set_id=reference_set_id,
            reference_chain_ids=list(reference_chain_ids),
            profile_invariance_status="underspecified",
            admissible_for_pairwise_dissimilarity=False,
            exploratory_mixed_context=exploratory_mixed_context,
            reason_if_blocked="missing reference_chain_ids or reference_set_id",
        )

    hashes = [measurement_protocol_hash(protocol) for protocol in protocols]
    unique_hashes = set(hashes)
    first_protocol = protocols[0]
    if len(unique_hashes) == 1:
        return ResponseProfileMetadata(
            profile_family_id=profile_family_id,
            profile_family_semantics=(
                "protocol-invariant multi-reference response profile"
            ),
            measurement_protocol_id=measurement_protocol_id(first_protocol),
            measurement_protocol_hash=hashes[0],
            measurement_protocol=first_protocol,
            reference_set_id=reference_set_id,
            reference_chain_ids=list(reference_chain_ids),
            profile_invariance_status="protocol_invariant",
            admissible_for_pairwise_dissimilarity=True,
            exploratory_mixed_context=exploratory_mixed_context,
            reason_if_blocked="",
        )

    return ResponseProfileMetadata(
        profile_family_id=profile_family_id,
        profile_family_semantics=(
            "protocol-mixed response profile; report-only unless explicitly "
            "exploratory"
        ),
        measurement_protocol_id="mixed",
        measurement_protocol_hash=";".join(sorted(unique_hashes)),
        measurement_protocol=first_protocol,
        reference_set_id=reference_set_id,
        reference_chain_ids=list(reference_chain_ids),
        profile_invariance_status="protocol_mixed",
        admissible_for_pairwise_dissimilarity=False,
        exploratory_mixed_context=exploratory_mixed_context,
        reason_if_blocked="measurement protocol varies inside one profile",
    )


def validate_profile_admissibility(
    metadata: ResponseProfileMetadata,
) -> dict[str, float | str]:
    """Return a compact profile-admissibility row."""

    protocol_hashes = {
        item
        for item in metadata.measurement_protocol_hash.split(";")
        if item and item != "mixed"
    }
    if metadata.profile_invariance_status == "protocol_invariant":
        mixed_count = 1
    elif metadata.profile_invariance_status == "protocol_mixed":
        mixed_count = max(2, len(protocol_hashes))
    else:
        mixed_count = 0
    return {
        "profile_family_id": metadata.profile_family_id,
        "profile_invariance_status": metadata.profile_invariance_status,
        "admissible_for_pairwise_dissimilarity": float(
            metadata.admissible_for_pairwise_dissimilarity
        ),
        "number_of_protocols_mixed": float(mixed_count),
        "number_of_references": float(len(metadata.reference_chain_ids)),
        "reason_if_blocked": metadata.reason_if_blocked,
    }


def profile_parameter_completeness_report(
    metadata: ResponseProfileMetadata,
) -> dict[str, float | str]:
    """Report parameter completeness for response-profile metadata."""

    missing = missing_required_protocol_parameters(metadata.measurement_protocol)
    complete = (
        not missing
        and bool(metadata.reference_chain_ids)
        and bool(metadata.reference_set_id)
    )
    reason = ""
    if missing:
        reason = "missing protocol parameters: " + ",".join(missing)
    elif not metadata.reference_chain_ids:
        reason = "missing reference_chain_ids"
    elif not metadata.reference_set_id:
        reason = "missing reference_set_id"
    return {
        "profile_family_id": metadata.profile_family_id,
        "parameter_complete": float(complete),
        "missing_parameter_reason": reason,
        "profile_invariance_status": metadata.profile_invariance_status,
        "admissible_for_pairwise_dissimilarity": float(
            metadata.admissible_for_pairwise_dissimilarity
        ),
    }


def _fallback_protocol() -> MeasurementProtocolSpec:
    return MeasurementProtocolSpec(
        echo_rule="earliest_return",
        emission_rule="fixed_position",
        gate_rule="none",
        reference_subsampling_rule="none",
        spectrum_type="full_transitive",
        normalization_rule="none",
        missing_policy="common_reachable",
        tie_policy="tie_as_unresolved",
        margin_policy="fixed_margin",
        protocol_label="underspecified",
    )
