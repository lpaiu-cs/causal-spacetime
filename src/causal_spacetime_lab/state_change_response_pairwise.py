"""Pairwise response-profile comparison protocols.

The dissimilarities here are chosen pre-metric protocol outputs. They are not
physical distances and do not define spatial geometry.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_response_profiles import (
    EchoResponseProfile,
    EchoResponseProfileWithMetadata,
)

PAIRWISE_METHODS = {
    "separation_fraction",
    "rank_gap_mean",
    "rank_gap_median",
    "reachability_mismatch",
    "combined_gap_and_mismatch",
}
MISSING_POLICIES = {
    "common_reachable",
    "penalize_mismatch",
    "require_all_reachable",
}


@dataclass(frozen=True)
class PairwiseResponseComparisonProtocol:
    """Configuration for a pre-metric pairwise response-profile comparison."""

    name: str
    method: str
    missing_policy: str = "common_reachable"
    min_common_protocols: int = 1
    normalize_by_protocol_range: bool = True
    reachability_mismatch_penalty: float = 1.0

    def __post_init__(self) -> None:
        if self.method not in PAIRWISE_METHODS:
            allowed = ", ".join(sorted(PAIRWISE_METHODS))
            raise ValueError(f"method must be one of: {allowed}")
        if self.missing_policy not in MISSING_POLICIES:
            allowed = ", ".join(sorted(MISSING_POLICIES))
            raise ValueError(f"missing_policy must be one of: {allowed}")
        if self.min_common_protocols < 1:
            raise ValueError("min_common_protocols must be at least 1")
        if self.reachability_mismatch_penalty < 0:
            raise ValueError("reachability_mismatch_penalty must be nonnegative")


@dataclass(frozen=True)
class PairwiseResponseDissimilarity:
    """Flattened pairwise response-profile dissimilarity values."""

    target_event_ids: NDArray[np.int_]
    pair_indices: NDArray[np.int_]
    dissimilarity_values: NDArray[np.float64]
    valid_pair_mask: NDArray[np.bool_]
    protocol_name: str
    method: str


def unordered_target_pairs(target_event_ids: ArrayLike) -> NDArray[np.int_]:
    """Return unordered target row-index pairs."""

    targets = np.asarray(target_event_ids)
    if targets.ndim != 1:
        raise ValueError("target_event_ids must be one-dimensional")
    pairs = [
        (i, j)
        for i in range(max(0, targets.size - 1))
        for j in range(i + 1, targets.size)
    ]
    return np.asarray(pairs, dtype=int).reshape((-1, 2))


def _protocol_ranges(profile: EchoResponseProfile) -> NDArray[np.float64]:
    ranges = np.ones(profile.delay_rank_matrix.shape[1], dtype=float)
    for column in range(profile.delay_rank_matrix.shape[1]):
        reachable = profile.reachable_matrix[:, column]
        if np.any(reachable):
            values = profile.delay_rank_matrix[reachable, column].astype(float)
            span = float(np.max(values) - np.min(values))
            ranges[column] = span if span > 0 else 1.0
    return ranges


def _gap_values(
    profile: EchoResponseProfile,
    i: int,
    j: int,
    common: NDArray[np.bool_],
    ranges: NDArray[np.float64],
    normalize: bool,
) -> NDArray[np.float64]:
    gaps = np.abs(
        profile.delay_rank_matrix[i, common].astype(float)
        - profile.delay_rank_matrix[j, common].astype(float)
    )
    if normalize:
        gaps = gaps / ranges[common]
    return gaps


def _pair_value(
    profile: EchoResponseProfile,
    protocol: PairwiseResponseComparisonProtocol,
    i: int,
    j: int,
    ranges: NDArray[np.float64],
) -> float:
    reachable_i = profile.reachable_matrix[i]
    reachable_j = profile.reachable_matrix[j]
    common = reachable_i & reachable_j
    mismatch = reachable_i ^ reachable_j
    common_count = int(np.sum(common))
    mismatch_fraction = float(np.mean(mismatch)) if mismatch.size else 0.0
    all_reachable = bool(np.all(reachable_i & reachable_j))

    if protocol.missing_policy == "require_all_reachable" and not all_reachable:
        return float("nan")

    if protocol.method == "reachability_mismatch":
        return mismatch_fraction

    if common_count >= protocol.min_common_protocols:
        if protocol.method == "separation_fraction":
            base = float(
                np.mean(
                    profile.delay_rank_matrix[i, common]
                    != profile.delay_rank_matrix[j, common]
                )
            )
        else:
            gaps = _gap_values(
                profile,
                i,
                j,
                common,
                ranges,
                protocol.normalize_by_protocol_range,
            )
            if protocol.method == "rank_gap_median":
                base = float(np.median(gaps))
            else:
                base = float(np.mean(gaps))
    elif protocol.missing_policy == "penalize_mismatch" and mismatch_fraction > 0:
        base = 0.0
    else:
        return float("nan")

    if protocol.missing_policy == "penalize_mismatch" or (
        protocol.method == "combined_gap_and_mismatch"
    ):
        if protocol.method == "combined_gap_and_mismatch":
            return base + protocol.reachability_mismatch_penalty * mismatch_fraction
        if common_count < protocol.min_common_protocols:
            return protocol.reachability_mismatch_penalty * mismatch_fraction
    return base


def pairwise_response_dissimilarity_matrix(
    profile: EchoResponseProfile,
    protocol: PairwiseResponseComparisonProtocol,
) -> NDArray[np.float64]:
    """Return a response-profile dissimilarity matrix.

    Values are protocol dissimilarities, not distances. Invalid comparisons
    are ``NaN`` and the diagonal is zero.
    """

    n_targets = profile.target_event_ids.size
    matrix = np.full((n_targets, n_targets), np.nan, dtype=float)
    np.fill_diagonal(matrix, 0.0)
    ranges = _protocol_ranges(profile)
    for i in range(max(0, n_targets - 1)):
        for j in range(i + 1, n_targets):
            value = _pair_value(profile, protocol, i, j, ranges)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix


def pairwise_response_dissimilarity(
    profile: EchoResponseProfile,
    protocol: PairwiseResponseComparisonProtocol,
) -> PairwiseResponseDissimilarity:
    """Return flattened pairwise response-profile dissimilarities."""

    matrix = pairwise_response_dissimilarity_matrix(profile, protocol)
    pairs = unordered_target_pairs(profile.target_event_ids)
    values = np.asarray([matrix[i, j] for i, j in pairs], dtype=float)
    valid = np.isfinite(values)
    return PairwiseResponseDissimilarity(
        target_event_ids=profile.target_event_ids.astype(int),
        pair_indices=pairs,
        dissimilarity_values=values,
        valid_pair_mask=valid,
        protocol_name=protocol.name,
        method=protocol.method,
    )


def validate_profile_for_pairwise_dissimilarity(
    profile_or_wrapped: EchoResponseProfile | EchoResponseProfileWithMetadata,
) -> dict[str, float | str]:
    """Validate profile metadata before production pairwise comparison."""

    metadata = getattr(profile_or_wrapped, "metadata", None)
    if isinstance(profile_or_wrapped, EchoResponseProfileWithMetadata):
        metadata = profile_or_wrapped.metadata
    if metadata is None:
        return {
            "profile_invariance_status": "underspecified",
            "admissible_for_pairwise_dissimilarity": 0.0,
            "exploratory_mixed_context": 0.0,
            "reason_if_blocked": "missing measurement protocol metadata",
        }
    admissible = (
        metadata.profile_invariance_status == "protocol_invariant"
        and metadata.admissible_for_pairwise_dissimilarity
    )
    reason = metadata.reason_if_blocked
    if not admissible and not reason:
        reason = "profile metadata is not production admissible"
    return {
        "profile_invariance_status": metadata.profile_invariance_status,
        "admissible_for_pairwise_dissimilarity": float(admissible),
        "exploratory_mixed_context": float(metadata.exploratory_mixed_context),
        "reason_if_blocked": reason,
    }


def pairwise_response_dissimilarity_checked(
    profile_or_wrapped: EchoResponseProfile | EchoResponseProfileWithMetadata,
    protocol: PairwiseResponseComparisonProtocol,
    *,
    allow_exploratory_mixed_context: bool = False,
) -> PairwiseResponseDissimilarity:
    """Return pairwise response-profile dissimilarity after metadata checks."""

    validation = validate_profile_for_pairwise_dissimilarity(profile_or_wrapped)
    status = str(validation["profile_invariance_status"])
    admissible = bool(float(validation["admissible_for_pairwise_dissimilarity"]))
    exploratory = bool(float(validation["exploratory_mixed_context"]))
    if not admissible:
        if not (
            allow_exploratory_mixed_context
            and exploratory
            and status == "protocol_mixed"
        ):
            reason = validation["reason_if_blocked"]
            raise ValueError(
                "response profile is inadmissible for production pairwise "
                f"response-profile dissimilarity: {reason}"
            )
    profile = (
        profile_or_wrapped.profile
        if isinstance(profile_or_wrapped, EchoResponseProfileWithMetadata)
        else profile_or_wrapped
    )
    return pairwise_response_dissimilarity(profile, protocol)


def response_pair_comparison_constraints(
    dissimilarity: PairwiseResponseDissimilarity,
    num_constraints: int,
    *,
    tolerance: float = 0.0,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Sample response-comparison quadruplets over target row indices."""

    if num_constraints < 0:
        raise ValueError("num_constraints must be nonnegative")
    if tolerance < 0:
        raise ValueError("tolerance must be nonnegative")
    valid_indices = np.flatnonzero(dissimilarity.valid_pair_mask)
    rows: list[tuple[int, int, int, int]] = []
    for left in valid_indices:
        for right in valid_indices:
            if left == right:
                continue
            left_value = dissimilarity.dissimilarity_values[left]
            right_value = dissimilarity.dissimilarity_values[right]
            if left_value + tolerance < right_value:
                i, j = dissimilarity.pair_indices[left]
                k, ell = dissimilarity.pair_indices[right]
                rows.append((int(i), int(j), int(k), int(ell)))
    if not rows or num_constraints == 0:
        return np.empty((0, 4), dtype=int)
    rng = np.random.default_rng(seed)
    if len(rows) > num_constraints:
        selected = rng.choice(len(rows), size=num_constraints, replace=False)
        rows = [rows[int(index)] for index in selected]
    return np.asarray(rows, dtype=int).reshape((-1, 4))


def _pair_value_by_event_id(
    dissimilarity: PairwiseResponseDissimilarity,
) -> dict[tuple[int, int], float]:
    values: dict[tuple[int, int], float] = {}
    for pair, value, valid in zip(
        dissimilarity.pair_indices,
        dissimilarity.dissimilarity_values,
        dissimilarity.valid_pair_mask,
        strict=True,
    ):
        if not valid:
            continue
        a = int(dissimilarity.target_event_ids[int(pair[0])])
        b = int(dissimilarity.target_event_ids[int(pair[1])])
        values[(min(a, b), max(a, b))] = float(value)
    return values


def pairwise_response_order_inversion_rate(
    baseline: PairwiseResponseDissimilarity,
    variant: PairwiseResponseDissimilarity,
) -> float:
    """Compare pairwise response-comparison orders over common valid pairs."""

    values_a = _pair_value_by_event_id(baseline)
    values_b = _pair_value_by_event_id(variant)
    common_pairs = sorted(set(values_a) & set(values_b))
    if len(common_pairs) < 2:
        return float("nan")
    total = 0
    inversions = 0
    for outer in range(len(common_pairs) - 1):
        for inner in range(outer + 1, len(common_pairs)):
            pair_a = common_pairs[outer]
            pair_b = common_pairs[inner]
            sign_a = np.sign(values_a[pair_a] - values_a[pair_b])
            sign_b = np.sign(values_b[pair_a] - values_b[pair_b])
            if sign_a == 0 or sign_b == 0:
                continue
            total += 1
            inversions += int(sign_a != sign_b)
    return float(inversions / total) if total else 0.0
