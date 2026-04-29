"""Null-baseline validation for response-comparison constraint pools."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_pool import (
    ResponseComparisonConstraintPool,
    evaluate_constraint_pool_on_dissimilarity,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
)
from causal_spacetime_lab.state_change_response_pairwise_nulls import (
    permute_target_profiles,
    random_profile_with_same_marginals,
    shuffle_profile_delays_within_protocols,
    shuffle_profile_reachability_masks,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _null_profiles(
    profile: EchoResponseProfile,
    seed: int | None,
) -> list[tuple[str, EchoResponseProfile]]:
    return [
        ("shuffle_delays", shuffle_profile_delays_within_protocols(profile, seed)),
        (
            "shuffle_reachability",
            shuffle_profile_reachability_masks(
                profile,
                None if seed is None else seed + 11,
            ),
        ),
        (
            "permute_profiles",
            permute_target_profiles(profile, seed=None if seed is None else seed + 22),
        ),
        (
            "random_same_marginals",
            random_profile_with_same_marginals(
                profile,
                None if seed is None else seed + 33,
            ),
        ),
    ]


def _z_score(value: float, null_values: np.ndarray) -> float:
    if null_values.size == 0:
        return float("nan")
    mean = float(np.mean(null_values))
    std = float(np.std(null_values))
    if std == 0:
        if value == mean:
            return 0.0
        return float("inf") if value > mean else float("-inf")
    return float((value - mean) / std)


def evaluate_constraint_pool_against_nulls(
    profile: EchoResponseProfile,
    comparison_protocol: PairwiseResponseComparisonProtocol,
    pool: ResponseComparisonConstraintPool,
    *,
    null_repetitions: int = 20,
    seed: int | None = None,
) -> dict[str, float]:
    """Evaluate a fixed pool on structured and null response profiles."""

    if null_repetitions < 1:
        raise ValueError("null_repetitions must be at least 1")
    structured = evaluate_constraint_pool_on_dissimilarity(
        pool,
        pairwise_response_dissimilarity(profile, comparison_protocol),
    )
    null_agreements: list[float] = []
    for repetition in range(null_repetitions):
        base_seed = None if seed is None else seed + 1000 * repetition
        for _, null_profile in _null_profiles(profile, base_seed):
            evaluation = evaluate_constraint_pool_on_dissimilarity(
                pool,
                pairwise_response_dissimilarity(null_profile, comparison_protocol),
            )
            null_agreements.append(float(evaluation["agreement_fraction"]))
    values = np.asarray(null_agreements, dtype=float)
    structured_agreement = float(structured["agreement_fraction"])
    return {
        "structured_agreement_fraction": structured_agreement,
        "mean_null_agreement_fraction": float(np.mean(values)),
        "std_null_agreement_fraction": float(np.std(values)),
        "null_z_score": _z_score(structured_agreement, values),
        "best_null_agreement_fraction": float(np.max(values)),
        "null_repetition_count": float(null_repetitions),
    }


def null_separation_by_type(
    profile: EchoResponseProfile,
    comparison_protocol: PairwiseResponseComparisonProtocol,
    pool: ResponseComparisonConstraintPool,
    *,
    null_repetitions: int = 20,
    seed: int | None = None,
) -> list[dict[str, float | str]]:
    """Return null-separation rows grouped by null type."""

    if null_repetitions < 1:
        raise ValueError("null_repetitions must be at least 1")
    structured = evaluate_constraint_pool_on_dissimilarity(
        pool,
        pairwise_response_dissimilarity(profile, comparison_protocol),
    )
    by_type: dict[str, list[float]] = {}
    for repetition in range(null_repetitions):
        base_seed = None if seed is None else seed + 1000 * repetition
        for null_type, null_profile in _null_profiles(profile, base_seed):
            evaluation = evaluate_constraint_pool_on_dissimilarity(
                pool,
                pairwise_response_dissimilarity(null_profile, comparison_protocol),
            )
            by_type.setdefault(null_type, []).append(
                float(evaluation["agreement_fraction"])
            )

    structured_agreement = float(structured["agreement_fraction"])
    rows: list[dict[str, float | str]] = []
    for null_type, agreements in sorted(by_type.items()):
        values = np.asarray(agreements, dtype=float)
        rows.append(
            {
                "null_type": null_type,
                "structured_agreement_fraction": structured_agreement,
                "mean_null_agreement_fraction": float(np.mean(values)),
                "std_null_agreement_fraction": float(np.std(values)),
                "null_z_score": _z_score(structured_agreement, values),
                "best_null_agreement_fraction": float(np.max(values)),
                "null_repetition_count": float(null_repetitions),
            }
        )
    return rows
