"""Validation diagnostics for response-comparison constraint pools."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_pool import (
    ResponseComparisonConstraintPool,
    build_constraint_pool_from_dissimilarity,
    evaluate_constraint_pool_on_dissimilarity,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    PairwiseResponseDissimilarity,
    pairwise_response_dissimilarity,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


@dataclass(frozen=True)
class ConstraintValidationGate:
    """Thresholds for pre-embedding response-comparison pool validation."""

    name: str
    min_constraint_count: int = 100
    min_evaluable_fraction: float = 0.5
    min_agreement_fraction: float = 0.7
    max_inversion_fraction: float = 0.2
    max_tie_or_unresolved_fraction: float = 0.5
    min_null_z_score: float = 1.0
    min_bootstrap_confidence: float = 0.7


def _profile_with_columns(
    profile: EchoResponseProfile,
    columns: np.ndarray,
) -> EchoResponseProfile:
    columns = np.asarray(columns, dtype=int)
    return EchoResponseProfile(
        target_event_ids=profile.target_event_ids.copy(),
        protocol_labels=[profile.protocol_labels[int(column)] for column in columns],
        delay_rank_matrix=profile.delay_rank_matrix[:, columns].copy(),
        reachable_matrix=profile.reachable_matrix[:, columns].copy(),
    )


def split_profile_protocol_columns(
    profile: EchoResponseProfile,
    train_fraction: float,
    seed: int | None = None,
) -> tuple[EchoResponseProfile, EchoResponseProfile]:
    """Split response-profile protocol columns into train and test profiles."""

    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must be between 0 and 1")
    protocol_count = len(profile.protocol_labels)
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(protocol_count)
    if protocol_count <= 1:
        train_count = protocol_count
    else:
        train_count = int(round(protocol_count * train_fraction))
        train_count = min(max(train_count, 1), protocol_count - 1)
    train_columns = np.sort(permutation[:train_count])
    test_columns = np.sort(permutation[train_count:])
    return (
        _profile_with_columns(profile, train_columns),
        _profile_with_columns(profile, test_columns),
    )


def build_pairwise_dissimilarity_for_profile(
    profile: EchoResponseProfile,
    protocol: PairwiseResponseComparisonProtocol,
) -> PairwiseResponseDissimilarity:
    """Wrap pairwise response-profile dissimilarity construction."""

    return pairwise_response_dissimilarity(profile, protocol)


def heldout_protocol_constraint_validation(
    profile: EchoResponseProfile,
    comparison_protocol: PairwiseResponseComparisonProtocol,
    *,
    train_fraction: float = 0.6,
    max_constraints: int = 5000,
    min_margin: float = 0.0,
    seed: int | None = None,
) -> dict[str, float]:
    """Build constraints on train protocol columns and evaluate held-out columns."""

    train_profile, test_profile = split_profile_protocol_columns(
        profile,
        train_fraction,
        seed,
    )
    train_dissimilarity = build_pairwise_dissimilarity_for_profile(
        train_profile,
        comparison_protocol,
    )
    test_dissimilarity = build_pairwise_dissimilarity_for_profile(
        test_profile,
        comparison_protocol,
    )
    pool = build_constraint_pool_from_dissimilarity(
        train_dissimilarity,
        max_constraints=max_constraints,
        min_margin=min_margin,
        seed=seed,
        source_label="heldout_train",
    )
    evaluation = evaluate_constraint_pool_on_dissimilarity(pool, test_dissimilarity)
    return {
        **evaluation,
        "train_protocol_count": float(len(train_profile.protocol_labels)),
        "test_protocol_count": float(len(test_profile.protocol_labels)),
        "train_valid_pair_fraction": float(np.mean(train_dissimilarity.valid_pair_mask))
        if train_dissimilarity.valid_pair_mask.size
        else 0.0,
        "test_valid_pair_fraction": float(np.mean(test_dissimilarity.valid_pair_mask))
        if test_dissimilarity.valid_pair_mask.size
        else 0.0,
    }


def validation_gate_pass_fail(
    validation_row: dict[str, float],
    gate: ConstraintValidationGate,
) -> dict[str, float | str]:
    """Apply gate thresholds to a validation row."""

    checks = [
        (
            "constraint_count",
            validation_row.get("constraint_count", float("nan")),
            ">=",
            float(gate.min_constraint_count),
        ),
        (
            "evaluable_fraction",
            validation_row.get("evaluable_fraction", float("nan")),
            ">=",
            gate.min_evaluable_fraction,
        ),
        (
            "agreement_fraction",
            validation_row.get("agreement_fraction", float("nan")),
            ">=",
            gate.min_agreement_fraction,
        ),
        (
            "inversion_fraction",
            validation_row.get("inversion_fraction", float("nan")),
            "<=",
            gate.max_inversion_fraction,
        ),
        (
            "tie_or_unresolved_fraction",
            validation_row.get("tie_or_unresolved_fraction", float("nan")),
            "<=",
            gate.max_tie_or_unresolved_fraction,
        ),
    ]
    optional_checks = [
        ("null_z_score", ">=", gate.min_null_z_score),
        ("mean_agreement_fraction", ">=", gate.min_bootstrap_confidence),
        ("stable_constraint_fraction", ">=", gate.min_bootstrap_confidence),
    ]
    for key, operator, threshold in optional_checks:
        if key in validation_row and np.isfinite(float(validation_row[key])):
            checks.append((key, validation_row[key], operator, threshold))

    failed: list[str] = []
    for key, value, operator, threshold in checks:
        value_float = float(value)
        if not np.isfinite(value_float):
            failed.append(f"{key}=nan")
        elif operator == ">=" and value_float < threshold:
            failed.append(f"{key}<{threshold:g}")
        elif operator == "<=" and value_float > threshold:
            failed.append(f"{key}>{threshold:g}")
    return {
        "gate_name": gate.name,
        "passed": float(not failed),
        "failed_reasons": ";".join(failed),
        "constraint_count": float(validation_row.get("constraint_count", np.nan)),
        "evaluable_fraction": float(
            validation_row.get("evaluable_fraction", np.nan),
        ),
        "agreement_fraction": float(
            validation_row.get("agreement_fraction", np.nan),
        ),
        "inversion_fraction": float(
            validation_row.get("inversion_fraction", np.nan),
        ),
        "tie_or_unresolved_fraction": float(
            validation_row.get("tie_or_unresolved_fraction", np.nan),
        ),
    }


def bootstrap_profile_protocol_columns(
    profile: EchoResponseProfile,
    *,
    sample_fraction: float = 0.8,
    seed: int | None = None,
) -> EchoResponseProfile:
    """Sample protocol columns with replacement."""

    if not sample_fraction > 0:
        raise ValueError("sample_fraction must be positive")
    protocol_count = len(profile.protocol_labels)
    sample_count = max(1, int(round(protocol_count * sample_fraction)))
    rng = np.random.default_rng(seed)
    columns = rng.choice(protocol_count, size=sample_count, replace=True)
    return EchoResponseProfile(
        target_event_ids=profile.target_event_ids.copy(),
        protocol_labels=[
            f"{profile.protocol_labels[int(column)]}#boot{index}"
            for index, column in enumerate(columns)
        ],
        delay_rank_matrix=profile.delay_rank_matrix[:, columns].copy(),
        reachable_matrix=profile.reachable_matrix[:, columns].copy(),
    )


def _constraint_agreement_flags(
    pool: ResponseComparisonConstraintPool,
    dissimilarity: PairwiseResponseDissimilarity,
    tolerance: float,
) -> np.ndarray:
    from causal_spacetime_lab.state_change_response_constraint_pool import (
        _dissimilarity_matrix_by_row,
        _pool_rows_to_dissimilarity_rows,
    )

    row_mapping = _pool_rows_to_dissimilarity_rows(pool, dissimilarity)
    matrix = _dissimilarity_matrix_by_row(dissimilarity)
    flags = np.full(pool.constraints.shape[0], np.nan, dtype=float)
    for index, (i, j, k, ell) in enumerate(pool.constraints):
        if (
            int(i) not in row_mapping
            or int(j) not in row_mapping
            or int(k) not in row_mapping
            or int(ell) not in row_mapping
        ):
            continue
        left = matrix[row_mapping[int(i)], row_mapping[int(j)]]
        right = matrix[row_mapping[int(k)], row_mapping[int(ell)]]
        if not np.isfinite(left) or not np.isfinite(right):
            continue
        flags[index] = 1.0 if left + tolerance < right else 0.0
    return flags


def bootstrap_constraint_stability(
    profile: EchoResponseProfile,
    comparison_protocol: PairwiseResponseComparisonProtocol,
    pool: ResponseComparisonConstraintPool,
    *,
    bootstrap_count: int = 50,
    sample_fraction: float = 0.8,
    tolerance: float = 0.0,
    seed: int | None = None,
) -> dict[str, float]:
    """Evaluate a fixed constraint pool across protocol-column bootstraps."""

    if bootstrap_count < 1:
        raise ValueError("bootstrap_count must be at least 1")
    aggregate_rows: list[dict[str, float]] = []
    flags_by_bootstrap: list[np.ndarray] = []
    for bootstrap_index in range(bootstrap_count):
        bootstrap_profile = bootstrap_profile_protocol_columns(
            profile,
            sample_fraction=sample_fraction,
            seed=None if seed is None else seed + bootstrap_index,
        )
        dissimilarity = build_pairwise_dissimilarity_for_profile(
            bootstrap_profile,
            comparison_protocol,
        )
        aggregate_rows.append(
            evaluate_constraint_pool_on_dissimilarity(
                pool,
                dissimilarity,
                tolerance=tolerance,
            )
        )
        flags_by_bootstrap.append(
            _constraint_agreement_flags(pool, dissimilarity, tolerance)
        )

    agreement_values = np.asarray(
        [row["agreement_fraction"] for row in aggregate_rows],
        dtype=float,
    )
    inversion_values = np.asarray(
        [row["inversion_fraction"] for row in aggregate_rows],
        dtype=float,
    )
    evaluable_values = np.asarray(
        [row["evaluable_fraction"] for row in aggregate_rows],
        dtype=float,
    )
    flags = np.vstack(flags_by_bootstrap) if flags_by_bootstrap else np.empty((0, 0))
    stable_count = 0
    for column in range(flags.shape[1]):
        finite = np.isfinite(flags[:, column])
        if np.any(finite) and float(np.mean(flags[finite, column])) >= 0.7:
            stable_count += 1
    stable_fraction = stable_count / flags.shape[1] if flags.shape[1] else 0.0
    return {
        "bootstrap_count": float(bootstrap_count),
        "mean_agreement_fraction": float(np.mean(agreement_values)),
        "std_agreement_fraction": float(np.std(agreement_values)),
        "mean_inversion_fraction": float(np.mean(inversion_values)),
        "mean_evaluable_fraction": float(np.mean(evaluable_values)),
        "stable_constraint_fraction": float(stable_fraction),
    }
