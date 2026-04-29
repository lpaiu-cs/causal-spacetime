"""Response-comparison constraint-pool utilities.

The constraints in this module compare pairwise response-profile
dissimilarities. They are pre-metric response-comparison constraints, not
distance-order constraints.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseDissimilarity,
)


@dataclass(frozen=True)
class ResponseComparisonConstraintPool:
    """Pool of response-comparison constraints over target row indices."""

    target_event_ids: NDArray[np.int_]
    constraints: NDArray[np.int_]
    margins: NDArray[np.float64]
    protocol_name: str
    method: str
    source_label: str


def _dissimilarity_matrix_by_row(
    dissimilarity: PairwiseResponseDissimilarity,
) -> NDArray[np.float64]:
    n_targets = dissimilarity.target_event_ids.size
    matrix = np.full((n_targets, n_targets), np.nan, dtype=float)
    np.fill_diagonal(matrix, 0.0)
    for pair, value, valid in zip(
        dissimilarity.pair_indices,
        dissimilarity.dissimilarity_values,
        dissimilarity.valid_pair_mask,
        strict=True,
    ):
        if not valid:
            continue
        i = int(pair[0])
        j = int(pair[1])
        matrix[i, j] = float(value)
        matrix[j, i] = float(value)
    return matrix


def build_constraint_pool_from_dissimilarity(
    dissimilarity: PairwiseResponseDissimilarity,
    *,
    max_constraints: int,
    min_margin: float = 0.0,
    seed: int | None = None,
    source_label: str = "baseline",
) -> ResponseComparisonConstraintPool:
    """Build a sampled high-margin response-comparison constraint pool."""

    if max_constraints < 0:
        raise ValueError("max_constraints must be nonnegative")
    if min_margin < 0:
        raise ValueError("min_margin must be nonnegative")
    valid_pair_indices = np.flatnonzero(dissimilarity.valid_pair_mask)
    if max_constraints == 0 or valid_pair_indices.size < 2:
        return ResponseComparisonConstraintPool(
            target_event_ids=dissimilarity.target_event_ids.astype(int),
            constraints=np.empty((0, 4), dtype=int),
            margins=np.empty(0, dtype=float),
            protocol_name=dissimilarity.protocol_name,
            method=dissimilarity.method,
            source_label=source_label,
        )
    if valid_pair_indices.size * valid_pair_indices.size > 250_000:
        return _sample_constraint_pool_from_dissimilarity(
            dissimilarity,
            valid_pair_indices,
            max_constraints=max_constraints,
            min_margin=min_margin,
            seed=seed,
            source_label=source_label,
        )
    rows: list[tuple[int, int, int, int]] = []
    margins: list[float] = []
    for left in valid_pair_indices:
        left_value = float(dissimilarity.dissimilarity_values[left])
        for right in valid_pair_indices:
            if left == right:
                continue
            right_value = float(dissimilarity.dissimilarity_values[right])
            margin = right_value - left_value
            if margin > min_margin:
                i, j = dissimilarity.pair_indices[left]
                k, ell = dissimilarity.pair_indices[right]
                rows.append((int(i), int(j), int(k), int(ell)))
                margins.append(margin)

    if not rows:
        constraints = np.empty((0, 4), dtype=int)
        margin_array = np.empty(0, dtype=float)
    else:
        rng = np.random.default_rng(seed)
        selected = np.arange(len(rows), dtype=int)
        if len(rows) > max_constraints:
            selected = rng.choice(len(rows), size=max_constraints, replace=False)
        constraints = np.asarray([rows[int(index)] for index in selected], dtype=int)
        margin_array = np.asarray(
            [margins[int(index)] for index in selected],
            dtype=float,
        )

    return ResponseComparisonConstraintPool(
        target_event_ids=dissimilarity.target_event_ids.astype(int),
        constraints=constraints.reshape((-1, 4)),
        margins=margin_array,
        protocol_name=dissimilarity.protocol_name,
        method=dissimilarity.method,
        source_label=source_label,
    )


def _sample_constraint_pool_from_dissimilarity(
    dissimilarity: PairwiseResponseDissimilarity,
    valid_pair_indices: NDArray[np.int_],
    *,
    max_constraints: int,
    min_margin: float,
    seed: int | None,
    source_label: str,
) -> ResponseComparisonConstraintPool:
    """Sample a large constraint pool without enumerating all pair pairs."""

    rng = np.random.default_rng(seed)
    valid_values = dissimilarity.dissimilarity_values[valid_pair_indices].astype(float)
    order = np.argsort(valid_values, kind="mergesort")
    sorted_values = valid_values[order]
    sorted_pair_indices = valid_pair_indices[order]
    constraints: list[tuple[int, int, int, int]] = []
    margins: list[float] = []
    seen: set[tuple[int, int, int, int]] = set()
    attempts = max(10_000, max_constraints * 50)
    for _ in range(attempts):
        if len(constraints) >= max_constraints:
            break
        left_order_index = int(rng.integers(0, sorted_values.size - 1))
        threshold = sorted_values[left_order_index] + min_margin
        first_right = int(np.searchsorted(sorted_values, threshold, side="right"))
        if first_right >= sorted_values.size:
            continue
        right_order_index = int(rng.integers(first_right, sorted_values.size))
        left_pair_index = int(sorted_pair_indices[left_order_index])
        right_pair_index = int(sorted_pair_indices[right_order_index])
        if left_pair_index == right_pair_index:
            continue
        i, j = dissimilarity.pair_indices[left_pair_index]
        k, ell = dissimilarity.pair_indices[right_pair_index]
        key = (int(i), int(j), int(k), int(ell))
        if key in seen:
            continue
        seen.add(key)
        constraints.append(key)
        margins.append(
            float(sorted_values[right_order_index] - sorted_values[left_order_index])
        )
    return ResponseComparisonConstraintPool(
        target_event_ids=dissimilarity.target_event_ids.astype(int),
        constraints=np.asarray(constraints, dtype=int).reshape((-1, 4)),
        margins=np.asarray(margins, dtype=float),
        protocol_name=dissimilarity.protocol_name,
        method=dissimilarity.method,
        source_label=source_label,
    )


def _pool_rows_to_dissimilarity_rows(
    pool: ResponseComparisonConstraintPool,
    dissimilarity: PairwiseResponseDissimilarity,
) -> dict[int, int]:
    target_lookup = {
        int(target_id): int(index)
        for index, target_id in enumerate(dissimilarity.target_event_ids)
    }
    mapping: dict[int, int] = {}
    for pool_row, target_id in enumerate(pool.target_event_ids):
        if int(target_id) in target_lookup:
            mapping[pool_row] = target_lookup[int(target_id)]
    return mapping


def evaluate_constraint_pool_on_dissimilarity(
    pool: ResponseComparisonConstraintPool,
    dissimilarity: PairwiseResponseDissimilarity,
    *,
    tolerance: float = 0.0,
) -> dict[str, float]:
    """Evaluate a response-comparison pool on another dissimilarity table."""

    if tolerance < 0:
        raise ValueError("tolerance must be nonnegative")
    constraint_count = int(pool.constraints.shape[0])
    if constraint_count == 0:
        return {
            "constraint_count": 0.0,
            "evaluable_count": 0.0,
            "evaluable_fraction": 0.0,
            "agreement_count": 0.0,
            "inversion_count": 0.0,
            "tie_or_unresolved_count": 0.0,
            "agreement_fraction": 0.0,
            "inversion_fraction": 0.0,
            "tie_or_unresolved_fraction": 0.0,
        }

    row_mapping = _pool_rows_to_dissimilarity_rows(pool, dissimilarity)
    matrix = _dissimilarity_matrix_by_row(dissimilarity)
    evaluable_count = 0
    agreement_count = 0
    inversion_count = 0
    tie_or_unresolved_count = 0
    for i, j, k, ell in pool.constraints:
        if (
            int(i) not in row_mapping
            or int(j) not in row_mapping
            or int(k) not in row_mapping
            or int(ell) not in row_mapping
        ):
            tie_or_unresolved_count += 1
            continue
        left = matrix[row_mapping[int(i)], row_mapping[int(j)]]
        right = matrix[row_mapping[int(k)], row_mapping[int(ell)]]
        if not np.isfinite(left) or not np.isfinite(right):
            tie_or_unresolved_count += 1
            continue
        evaluable_count += 1
        if left + tolerance < right:
            agreement_count += 1
        elif left > right + tolerance:
            inversion_count += 1
        else:
            tie_or_unresolved_count += 1

    return {
        "constraint_count": float(constraint_count),
        "evaluable_count": float(evaluable_count),
        "evaluable_fraction": float(evaluable_count / constraint_count),
        "agreement_count": float(agreement_count),
        "inversion_count": float(inversion_count),
        "tie_or_unresolved_count": float(tie_or_unresolved_count),
        "agreement_fraction": float(agreement_count / constraint_count),
        "inversion_fraction": float(inversion_count / constraint_count),
        "tie_or_unresolved_fraction": float(
            tie_or_unresolved_count / constraint_count
        ),
    }


def filter_constraint_pool_by_margin(
    pool: ResponseComparisonConstraintPool,
    min_margin: float,
) -> ResponseComparisonConstraintPool:
    """Return constraints with margin at least ``min_margin``."""

    if min_margin < 0:
        raise ValueError("min_margin must be nonnegative")
    keep = pool.margins >= min_margin
    return ResponseComparisonConstraintPool(
        target_event_ids=pool.target_event_ids.copy(),
        constraints=pool.constraints[keep].copy(),
        margins=pool.margins[keep].copy(),
        protocol_name=pool.protocol_name,
        method=pool.method,
        source_label=pool.source_label,
    )


def merge_constraint_pools(
    pools: list[ResponseComparisonConstraintPool],
    source_label: str = "merged",
) -> ResponseComparisonConstraintPool:
    """Merge pools and deduplicate identical constraints."""

    if not pools:
        return ResponseComparisonConstraintPool(
            target_event_ids=np.empty(0, dtype=int),
            constraints=np.empty((0, 4), dtype=int),
            margins=np.empty(0, dtype=float),
            protocol_name="mixed",
            method="mixed",
            source_label=source_label,
        )

    target_ids = pools[0].target_event_ids.astype(int)
    for pool in pools[1:]:
        if not np.array_equal(target_ids, pool.target_event_ids.astype(int)):
            raise ValueError("all pools must use the same target_event_ids")

    protocol_names = {pool.protocol_name for pool in pools}
    methods = {pool.method for pool in pools}
    best_margins: dict[tuple[int, int, int, int], float] = {}
    for pool in pools:
        for constraint, margin in zip(pool.constraints, pool.margins, strict=True):
            key = tuple(int(value) for value in constraint)
            best_margins[key] = max(float(margin), best_margins.get(key, -np.inf))

    constraints = np.asarray(list(best_margins), dtype=int).reshape((-1, 4))
    margins = np.asarray([best_margins[key] for key in best_margins], dtype=float)
    protocol_name = next(iter(protocol_names)) if len(protocol_names) == 1 else "mixed"
    return ResponseComparisonConstraintPool(
        target_event_ids=target_ids.copy(),
        constraints=constraints,
        margins=margins,
        protocol_name=protocol_name,
        method=next(iter(methods)) if len(methods) == 1 else "mixed",
        source_label=source_label,
    )
