"""Coverage diagnostics for response-comparison constraint pools."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_pool import (
    ResponseComparisonConstraintPool,
)


def constraint_target_coverage(
    pool: ResponseComparisonConstraintPool,
    target_count: int,
) -> dict[str, float]:
    """Summarize how many targets participate in a constraint pool."""

    if target_count < 0:
        raise ValueError("target_count must be nonnegative")
    counts = np.zeros(target_count, dtype=int)
    for constraint in pool.constraints:
        for target in constraint:
            if 0 <= int(target) < target_count:
                counts[int(target)] += 1
    touched = counts > 0
    touched_count = int(np.sum(touched))
    return {
        "target_count": float(target_count),
        "touched_target_count": float(touched_count),
        "touched_target_fraction": float(touched_count / target_count)
        if target_count
        else 0.0,
        "mean_constraints_per_touched_target": float(np.mean(counts[touched]))
        if touched_count
        else 0.0,
        "min_constraints_per_touched_target": float(np.min(counts[touched]))
        if touched_count
        else 0.0,
        "max_constraints_per_target": float(np.max(counts)) if target_count else 0.0,
    }


def _pair_key(i: int, j: int) -> tuple[int, int]:
    return (min(int(i), int(j)), max(int(i), int(j)))


def constraint_pair_node_coverage(
    pool: ResponseComparisonConstraintPool,
    target_count: int,
) -> dict[str, float]:
    """Treat unordered target pairs as nodes and summarize coverage."""

    if target_count < 0:
        raise ValueError("target_count must be nonnegative")
    possible = target_count * (target_count - 1) // 2
    counts: dict[tuple[int, int], int] = {}
    for i, j, k, ell in pool.constraints:
        for pair in [_pair_key(int(i), int(j)), _pair_key(int(k), int(ell))]:
            valid_pair = (
                pair[0] != pair[1]
                and 0 <= pair[0] < target_count
                and pair[1] < target_count
            )
            if valid_pair:
                counts[pair] = counts.get(pair, 0) + 1
    touched = len(counts)
    return {
        "possible_pair_node_count": float(possible),
        "touched_pair_node_count": float(touched),
        "touched_pair_node_fraction": float(touched / possible) if possible else 0.0,
        "mean_constraints_per_touched_pair_node": float(np.mean(list(counts.values())))
        if counts
        else 0.0,
    }


def constraint_pool_summary(
    pool: ResponseComparisonConstraintPool,
    target_count: int,
) -> dict[str, float | str]:
    """Combine margin, target, and pair-node coverage diagnostics."""

    margin_count = pool.margins.size
    if margin_count:
        mean_margin = float(np.mean(pool.margins))
        median_margin = float(np.median(pool.margins))
        min_margin = float(np.min(pool.margins))
        max_margin = float(np.max(pool.margins))
    else:
        mean_margin = median_margin = min_margin = max_margin = 0.0
    return {
        "protocol_name": pool.protocol_name,
        "method": pool.method,
        "source_label": pool.source_label,
        "constraint_count": float(pool.constraints.shape[0]),
        "mean_margin": mean_margin,
        "median_margin": median_margin,
        "min_margin": min_margin,
        "max_margin": max_margin,
        **constraint_target_coverage(pool, target_count),
        **constraint_pair_node_coverage(pool, target_count),
    }
