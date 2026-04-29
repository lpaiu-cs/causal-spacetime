from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
    evaluate_constraint_pool_on_dissimilarity,
    filter_constraint_pool_by_margin,
    merge_constraint_pools,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseDissimilarity,
)


def _dissimilarity(values: list[float]) -> PairwiseResponseDissimilarity:
    return PairwiseResponseDissimilarity(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        pair_indices=np.asarray(
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            dtype=int,
        ),
        dissimilarity_values=np.asarray(values, dtype=float),
        valid_pair_mask=np.isfinite(np.asarray(values, dtype=float)),
        protocol_name="gap",
        method="rank_gap_mean",
    )


def test_build_constraint_pool_from_dissimilarity_deterministic() -> None:
    pool = build_constraint_pool_from_dissimilarity(
        _dissimilarity([0.1, 0.2, 0.4, 0.5, np.nan, 0.9]),
        max_constraints=100,
        min_margin=0.3,
        seed=0,
    )

    assert pool.constraints.shape[1] == 4
    assert pool.constraints.shape[0] == pool.margins.size
    assert np.all(pool.margins > 0.3)


def test_evaluate_constraint_pool_identical_agreement() -> None:
    dissimilarity = _dissimilarity([0.1, 0.2, 0.4, 0.5, np.nan, 0.9])
    pool = build_constraint_pool_from_dissimilarity(
        dissimilarity,
        max_constraints=100,
        min_margin=0.0,
        seed=0,
    )

    report = evaluate_constraint_pool_on_dissimilarity(pool, dissimilarity)

    assert report["agreement_fraction"] == 1.0
    assert report["inversion_fraction"] == 0.0


def test_evaluate_constraint_pool_inverted_case() -> None:
    baseline = _dissimilarity([0.1, 0.2, 0.4, 0.5, 0.7, 0.9])
    inverted = _dissimilarity([0.9, 0.7, 0.5, 0.4, 0.2, 0.1])
    pool = build_constraint_pool_from_dissimilarity(
        baseline,
        max_constraints=100,
        min_margin=0.1,
        seed=0,
    )

    report = evaluate_constraint_pool_on_dissimilarity(pool, inverted)

    assert report["agreement_fraction"] == 0.0
    assert report["inversion_fraction"] == 1.0


def test_filter_constraint_pool_by_margin() -> None:
    pool = build_constraint_pool_from_dissimilarity(
        _dissimilarity([0.1, 0.2, 0.4, 0.5, 0.7, 0.9]),
        max_constraints=100,
        min_margin=0.0,
        seed=0,
    )

    filtered = filter_constraint_pool_by_margin(pool, 0.5)

    assert filtered.constraints.shape[0] < pool.constraints.shape[0]
    assert np.all(filtered.margins >= 0.5)


def test_merge_constraint_pools() -> None:
    dissimilarity = _dissimilarity([0.1, 0.2, 0.4, 0.5, 0.7, 0.9])
    pool_a = build_constraint_pool_from_dissimilarity(
        dissimilarity,
        max_constraints=5,
        seed=0,
        source_label="a",
    )
    pool_b = build_constraint_pool_from_dissimilarity(
        dissimilarity,
        max_constraints=5,
        seed=1,
        source_label="b",
    )

    merged = merge_constraint_pools([pool_a, pool_b])

    assert merged.source_label == "merged"
    input_count = pool_a.constraints.shape[0] + pool_b.constraints.shape[0]
    assert merged.constraints.shape[0] <= input_count
