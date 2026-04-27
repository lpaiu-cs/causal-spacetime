from __future__ import annotations

import numpy as np

from causal_spacetime_lab.embedding_stability import (
    fit_embeddings_on_constraint_subsets,
    pairwise_order_stability,
    pairwise_procrustes_stability,
)
from causal_spacetime_lab.ordinal_embedding import (
    sample_quadruplet_constraints_from_coords,
)


def test_fit_embeddings_on_constraint_subsets_returns_expected_count() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0], [6.0]])
    constraints = sample_quadruplet_constraints_from_coords(coords, 30, seed=1)

    embeddings = fit_embeddings_on_constraint_subsets(
        4,
        1,
        constraints,
        num_subsets=3,
        subset_size=10,
        seed=2,
        steps=10,
        restarts=1,
    )

    assert len(embeddings) == 3
    assert all(embedding.shape == (4, 1) for embedding in embeddings)


def test_pairwise_procrustes_stability_identical_embeddings_is_zero() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0]])

    report = pairwise_procrustes_stability([coords, coords.copy(), coords.copy()])

    assert report["mean_procrustes_rmse"] < 1e-12
    assert report["pair_count"] == 3.0


def test_pairwise_order_stability_identical_embeddings_is_zero() -> None:
    coords = np.asarray([[0.0], [1.0], [3.0], [6.0]])

    report = pairwise_order_stability([coords, coords.copy()], seed=3)

    assert report["mean_order_disagreement"] == 0.0
    assert report["pair_count"] == 1.0
