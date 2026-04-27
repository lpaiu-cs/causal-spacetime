"""Stability diagnostics for ordinal effective metric representations."""

from __future__ import annotations

import itertools

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.ordinal_embedding import (
    fit_ordinal_embedding_gradient_descent,
    procrustes_align,
    squared_distance_matrix,
)


def _as_constraints(constraints: ArrayLike) -> NDArray[np.int_]:
    array = np.asarray(constraints, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    return array


def fit_embeddings_on_constraint_subsets(
    n_points: int,
    embedding_dim: int,
    constraints: ArrayLike,
    num_subsets: int = 5,
    subset_size: int | None = None,
    seed: int | None = None,
    steps: int = 1000,
    restarts: int = 2,
    learning_rate: float = 0.05,
) -> list[NDArray[np.float64]]:
    """Fit embeddings on independent subsets of the same constraint pool."""

    constraint_array = _as_constraints(constraints)
    subset_count = int(num_subsets)
    if subset_count <= 0:
        raise ValueError("num_subsets must be positive")
    size = constraint_array.shape[0] if subset_size is None else int(subset_size)
    if size <= 0:
        raise ValueError("subset_size must be positive")
    replace = size > constraint_array.shape[0]
    rng = np.random.default_rng(seed)
    embeddings: list[NDArray[np.float64]] = []
    for subset_index in range(subset_count):
        selected = rng.choice(
            constraint_array.shape[0],
            size=size,
            replace=replace,
        )
        subset = constraint_array[selected]
        embedding, _ = fit_ordinal_embedding_gradient_descent(
            n_points,
            embedding_dim,
            subset,
            steps=steps,
            restarts=restarts,
            learning_rate=learning_rate,
            seed=None if seed is None else seed + 10_000 * subset_index,
            batch_size=min(2048, max(1, subset.shape[0])),
        )
        embeddings.append(embedding)
    return embeddings


def pairwise_procrustes_stability(
    embeddings: list[ArrayLike],
) -> dict[str, float]:
    """Return Procrustes RMSE statistics across fitted embeddings."""

    if len(embeddings) < 2:
        return {
            "mean_procrustes_rmse": float("nan"),
            "std_procrustes_rmse": float("nan"),
            "max_procrustes_rmse": float("nan"),
            "pair_count": 0.0,
        }
    values: list[float] = []
    for first, second in itertools.combinations(embeddings, 2):
        _, diagnostics = procrustes_align(first, second)
        values.append(diagnostics["rmse"])
    array = np.asarray(values, dtype=float)
    return {
        "mean_procrustes_rmse": float(np.mean(array)),
        "std_procrustes_rmse": float(np.std(array)),
        "max_procrustes_rmse": float(np.max(array)),
        "pair_count": float(array.size),
    }


def _order_disagreement(
    first: NDArray[np.float64],
    second: NDArray[np.float64],
    comparisons_a: NDArray[np.int_],
    comparisons_b: NDArray[np.int_],
) -> float:
    distances_first = squared_distance_matrix(first)
    distances_second = squared_distance_matrix(second)
    first_delta = distances_first[comparisons_a[:, 0], comparisons_a[:, 1]] - (
        distances_first[comparisons_b[:, 0], comparisons_b[:, 1]]
    )
    second_delta = distances_second[comparisons_a[:, 0], comparisons_a[:, 1]] - (
        distances_second[comparisons_b[:, 0], comparisons_b[:, 1]]
    )
    first_sign = np.sign(first_delta)
    second_sign = np.sign(second_delta)
    comparable = first_sign != 0.0
    if not np.any(comparable):
        return float("nan")
    return float(np.mean(first_sign[comparable] != second_sign[comparable]))


def pairwise_order_stability(
    embeddings: list[ArrayLike],
    num_pair_comparisons: int = 5000,
    seed: int | None = None,
) -> dict[str, float]:
    """Compare distance-order disagreement across fitted embeddings."""

    if len(embeddings) < 2:
        return {
            "mean_order_disagreement": float("nan"),
            "std_order_disagreement": float("nan"),
            "max_order_disagreement": float("nan"),
            "pair_count": 0.0,
        }
    arrays = [np.asarray(embedding, dtype=float) for embedding in embeddings]
    n_points = arrays[0].shape[0]
    if any(array.ndim != 2 or array.shape[0] != n_points for array in arrays):
        raise ValueError("all embeddings must have shape (n_points, embedding_dim)")
    count = int(num_pair_comparisons)
    if count <= 0:
        raise ValueError("num_pair_comparisons must be positive")
    rng = np.random.default_rng(seed)
    first_pairs = rng.integers(0, n_points, size=(count, 2))
    second_pairs = rng.integers(0, n_points, size=(count, 2))
    valid = (first_pairs[:, 0] != first_pairs[:, 1]) & (
        second_pairs[:, 0] != second_pairs[:, 1]
    )
    first_pairs = first_pairs[valid]
    second_pairs = second_pairs[valid]
    values = [
        _order_disagreement(a, b, first_pairs, second_pairs)
        for a, b in itertools.combinations(arrays, 2)
    ]
    array = np.asarray(values, dtype=float)
    return {
        "mean_order_disagreement": float(np.nanmean(array)),
        "std_order_disagreement": float(np.nanstd(array)),
        "max_order_disagreement": float(np.nanmax(array)),
        "pair_count": float(array.size),
    }
