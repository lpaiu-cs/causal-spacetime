"""Sampling utilities for validation experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")
    return array


def _as_square_bool_matrix(matrix: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(matrix, dtype=bool)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("causal_matrix must be a square boolean matrix")
    return array


def _interval_counts_for_candidates(
    causal_matrix: NDArray[np.bool_],
    sources: NDArray[np.int_],
    targets: NDArray[np.int_],
    chunk_size: int = 20_000,
) -> NDArray[np.int64]:
    counts = np.empty(sources.shape[0], dtype=np.int64)
    for start in range(0, sources.shape[0], chunk_size):
        stop = min(start + chunk_size, sources.shape[0])
        source_chunk = sources[start:stop]
        target_chunk = targets[start:stop]
        counts[start:stop] = np.count_nonzero(
            causal_matrix[source_chunk] & causal_matrix[:, target_chunk].T,
            axis=1,
        )
    return counts


def sample_timelike_pairs(
    events: ArrayLike,
    causal_matrix: ArrayLike,
    num_pairs: int,
    seed: int | None = None,
    min_tau: float | None = None,
    min_interval_count: int | None = None,
) -> list[tuple[int, int]]:
    """Sample pairs ``(i, j)`` where event ``i`` causally precedes event ``j``.

    Causality is determined only from ``causal_matrix``. If ``min_tau`` is
    supplied, true coordinates are used only for that explicit optional filter,
    which is useful for avoiding near-null or tiny intervals in validation
    studies where hidden Minkowski coordinates are available.
    """

    event_array = _as_event_array(events)
    causal_order = _as_square_bool_matrix(causal_matrix)
    if event_array.shape[0] != causal_order.shape[0]:
        raise ValueError("events and causal_matrix must contain the same count")

    requested = int(num_pairs)
    if requested < 0:
        raise ValueError("num_pairs must be non-negative")
    if requested == 0:
        return []

    sources, targets = np.nonzero(causal_order)

    if min_tau is not None:
        threshold = float(min_tau)
        if threshold < 0:
            raise ValueError("min_tau must be non-negative")
        dt = event_array[targets, 0] - event_array[sources, 0]
        dx = event_array[targets, 1] - event_array[sources, 1]
        tau_squared = np.maximum(dt * dt - dx * dx, 0.0)
        keep = tau_squared >= threshold * threshold
        sources = sources[keep]
        targets = targets[keep]

    if min_interval_count is not None:
        min_count = int(min_interval_count)
        if min_count < 0:
            raise ValueError("min_interval_count must be non-negative")
        counts = _interval_counts_for_candidates(causal_order, sources, targets)
        keep = counts >= min_count
        sources = sources[keep]
        targets = targets[keep]

    candidate_count = sources.shape[0]
    if candidate_count == 0:
        return []

    rng = np.random.default_rng(seed)
    sample_size = min(requested, candidate_count)
    chosen = rng.choice(candidate_count, size=sample_size, replace=False)
    return [(int(sources[index]), int(targets[index])) for index in chosen]

