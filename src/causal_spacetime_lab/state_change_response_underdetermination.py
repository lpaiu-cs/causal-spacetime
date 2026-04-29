"""Underdetermination diagnostics for scalar response orders.

The helpers in this module construct rank-preserving scalar layouts as
representation witnesses for a response order. They are not metric
coordinates and do not define physical target-target distances.
"""

from __future__ import annotations

import itertools

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_response_signature import (
    response_order_sign_matrix,
)


def _as_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    vector = np.asarray(values, dtype=float)
    if vector.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return vector


def _reachable_delay_data(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.bool_], NDArray[np.float64]]:
    delays = _as_vector(delay_ranks, "delay_ranks")
    reachable = np.asarray(reachable_mask, dtype=bool)
    if reachable.ndim != 1 or reachable.shape != delays.shape:
        raise ValueError("reachable_mask must match delay_ranks")
    unique = np.unique(delays[reachable])
    return delays, reachable, unique


def rank_preserving_scalar_values(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    *,
    spacing: str = "linear",
    seed: int | None = None,
) -> NDArray[np.float64]:
    """Return scalar values preserving the reachable response-rank order.

    The returned scalars are order witnesses only. They are not metric
    coordinates, calibrated time values, or target-target distances.
    """

    delays, reachable, unique = _reachable_delay_data(delay_ranks, reachable_mask)
    values = np.full(delays.shape, np.nan, dtype=float)
    if unique.size == 0:
        return values

    if spacing == "linear":
        coordinates = np.arange(unique.size, dtype=float)
    elif spacing == "exponential":
        coordinates = np.power(2.0, np.arange(unique.size, dtype=float)) - 1.0
    elif spacing == "random_monotone":
        rng = np.random.default_rng(seed)
        coordinates = np.cumsum(rng.uniform(0.25, 3.0, size=unique.size))
    else:
        raise ValueError(
            "spacing must be one of: linear, exponential, random_monotone"
        )

    by_rank = {
        float(rank): float(coord)
        for rank, coord in zip(unique, coordinates, strict=True)
    }
    for index, delay in enumerate(delays):
        if reachable[index]:
            values[index] = by_rank[float(delay)]
    return values


def response_order_preserved(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    scalar_values: ArrayLike,
) -> bool:
    """Return whether scalar values preserve the response-order sign matrix."""

    delays, reachable, _ = _reachable_delay_data(delay_ranks, reachable_mask)
    scalars = _as_vector(scalar_values, "scalar_values")
    if scalars.shape != delays.shape:
        raise ValueError("scalar_values must match delay_ranks")
    if np.any(~np.isfinite(scalars[reachable])):
        return False

    delay_signs = response_order_sign_matrix(delays.astype(int), reachable)
    scalar_signs = np.zeros(delay_signs.shape, dtype=int)
    reachable_pairs = reachable[:, None] & reachable[None, :]
    scalar_signs[reachable_pairs & (scalars[:, None] < scalars[None, :])] = -1
    scalar_signs[reachable_pairs & (scalars[:, None] > scalars[None, :])] = 1
    return bool(np.array_equal(delay_signs, scalar_signs))


def generate_rank_preserving_1d_layouts(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    layout_count: int = 10,
    seed: int | None = None,
) -> list[NDArray[np.float64]]:
    """Generate distinct 1D layouts preserving one scalar response order.

    The layouts can induce different synthetic target-target distance orders
    while preserving the same single-reference response order.
    """

    if layout_count < 1:
        raise ValueError("layout_count must be at least 1")
    delays, reachable, unique = _reachable_delay_data(delay_ranks, reachable_mask)
    rng = np.random.default_rng(seed)
    layouts: list[NDArray[np.float64]] = [
        rank_preserving_scalar_values(delays, reachable, spacing="linear"),
        rank_preserving_scalar_values(delays, reachable, spacing="exponential"),
    ]
    if unique.size:
        cluster = np.full(delays.shape, np.nan, dtype=float)
        cluster_values = np.cumsum(np.linspace(0.1, 1.0, unique.size))
        for rank, value in zip(unique, cluster_values, strict=True):
            cluster[reachable & (delays == rank)] = float(value)
        layouts.append(cluster)
    while len(layouts) < layout_count:
        layouts.append(
            rank_preserving_scalar_values(
                delays,
                reachable,
                spacing="random_monotone",
                seed=int(rng.integers(0, 2**32 - 1)),
            )
        )
    return layouts[:layout_count]


def _target_pairs(n_targets: int) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(max(0, n_targets - 1))
        for j in range(i + 1, n_targets)
    ]


def pairwise_distance_order_signature_from_layout(
    layout: ArrayLike,
    reachable_mask: ArrayLike,
) -> NDArray[np.int_]:
    """Return pair-pair distance-order signs induced by a synthetic layout.

    This is only an underdetermination diagnostic. The layout is not a
    physical space.
    """

    points = _as_vector(layout, "layout")
    reachable = np.asarray(reachable_mask, dtype=bool)
    if reachable.ndim != 1 or reachable.shape != points.shape:
        raise ValueError("reachable_mask must match layout")
    pairs = _target_pairs(points.size)
    distances = np.full(len(pairs), np.nan, dtype=float)
    for index, (i, j) in enumerate(pairs):
        if reachable[i] and reachable[j] and np.isfinite(points[i]) and np.isfinite(
            points[j]
        ):
            distances[index] = abs(points[i] - points[j])
    signs = np.zeros((len(pairs), len(pairs)), dtype=int)
    valid = np.isfinite(distances)
    valid_pairs = valid[:, None] & valid[None, :]
    signs[valid_pairs & (distances[:, None] < distances[None, :])] = -1
    signs[valid_pairs & (distances[:, None] > distances[None, :])] = 1
    return signs


def pairwise_distance_order_disagreement(
    layout_a: ArrayLike,
    layout_b: ArrayLike,
    reachable_mask: ArrayLike,
) -> float:
    """Compare synthetic pair-distance orders from two rank-preserving layouts."""

    signs_a = pairwise_distance_order_signature_from_layout(layout_a, reachable_mask)
    signs_b = pairwise_distance_order_signature_from_layout(layout_b, reachable_mask)
    if signs_a.shape != signs_b.shape:
        raise ValueError("pairwise signatures must have the same shape")
    pair_indices = [
        (i, j)
        for i, j in itertools.combinations(range(signs_a.shape[0]), 2)
        if signs_a[i, j] != 0 or signs_b[i, j] != 0
    ]
    if not pair_indices:
        return 0.0
    disagreements = [signs_a[i, j] != signs_b[i, j] for i, j in pair_indices]
    return float(np.mean(disagreements))


def single_reference_underdetermination_report(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    layout_count: int = 20,
    seed: int | None = None,
) -> dict[str, float]:
    """Report target-target order variation under one response order."""

    delays, reachable, unique = _reachable_delay_data(delay_ranks, reachable_mask)
    layouts = generate_rank_preserving_1d_layouts(
        delays,
        reachable,
        layout_count=layout_count,
        seed=seed,
    )
    preserved = [
        response_order_preserved(delays, reachable, layout) for layout in layouts
    ]
    compared_layouts = layouts[: min(6, len(layouts))]
    disagreements = [
        pairwise_distance_order_disagreement(a, b, reachable)
        for a, b in zip(compared_layouts[:-1], compared_layouts[1:], strict=True)
    ]
    return {
        "reachable_count": float(np.sum(reachable)),
        "unique_delay_rank_count": float(unique.size),
        "generated_layout_count": float(len(layouts)),
        "response_order_preserved_fraction": float(np.mean(preserved))
        if preserved
        else float("nan"),
        "mean_pair_distance_order_disagreement": float(np.mean(disagreements))
        if disagreements
        else 0.0,
        "max_pair_distance_order_disagreement": float(np.max(disagreements))
        if disagreements
        else 0.0,
    }
