"""Geometry-free control orders for Stage C (PC-V1 control C-RANDOM-ORDER).

The epsilon = 1 endpoint: a time-respecting random partial order with
Bernoulli edge density matched to the geometric scene's bulk relation
density before transitive closure. Chain-internal tick order is preserved so
the reference protocols remain valid chains. The graded epsilon sweep is
deferred to the P1 preregistration.
"""

from __future__ import annotations

import dataclasses
import hashlib

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.positive_control.scene import (
    PositiveControlScene,
    SceneValidityError,
    two_sided_bracket_counts,
)


def transitive_closure(relation: NDArray[np.bool_]) -> NDArray[np.bool_]:
    """Return the transitive closure of an acyclic boolean relation.

    Reachability doubling: closure |= closure @ closure until a fixpoint. The
    boolean matmul is done in float32 (BLAS-backed); sums are bounded by the
    node count, well below float32's exact-integer range, so ``> 0`` is exact.
    """

    closure = np.array(relation, dtype=bool)
    while True:
        reach = (
            closure.astype(np.float32) @ closure.astype(np.float32)
        ) > 0.0
        expanded = closure | reach
        if bool(np.array_equal(expanded, closure)):
            return closure
        closure = expanded


def bulk_relation_density(scene: PositiveControlScene) -> float:
    """Relation density among time-ordered bulk event pairs."""

    n_bulk = scene.bulk_count
    times = scene.events[:n_bulk, 0]
    ordered = times[:, None] < times[None, :]
    related = scene.causal[:n_bulk, :n_bulk] & ordered
    ordered_count = int(np.sum(ordered))
    if ordered_count == 0:
        raise ValueError("no time-ordered bulk pairs")
    return float(np.sum(related) / ordered_count)


def geometric_post_closure_density(
    scene: PositiveControlScene,
) -> tuple[float, NDArray[np.bool_], int]:
    """Return the geometric causal density among time-ordered event pairs.

    This is the density the geometry-free control is tuned to match
    (PC-V1 Section 9, amended). Also returns the time-order mask and its
    pair count so callers avoid recomputing them.
    """

    times = scene.events[:, 0]
    ordered = times[:, None] < times[None, :]
    ordered_count = int(np.sum(ordered))
    if ordered_count == 0:
        raise ValueError("no time-ordered event pairs")
    density = float(np.sum(scene.causal & ordered) / ordered_count)
    return density, ordered, ordered_count


def _closed_order_at_q(
    ordered: NDArray[np.bool_],
    uniform: NDArray[np.float64],
    q: float,
    chain_index_arrays: tuple[NDArray[np.int_], ...],
) -> NDArray[np.bool_]:
    relation = ordered & (uniform < q)
    for chain in chain_index_arrays:
        relation[chain[:-1], chain[1:]] = True
    return transitive_closure(relation)


def random_time_respecting_order(
    scene: PositiveControlScene,
    seed: int,
    target_density: float | None = None,
    tolerance: float = 0.03,
    max_iters: int = 14,
) -> tuple[NDArray[np.bool_], float, float]:
    """Density-matched geometry-free order (PC-V1 Section 9, amended).

    A time-respecting transitive percolation whose pre-closure edge
    probability ``q`` is tuned by bisection so the achieved post-closure
    relation density among time-ordered pairs matches the geometric scene's
    causal density (chain-internal edges preserved). A single shared uniform
    draw makes density monotone in ``q``, so bisection is exact and
    deterministic. Returns the closed relation, achieved density, and ``q``.
    """

    geometric_density, ordered, ordered_count = geometric_post_closure_density(
        scene
    )
    target = geometric_density if target_density is None else float(target_density)
    rng = np.random.default_rng(seed)
    n = scene.events.shape[0]
    uniform = rng.uniform(size=(n, n))

    low, high = 0.0, 1.0
    best_closed: NDArray[np.bool_] | None = None
    best_density = float("inf")
    best_q = float("nan")
    for _ in range(max_iters):
        mid = 0.5 * (low + high)
        closed = _closed_order_at_q(
            ordered, uniform, mid, scene.chain_index_arrays
        )
        density = float(np.sum(closed & ordered) / ordered_count)
        if abs(density - target) < abs(best_density - target):
            best_closed, best_density, best_q = closed, density, mid
        if abs(density - target) <= tolerance:
            break
        if density < target:
            low = mid
        else:
            high = mid

    if best_closed is None:
        raise RuntimeError("bisection produced no candidate order")
    np.fill_diagonal(best_closed, False)
    return best_closed, best_density, best_q


def geometry_free_scene(
    scene: PositiveControlScene,
    seed: int,
) -> tuple[PositiveControlScene, float]:
    """Return a scene with the causal order replaced by a random order.

    Event coordinates are retained only as labels; ground-truth spatial
    metrics are undefined for this scene. The preregistered target-selection
    policy is re-applied under the random order.
    """

    closed, achieved_density, _ = random_time_respecting_order(scene, seed)
    config = scene.config

    bulk = scene.events[: scene.bulk_count]
    in_band = (np.abs(bulk[:, 0]) <= config.target_band_t) & (
        np.abs(bulk[:, 1]) <= config.target_band_x
    )
    candidates = np.flatnonzero(in_band)
    bracket_counts = two_sided_bracket_counts(
        closed,
        scene.chain_index_arrays,
        scene.tick_ranks,
        candidates,
    )
    eligible = candidates[bracket_counts >= config.min_bracketing_chains]
    if eligible.size < config.min_targets:
        raise SceneValidityError(
            f"random-order control: only {eligible.size} eligible targets"
        )
    if eligible.size > config.max_targets:
        subsample_rng = np.random.default_rng(config.seed + 1)
        chosen = subsample_rng.choice(
            eligible.size,
            size=config.max_targets,
            replace=False,
        )
        eligible = np.sort(eligible[chosen])

    causal_digest = hashlib.sha256(
        np.ascontiguousarray(closed).tobytes()
    ).hexdigest()[:16]
    control_scene = dataclasses.replace(
        scene,
        causal=closed,
        target_indices=eligible.astype(int),
        causal_digest=causal_digest,
    )
    return control_scene, achieved_density
