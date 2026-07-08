"""Epsilon-diluted geometry orders for the P1 dose-response sweep.

Geometry is the manipulated variable: epsilon = 0 is the measured Minkowski
causal order; epsilon = 1 is a density-matched geometry-free order; in between,
a fraction epsilon of the geometric COVERING edges is rewired to random
time-respecting edges and the relation is re-closed, holding the post-closure
relation density fixed (the D2 fairness lesson) so the sweep isolates geometry
from density. See docs/prereg/p1_epsilon_sweep.md.
"""

from __future__ import annotations

import dataclasses
import hashlib

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.positive_control.rewire import transitive_closure
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlScene,
    SceneValidityError,
    two_sided_bracket_counts,
)


def transitive_reduction(closed: NDArray[np.bool_]) -> NDArray[np.bool_]:
    """Return the covering relation (Hasse edges) of a closed acyclic order.

    An edge ``i -> j`` is covering iff no intermediate ``k`` has
    ``i -> k -> j``. Computed as ``closed & ~(closed @ closed)`` over the
    strict relation (float32 BLAS matmul; ``> 0`` is exact below the node
    count). Removing a covering edge genuinely reduces reachability, so
    rewiring covering edges dilutes geometry effectively.
    """

    strict = np.array(closed, dtype=bool)
    np.fill_diagonal(strict, False)
    two_hop = (strict.astype(np.float32) @ strict.astype(np.float32)) > 0.0
    return strict & ~two_hop


def _chain_edge_mask(scene: PositiveControlScene) -> NDArray[np.bool_]:
    mask = np.zeros_like(scene.causal, dtype=bool)
    for chain in scene.chain_index_arrays:
        mask[chain[:-1], chain[1:]] = True
    return mask


def epsilon_diluted_order(
    scene: PositiveControlScene,
    epsilon: float,
    seed: int,
    target_density: float,
    tolerance: float = 0.02,
    max_iters: int = 16,
) -> tuple[NDArray[np.bool_], float, float]:
    """Return a density-held order that dilutes geometry by fraction epsilon.

    Keep each non-chain geometric covering edge with probability
    ``1 - epsilon``, add random time-respecting edges at a bisection-tuned
    probability so the achieved post-closure density matches
    ``target_density``, always keep chain-internal covering edges, then
    transitively close. Returns the closed relation, achieved density, and the
    tuned random-edge probability.
    """

    if not 0.0 <= epsilon <= 1.0:
        raise ValueError("epsilon must be in [0, 1]")
    times = scene.events[:, 0]
    ordered = times[:, None] < times[None, :]
    ordered_count = int(np.sum(ordered))
    reduction = transitive_reduction(scene.causal)
    chains = _chain_edge_mask(scene)
    geometric_links = reduction & ~chains

    rng = np.random.default_rng(seed)
    kept = geometric_links & (rng.uniform(size=scene.causal.shape) < 1.0 - epsilon)
    uniform = rng.uniform(size=scene.causal.shape)

    def build(probability: float) -> NDArray[np.bool_]:
        pre = kept | (ordered & (uniform < probability)) | chains
        return transitive_closure(pre)

    low, high = 0.0, 0.15
    best: NDArray[np.bool_] | None = None
    best_density = float("inf")
    best_probability = 0.0
    for _ in range(max_iters):
        mid = 0.5 * (low + high)
        closed = build(mid)
        density = float(np.sum(closed & ordered) / ordered_count)
        if abs(density - target_density) < abs(best_density - target_density):
            best, best_density, best_probability = closed, density, mid
        if abs(density - target_density) <= tolerance:
            break
        if density < target_density:
            low = mid
        else:
            high = mid

    if best is None:
        raise RuntimeError("bisection produced no candidate order")
    np.fill_diagonal(best, False)
    return best, best_density, best_probability


def build_epsilon_scene(
    scene: PositiveControlScene,
    epsilon: float,
    seed: int,
    target_density: float,
) -> tuple[PositiveControlScene, float]:
    """Return a scene whose order is geometry diluted by fraction epsilon.

    Epsilon = 0 returns the pristine geometric scene unchanged (continuity
    with PC-V1). For epsilon > 0 the target set is re-selected under the
    diluted order using the same preregistered policy; too few eligible
    targets raises ``SceneValidityError`` (recorded, never silently dropped).
    True event coordinates are retained as validation labels.
    """

    if epsilon == 0.0:
        return scene, target_density

    closed, achieved_density, _ = epsilon_diluted_order(
        scene, epsilon, seed, target_density
    )
    config = scene.config
    bulk = scene.events[: scene.bulk_count]
    in_band = (np.abs(bulk[:, 0]) <= config.target_band_t) & (
        np.abs(bulk[:, 1]) <= config.target_band_x
    )
    candidates = np.flatnonzero(in_band)
    bracket_counts = two_sided_bracket_counts(
        closed, scene.chain_index_arrays, scene.tick_ranks, candidates
    )
    eligible = candidates[bracket_counts >= config.min_bracketing_chains]
    if eligible.size < config.min_targets:
        raise SceneValidityError(
            f"epsilon={epsilon}: only {eligible.size} eligible targets"
        )
    if eligible.size > config.max_targets:
        subsample_rng = np.random.default_rng(config.seed + 1)
        chosen = subsample_rng.choice(
            eligible.size, size=config.max_targets, replace=False
        )
        eligible = np.sort(eligible[chosen])

    causal_digest = hashlib.sha256(
        np.ascontiguousarray(closed).tobytes()
    ).hexdigest()[:16]
    epsilon_scene = dataclasses.replace(
        scene,
        causal=closed,
        target_indices=eligible.astype(int),
        causal_digest=causal_digest,
    )
    return epsilon_scene, achieved_density
