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
    """Return the transitive closure of an acyclic boolean relation."""

    closure = np.asarray(relation, dtype=bool).copy()
    while True:
        expanded = closure | (
            np.matmul(closure.astype(np.uint8), closure.astype(np.uint8)) > 0
        )
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


def random_time_respecting_order(
    scene: PositiveControlScene,
    seed: int,
) -> tuple[NDArray[np.bool_], float]:
    """Random time-respecting order with matched pre-closure edge density.

    Returns the transitively closed relation and the achieved post-closure
    density among time-ordered pairs (recorded, per preregistration).
    """

    rng = np.random.default_rng(seed)
    density = bulk_relation_density(scene)
    times = scene.events[:, 0]
    n = times.size
    ordered = times[:, None] < times[None, :]
    relation = ordered & (rng.uniform(size=(n, n)) < density)

    for chain in scene.chain_index_arrays:
        for offset in range(chain.size - 1):
            relation[chain[offset], chain[offset + 1]] = True

    closed = transitive_closure(relation)
    np.fill_diagonal(closed, False)
    achieved = float(np.sum(closed & ordered) / max(1, int(np.sum(ordered))))
    return closed, achieved


def geometry_free_scene(
    scene: PositiveControlScene,
    seed: int,
) -> tuple[PositiveControlScene, float]:
    """Return a scene with the causal order replaced by a random order.

    Event coordinates are retained only as labels; ground-truth spatial
    metrics are undefined for this scene. The preregistered target-selection
    policy is re-applied under the random order.
    """

    closed, achieved_density = random_time_respecting_order(scene, seed)
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
