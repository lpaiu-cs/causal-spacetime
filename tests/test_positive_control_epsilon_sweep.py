"""Tests for the P1 epsilon-diluted geometry order generator."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
)
from causal_spacetime_lab.positive_control.dissimilarity import (
    build_constraint_split,
    margin_from_probe_quantile,
    profile_dissimilarity_matrix,
)
from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
)
from causal_spacetime_lab.positive_control.epsilon_sweep import (
    build_epsilon_scene,
    epsilon_diluted_order,
    transitive_reduction,
)
from causal_spacetime_lab.positive_control.rewire import (
    geometric_post_closure_density,
    transitive_closure,
)
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlSceneConfig,
    build_positive_control_scene,
)

SMOKE = PositiveControlSceneConfig(
    n_events=400,
    chain_positions=(-0.2, -0.07, 0.07, 0.2),
    ticks_per_chain=32,
    target_band_t=0.10,
    target_band_x=0.20,
    max_targets=16,
    min_targets=8,
    min_bracketing_chains=4,
    seed=3,
)


def test_transitive_reduction_recloses_to_original() -> None:
    scene = build_positive_control_scene(SMOKE)
    reduction = transitive_reduction(scene.causal)
    # the reduction has no self loops and re-closes to the original order
    assert not np.any(np.diag(reduction))
    assert np.array_equal(transitive_closure(reduction), scene.causal)
    # a covering relation is a strict subset of the closed relation
    assert np.all(scene.causal[reduction])
    assert int(reduction.sum()) < int(scene.causal.sum())


def test_epsilon_zero_returns_pristine_scene() -> None:
    scene = build_positive_control_scene(SMOKE)
    same, density = build_epsilon_scene(scene, 0.0, seed=41, target_density=0.5)
    assert same is scene
    assert density == 0.5


def test_epsilon_order_holds_density_and_stays_transitive() -> None:
    scene = build_positive_control_scene(SMOKE)
    target, _, _ = geometric_post_closure_density(scene)
    for epsilon in (0.25, 0.5, 1.0):
        closed, density, _ = epsilon_diluted_order(scene, epsilon, 41, target)
        assert np.array_equal(closed, transitive_closure(closed))  # transitive
        assert not np.any(np.diag(closed))  # irreflexive
        assert abs(density - target) <= 0.05  # density held
        for chain in scene.chain_index_arrays:  # chains preserved
            assert bool(np.all(closed[chain[:-1], chain[1:]]))


def _truth_error_at(scene, seed) -> float:
    profiles = measure_bracket_echo_profiles(scene)
    matrix = profile_dissimilarity_matrix(profiles, 4)
    margin = margin_from_probe_quantile(matrix, seed=seed + 3)
    split = build_constraint_split(matrix, 1500, 400, margin, seed=seed + 5)
    coords, _ = fit_ordinal_embedding_gradient_descent(
        profiles.target_count, 1, split.train, steps=800, learning_rate=0.05,
        seed=seed + 100, restarts=3,
    )
    return embedding_distance_order_error(
        coords, scene.target_coords[:, 1].reshape(-1, 1),
        num_pair_comparisons=6000, seed=seed + 9,
    )


def test_geometry_recovery_degrades_from_eps0_to_eps1() -> None:
    scene = build_positive_control_scene(SMOKE)
    target, _, _ = geometric_post_closure_density(scene)
    err0 = _truth_error_at(scene, seed=3)
    eps1_scene, _ = build_epsilon_scene(scene, 1.0, seed=44, target_density=target)
    err1 = _truth_error_at(eps1_scene, seed=3)
    # geometric recovers true x-order; fully diluted does not
    assert err0 < 0.20
    assert err1 > err0 + 0.15
