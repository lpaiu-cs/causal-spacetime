"""Exact-case tests for bracket-width echo profile measurement."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer import (
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
)
from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
    relabel_targets,
    shuffle_profile_columns,
)
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlScene,
    PositiveControlSceneConfig,
)


def _hand_scene(targets: np.ndarray) -> PositiveControlScene:
    """One stationary chain at x=0 with ticks at -0.8, -0.7, ..., 0.8."""

    chain_events, _ = make_stationary_observer_chain_1p1(1.6, 17, x=0.0)
    events = np.vstack([targets, chain_events])
    return PositiveControlScene(
        config=PositiveControlSceneConfig(seed=0),
        events=events,
        causal=causal_matrix_1p1(events),
        bulk_count=targets.shape[0],
        chain_index_arrays=(observer_chain_indices(targets.shape[0], 17),),
        tick_ranks=np.arange(17, dtype=np.float64),
        target_indices=np.arange(targets.shape[0], dtype=int),
        events_digest="test",
        causal_digest="test",
    )


def test_bracket_widths_match_hand_computation() -> None:
    # Tick index i sits at time -0.8 + 0.1 * i on the x = 0 worldline.
    targets = np.array([[0.0, 0.1], [0.0, 0.2], [0.05, 0.1]])
    profiles = measure_bracket_echo_profiles(_hand_scene(targets))

    # (0.0, 0.1): tau- = -0.1 (index 7, null-related), tau+ = 0.1 (index 9).
    # (0.0, 0.2): tau- = -0.2 (index 6), tau+ = 0.2 (index 10).
    # (0.05, 0.1): last predecessor index 7 (-0.1), first successor index 10.
    expected = np.array([[2.0], [4.0], [3.0]])
    assert np.array_equal(profiles.delay_ranks, expected)
    assert bool(np.all(profiles.reachable))


def test_bracket_width_is_monotone_in_distance() -> None:
    chain_events, _ = make_stationary_observer_chain_1p1(1.8, 181, x=0.0)
    offsets = np.array([0.05, 0.15, 0.25, 0.35])
    targets = np.column_stack([np.zeros_like(offsets), offsets])
    events = np.vstack([targets, chain_events])
    scene = PositiveControlScene(
        config=PositiveControlSceneConfig(seed=0),
        events=events,
        causal=causal_matrix_1p1(events),
        bulk_count=targets.shape[0],
        chain_index_arrays=(observer_chain_indices(targets.shape[0], 181),),
        tick_ranks=np.arange(181, dtype=np.float64),
        target_indices=np.arange(targets.shape[0], dtype=int),
        events_digest="test",
        causal_digest="test",
    )
    widths = measure_bracket_echo_profiles(scene).delay_ranks[:, 0]
    assert bool(np.all(np.diff(widths) > 0))


def test_out_of_bracket_target_is_missing_not_zero() -> None:
    targets = np.array([[0.85, 0.0]])  # strictly after the last tick at 0.8
    profiles = measure_bracket_echo_profiles(_hand_scene(targets))
    assert not bool(profiles.reachable[0, 0])
    assert bool(np.isnan(profiles.delay_ranks[0, 0]))


def test_column_shuffle_preserves_marginals_and_breaks_rows() -> None:
    rng = np.random.default_rng(0)
    targets = np.column_stack(
        [rng.uniform(-0.05, 0.05, 12), rng.uniform(-0.3, 0.3, 12)]
    )
    profiles = measure_bracket_echo_profiles(_hand_scene(targets))
    shuffled = shuffle_profile_columns(profiles, seed=1)
    original = np.sort(profiles.delay_ranks[:, 0])
    permuted = np.sort(shuffled.delay_ranks[:, 0])
    assert np.array_equal(original, permuted)


def test_relabeling_permutes_profiles_and_truth_together() -> None:
    targets = np.array([[0.0, 0.05], [0.0, 0.15], [0.0, 0.25], [0.0, 0.35]])
    profiles = measure_bracket_echo_profiles(_hand_scene(targets))
    truth = targets[:, 1]
    relabeled, permuted_truth = relabel_targets(profiles, truth, seed=5)
    order = np.argsort(permuted_truth)
    assert np.array_equal(
        relabeled.delay_ranks[order],
        profiles.delay_ranks[np.argsort(truth)],
    )
