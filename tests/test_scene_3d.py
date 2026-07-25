"""Regressions for the 3+1D scene builder (P8).

The builder is the only 3+1D-specific piece of the arm; everything
downstream is the frozen PC-V1 pipeline. So these pin the preconditions
the preregistration relies on -- that the observer shell spans, that
chains stay inside the diamond, that the target band is what it claims,
and that a seed reproduces a scene exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.positive_control.scene import SceneValidityError
from causal_spacetime_lab.positive_control.scene_3d import (
    Scene3DConfig,
    affine_rank,
    build_scene_3plus1d,
    fibonacci_shell,
    target_positions_3d,
)


def test_the_observer_shell_spans_three_dimensions():
    """A coplanar observer set is exactly the degenerate case for
    multilateration, so spanning is checked rather than assumed."""

    for count in (4, 6, 8, 12, 16, 20):
        shell = fibonacci_shell(count, 0.25)
        assert shell.shape == (count, 3)
        assert affine_rank(shell) == 3, count
        np.testing.assert_allclose(np.linalg.norm(shell, axis=1), 0.25)


def test_the_shell_refuses_a_count_that_cannot_span():
    with pytest.raises(ValueError, match="at least 4"):
        fibonacci_shell(3, 0.25)


def test_the_shell_is_deterministic_without_a_seed():
    """A frozen configuration should be reproducible from the count
    alone, with no RNG state in the way."""

    np.testing.assert_array_equal(
        fibonacci_shell(12, 0.25), fibonacci_shell(12, 0.25)
    )


@pytest.mark.slow
def test_a_default_scene_builds_and_respects_its_own_band():
    config = Scene3DConfig(seed=100)
    scene = build_scene_3plus1d(config)

    assert config.min_targets <= scene.target_indices.size <= config.max_targets
    assert len(scene.chain_index_arrays) == config.chain_count

    targets = scene.events[scene.target_indices]
    assert np.all(np.abs(targets[:, 0]) <= config.target_band_t)
    radius = np.linalg.norm(targets[:, 1:4], axis=1)
    assert np.all(radius <= config.target_band_radius + 1e-12)

    positions = target_positions_3d(scene)
    assert positions.shape == (scene.target_indices.size, 3)
    np.testing.assert_array_equal(positions, targets[:, 1:4])


@pytest.mark.slow
def test_every_chain_stays_inside_the_diamond():
    config = Scene3DConfig(seed=100)
    scene = build_scene_3plus1d(config)
    half_t = config.diamond_T / 2.0
    for chain in scene.chain_index_arrays:
        block = scene.events[chain]
        radius = np.linalg.norm(block[:, 1:4], axis=1)
        assert np.allclose(radius, config.chain_shell_radius)
        # the diamond condition itself: ||x|| <= T/2 - |t|
        assert np.all(radius <= half_t - np.abs(block[:, 0]) + 1e-12)


@pytest.mark.slow
def test_the_same_seed_reproduces_the_same_scene():
    first = build_scene_3plus1d(Scene3DConfig(seed=100))
    second = build_scene_3plus1d(Scene3DConfig(seed=100))
    assert first.events_digest == second.events_digest
    assert first.causal_digest == second.causal_digest
    np.testing.assert_array_equal(first.target_indices, second.target_indices)


@pytest.mark.slow
def test_a_shell_that_leaves_the_diamond_is_refused_rather_than_clipped():
    """The chains span t in [-0.7, 0.7] inside a T = 2 diamond, so any
    shell radius above 0.3 puts a tick outside. That has to raise, not
    silently produce a scene whose observers are not there."""

    with pytest.raises(SceneValidityError, match="leaves the diamond"):
        build_scene_3plus1d(Scene3DConfig(seed=100, chain_shell_radius=0.45))


@pytest.mark.slow
def test_too_few_events_is_reported_as_an_invalid_scene():
    """P8-B records invalid seeds and excludes them; it never tops the
    denominator up. So the failure has to be typed, not a crash."""

    with pytest.raises(SceneValidityError, match="eligible targets"):
        build_scene_3plus1d(Scene3DConfig(seed=100, n_events=400))
