"""Tests for the P2 2+1D scene builder."""

from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
)
from causal_spacetime_lab.positive_control.scene import SceneValidityError
from causal_spacetime_lab.positive_control.scene_2d import (
    Scene2DConfig,
    build_scene_2plus1d,
    target_positions_2d,
)

SMALL = Scene2DConfig(
    n_events=1200,
    chain_count=8,
    chain_ring_radius=0.25,
    ticks_per_chain=48,
    target_band_t=0.10,
    target_band_radius=0.22,
    max_targets=30,
    min_targets=8,
    min_bracketing_chains=8,
    seed=1,
)


def test_scene_is_2plus1d_and_content_addressed() -> None:
    scene = build_scene_2plus1d(SMALL)
    assert scene.events.shape[1] == 3  # (t, x, y)
    assert scene.causal.shape[0] == scene.events.shape[0]
    again = build_scene_2plus1d(SMALL)
    assert scene.events_digest == again.events_digest
    assert scene.causal_digest == again.causal_digest


def test_chains_form_causal_chains_on_a_ring() -> None:
    scene = build_scene_2plus1d(SMALL)
    assert len(scene.chain_index_arrays) == SMALL.chain_count
    for chain in scene.chain_index_arrays:
        for earlier in range(chain.size - 1):
            assert scene.causal[chain[earlier], chain[earlier + 1]]
    # chain spatial positions lie on the ring
    firsts = np.array([scene.events[chain[0]] for chain in scene.chain_index_arrays])
    radii = np.hypot(firsts[:, 1], firsts[:, 2])
    assert np.allclose(radii, SMALL.chain_ring_radius)


def test_targets_in_2d_band_and_bracketed() -> None:
    scene = build_scene_2plus1d(SMALL)
    xy = target_positions_2d(scene)
    assert xy.shape[1] == 2
    coords = scene.target_coords
    assert np.all(np.abs(coords[:, 0]) <= SMALL.target_band_t)
    radii = np.hypot(coords[:, 1], coords[:, 2])
    assert np.all(radii <= SMALL.target_band_radius + 1e-9)
    assert scene.target_indices.size >= SMALL.min_targets
    profiles = measure_bracket_echo_profiles(scene)
    # every target reached by every chain (two-sided bracketed by construction)
    assert bool(np.all(profiles.reachable))


def test_oversized_ring_leaves_diamond_is_rejected() -> None:
    import dataclasses

    bad = dataclasses.replace(SMALL, chain_ring_radius=0.9)
    with pytest.raises(SceneValidityError):
        build_scene_2plus1d(bad)


def test_insufficient_targets_recorded_as_validity_error() -> None:
    import dataclasses

    bad = dataclasses.replace(SMALL, min_targets=100000)
    with pytest.raises(SceneValidityError):
        build_scene_2plus1d(bad)
