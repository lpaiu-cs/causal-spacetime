"""Scene construction tests for the PC-V1 positive-control pipeline."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from causal_spacetime_lab.positive_control.scene import (
    PositiveControlSceneConfig,
    SceneValidityError,
    build_positive_control_scene,
)

SMALL_CONFIG = PositiveControlSceneConfig(
    n_events=400,
    chain_positions=(-0.2, -0.07, 0.07, 0.2),
    ticks_per_chain=32,
    target_band_t=0.10,
    target_band_x=0.20,
    max_targets=16,
    min_targets=6,
    min_bracketing_chains=4,
    seed=3,
)


def test_scene_is_deterministic_and_content_addressed() -> None:
    first = build_positive_control_scene(SMALL_CONFIG)
    second = build_positive_control_scene(SMALL_CONFIG)
    assert first.events_digest == second.events_digest
    assert first.causal_digest == second.causal_digest
    assert np.array_equal(first.target_indices, second.target_indices)

    other_seed = build_positive_control_scene(
        dataclasses.replace(SMALL_CONFIG, seed=SMALL_CONFIG.seed + 1)
    )
    assert other_seed.events_digest != first.events_digest


def test_chain_events_form_causal_chains() -> None:
    scene = build_positive_control_scene(SMALL_CONFIG)
    for chain in scene.chain_index_arrays:
        for earlier in range(chain.size - 1):
            assert scene.causal[chain[earlier], chain[earlier + 1]]
        assert not scene.causal[chain[-1], chain[0]]


def test_targets_respect_band_and_bracketing() -> None:
    scene = build_positive_control_scene(SMALL_CONFIG)
    coords = scene.target_coords
    assert np.all(np.abs(coords[:, 0]) <= SMALL_CONFIG.target_band_t)
    assert np.all(np.abs(coords[:, 1]) <= SMALL_CONFIG.target_band_x)
    assert scene.target_indices.size >= SMALL_CONFIG.min_targets
    assert scene.target_indices.size <= SMALL_CONFIG.max_targets
    assert np.all(scene.target_indices < scene.bulk_count)


def test_chain_leaving_diamond_is_rejected() -> None:
    config = dataclasses.replace(
        SMALL_CONFIG, chain_positions=(0.9,), min_bracketing_chains=1
    )
    with pytest.raises(SceneValidityError):
        build_positive_control_scene(config)


def test_insufficient_targets_is_a_recorded_validity_error() -> None:
    config = dataclasses.replace(SMALL_CONFIG, min_targets=10_000)
    with pytest.raises(SceneValidityError):
        build_positive_control_scene(config)
