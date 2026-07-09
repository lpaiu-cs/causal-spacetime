"""2+1D scene construction for the P2 robustness experiment.

A 2+1D scene sprinkles a 3D Minkowski causal diamond and places stationary
observer reference chains on a ring of 2D spatial positions (non-collinear, so
cross-observer parallax determines a 2D position). It reuses the dimension-
agnostic PositiveControlScene container and the frozen echo/dissimilarity/fit
pipeline; only scene construction is 2+1D-specific. See
docs/prereg/p2_2plus1d.md.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.causal import causal_matrix_minkowski
from causal_spacetime_lab.discrete_radar import find_radar_ticks_from_order
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlScene,
    PositiveControlSceneConfig,
    SceneValidityError,
)
from causal_spacetime_lab.sprinkling import sprinkle_minkowski_causal_diamond


@dataclass(frozen=True)
class Scene2DConfig:
    """Preregistered 2+1D configuration (P2)."""

    n_events: int = 2600
    diamond_T: float = 2.0
    chain_count: int = 8
    chain_ring_radius: float = 0.25
    ticks_per_chain: int = 96
    chain_span: float = 1.4
    target_band_t: float = 0.10
    target_band_radius: float = 0.22
    max_targets: int = 44
    min_targets: int = 30
    min_bracketing_chains: int = 8
    seed: int = 0

    def as_scene_config(self) -> PositiveControlSceneConfig:
        """A 1+1D-style config object carrying the shared fields for reuse."""

        return PositiveControlSceneConfig(
            n_events=self.n_events,
            diamond_T=self.diamond_T,
            ticks_per_chain=self.ticks_per_chain,
            chain_span=self.chain_span,
            target_band_t=self.target_band_t,
            target_band_x=self.target_band_radius,
            max_targets=self.max_targets,
            min_targets=self.min_targets,
            min_bracketing_chains=self.min_bracketing_chains,
            seed=self.seed,
        )


def _digest(array: NDArray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()[:16]


def _ring_positions(config: Scene2DConfig) -> NDArray[np.float64]:
    angles = np.linspace(0.0, 2.0 * np.pi, config.chain_count, endpoint=False)
    return config.chain_ring_radius * np.column_stack(
        (np.cos(angles), np.sin(angles))
    )


def build_scene_2plus1d(config: Scene2DConfig) -> PositiveControlScene:
    """Build a 2+1D scene, enforcing the preregistered validity preconditions."""

    if config.min_bracketing_chains > config.chain_count:
        raise ValueError("min_bracketing_chains exceeds chain count")

    bulk = sprinkle_minkowski_causal_diamond(
        config.n_events, spacetime_dim=3, T=config.diamond_T, seed=config.seed
    )
    clock = np.linspace(
        -config.chain_span / 2.0, config.chain_span / 2.0, config.ticks_per_chain
    )
    positions = _ring_positions(config)
    half_t = config.diamond_T / 2.0
    chain_blocks: list[NDArray[np.float64]] = []
    for x0, y0 in positions:
        radius = float(np.hypot(x0, y0))
        # stationary chain worldline stays inside the diamond
        if radius > half_t - np.max(np.abs(clock)) + 1e-12:
            raise SceneValidityError(
                f"chain at radius {radius} leaves the diamond; reduce ring radius"
            )
        chain_blocks.append(
            np.column_stack(
                (clock, np.full(clock.size, x0), np.full(clock.size, y0))
            )
        )

    events = np.vstack([bulk, *chain_blocks])
    causal = causal_matrix_minkowski(events)
    tick_ranks = np.arange(config.ticks_per_chain, dtype=np.float64)
    chain_index_arrays = []
    start = bulk.shape[0]
    for block in chain_blocks:
        chain_index_arrays.append(np.arange(start, start + block.shape[0], dtype=int))
        start += block.shape[0]

    radius = np.hypot(bulk[:, 1], bulk[:, 2])
    in_band = (np.abs(bulk[:, 0]) <= config.target_band_t) & (
        radius <= config.target_band_radius
    )
    candidates = np.flatnonzero(in_band)
    eligible = [
        int(target)
        for target in candidates
        if all(
            find_radar_ticks_from_order(causal, chain, int(target), tick_ranks)
            is not None
            for chain in chain_index_arrays
        )
    ]
    eligible = np.array(eligible, dtype=int)
    if eligible.size < config.min_targets:
        raise SceneValidityError(
            f"2+1D scene: only {eligible.size} eligible targets "
            f"(minimum {config.min_targets})"
        )
    if eligible.size > config.max_targets:
        subsample_rng = np.random.default_rng(config.seed + 1)
        chosen = subsample_rng.choice(
            eligible.size, size=config.max_targets, replace=False
        )
        eligible = np.sort(eligible[chosen])

    return PositiveControlScene(
        config=config.as_scene_config(),
        events=events,
        causal=causal,
        bulk_count=bulk.shape[0],
        chain_index_arrays=tuple(chain_index_arrays),
        tick_ranks=tick_ranks,
        target_indices=eligible.astype(int),
        events_digest=_digest(events),
        causal_digest=_digest(causal),
    )


def target_positions_2d(scene: PositiveControlScene) -> NDArray[np.float64]:
    """Return true 2D spatial coordinates (x, y) of the scene's targets."""

    return scene.events[scene.target_indices][:, 1:3]
