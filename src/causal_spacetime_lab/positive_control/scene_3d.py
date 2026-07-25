"""3+1D scene construction for the P8 dimension-selection experiment.

A 3+1D scene sprinkles a 4D Minkowski causal diamond and places stationary
observer reference chains on a sphere of 3D spatial positions. It reuses the
dimension-agnostic PositiveControlScene container and the frozen
echo/dissimilarity/fit pipeline; only scene construction is 3+1D-specific,
exactly as `scene_2d` is for P2. See docs/prereg/p8_3plus1d.md.

Two design choices differ from P2 and both are deliberate.

*Observer placement* is a Fibonacci sphere rather than a ring. In three
spatial dimensions the labeled decoder is multilateration, which needs
observers that affinely span -- a coplanar set is exactly the degenerate
case -- so the builder asserts affine rank 3 rather than trusting the
construction. A Fibonacci spiral is used in preference to a highly
symmetric shell because symmetric configurations are non-generic, and a
rank measured on one understates what a typical layout achieves.

*The target band is wider*, 0.30 against P2's 0.22. This is not a
loosening: the band is a solid region, so at fixed radius it captures a
smaller fraction of a 4D diamond than of a 3D one, and 0.30 restores the
candidate count to P2's range at a comparable event budget. The
bracketing filter turns out not to bind at all here -- every in-band
candidate is two-sided bracketed by every chain, measured across the
whole design grid -- so the candidate count is the only lever.

*And there are twelve chains, not P2's eight.* This one was measured
rather than assumed, because transplanting P2's layout into 3+1D leaves
almost no room between the correct dimension and one short of it. On an
exploratory sweep (seeds 0-3, four seeds, no gates set from it) the
separation between the `d = 3` and `d = 2` truth-error clusters ran

    R  =    6      8     12     16     20
    gap = 0.015  0.016  0.035  0.038  0.031

against P2's 0.074 in its own dimension. There is a step of about two
between `R <= 8` and `R >= 12`, and a turnover by 20; four seeds cannot
separate 12 from 16, and this note does not claim they differ. Twelve is
chosen as the cheap end of the plateau. It also happens to be where the
independent conditioning study (Section 5c of the theory document) put
the margin optimum for three spatial dimensions, which is corroboration
rather than evidence -- that study measured a different quantity.

Two levers that did NOT work, recorded so they are not retried: doubling
the tick count moved the separation not at all (0.035 -> 0.034), and
raising the target count to 80 halved it (0.035 -> 0.016), because the
frozen fit policy spends a fixed constraint budget and more targets
simply divide it more thinly.
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
class Scene3DConfig:
    """Preregistered 3+1D configuration (P8)."""

    n_events: int = 7000
    diamond_T: float = 2.0
    chain_count: int = 12
    chain_shell_radius: float = 0.25
    ticks_per_chain: int = 96
    chain_span: float = 1.4
    target_band_t: float = 0.10
    target_band_radius: float = 0.30
    max_targets: int = 44
    min_targets: int = 30
    min_bracketing_chains: int = 12
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


def fibonacci_shell(count: int, radius: float) -> NDArray[np.float64]:
    """``count`` quasi-uniform points on the sphere of the given radius.

    Deterministic and seed-free, which matters for a configuration that is
    going to be frozen: the layout is reproducible from the count alone.
    """

    if count < 4:
        raise ValueError("a 3D observer shell needs at least 4 chains to span")
    indices = np.arange(count, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * indices / count
    rho = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    golden = np.pi * (1.0 + 5.0**0.5)
    phi = golden * indices
    return radius * np.column_stack((rho * np.cos(phi), rho * np.sin(phi), z))


def affine_rank(points: NDArray[np.float64]) -> int:
    centered = points - points.mean(axis=0, keepdims=True)
    spectrum = np.linalg.svd(centered, compute_uv=False)
    if spectrum.size == 0 or spectrum[0] == 0.0:
        return 0
    return int(np.sum(spectrum > spectrum[0] * 1e-9))


def build_scene_3plus1d(config: Scene3DConfig) -> PositiveControlScene:
    """Build a 3+1D scene, enforcing the preregistered validity preconditions."""

    if config.min_bracketing_chains > config.chain_count:
        raise ValueError("min_bracketing_chains exceeds chain count")

    bulk = sprinkle_minkowski_causal_diamond(
        config.n_events, spacetime_dim=4, T=config.diamond_T, seed=config.seed
    )
    clock = np.linspace(
        -config.chain_span / 2.0, config.chain_span / 2.0, config.ticks_per_chain
    )
    positions = fibonacci_shell(config.chain_count, config.chain_shell_radius)
    if affine_rank(positions) != 3:
        raise SceneValidityError(
            "observer shell is degenerate: the chains do not span three "
            "spatial dimensions, so multilateration is not determined"
        )

    half_t = config.diamond_T / 2.0
    chain_blocks: list[NDArray[np.float64]] = []
    for position in positions:
        radius = float(np.linalg.norm(position))
        # stationary chain worldline stays inside the diamond
        if radius > half_t - np.max(np.abs(clock)) + 1e-12:
            raise SceneValidityError(
                f"chain at radius {radius} leaves the diamond; reduce shell radius"
            )
        chain_blocks.append(
            np.column_stack(
                (clock, *[np.full(clock.size, c) for c in position])
            )
        )

    events = np.vstack([bulk, *chain_blocks])
    causal = causal_matrix_minkowski(events)
    tick_ranks = np.arange(config.ticks_per_chain, dtype=np.float64)
    chain_index_arrays = []
    start = bulk.shape[0]
    for block in chain_blocks:
        chain_index_arrays.append(
            np.arange(start, start + block.shape[0], dtype=int)
        )
        start += block.shape[0]

    radius = np.linalg.norm(bulk[:, 1:], axis=1)
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
            f"3+1D scene: only {eligible.size} eligible targets "
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


def target_positions_3d(scene: PositiveControlScene) -> NDArray[np.float64]:
    """Return true 3D spatial coordinates (x, y, z) of the scene's targets."""

    return scene.events[scene.target_indices][:, 1:4]
