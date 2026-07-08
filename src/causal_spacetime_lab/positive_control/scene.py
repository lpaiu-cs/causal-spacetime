"""Scene construction for the PC-V1 positive-control experiment.

A scene is a sprinkled 1+1D Minkowski causal set plus supplied stationary
observer tick chains, with a declared target set and content-addressed
provenance. See docs/prereg/pc_v1_positive_control.md Section 4.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.discrete_radar import find_radar_ticks_from_order
from causal_spacetime_lab.observer import (
    inertial_chain_inside_diamond_mask,
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond


class SceneValidityError(ValueError):
    """Raised when a seed fails a preregistered scene-validity precondition."""


@dataclass(frozen=True)
class PositiveControlSceneConfig:
    """Preregistered primary configuration (PC-V1 Section 4)."""

    n_events: int = 900
    diamond_T: float = 2.0
    chain_positions: tuple[float, ...] = (-0.25, -0.15, -0.05, 0.05, 0.15, 0.25)
    ticks_per_chain: int = 96
    chain_span: float = 1.4
    target_band_t: float = 0.10
    target_band_x: float = 0.25
    max_targets: int = 40
    min_targets: int = 30
    min_bracketing_chains: int = 6
    seed: int = 0

    def config_hash(self) -> str:
        """Return a stable hash of the configuration."""

        payload = json.dumps(asdict(self), sort_keys=True).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()[:16]


@dataclass(frozen=True)
class PositiveControlScene:
    """Combined event set, causal order, references, targets, provenance."""

    config: PositiveControlSceneConfig
    events: NDArray[np.float64]
    causal: NDArray[np.bool_]
    bulk_count: int
    chain_index_arrays: tuple[NDArray[np.int_], ...]
    tick_ranks: NDArray[np.float64]
    target_indices: NDArray[np.int_]
    events_digest: str
    causal_digest: str

    @property
    def target_coords(self) -> NDArray[np.float64]:
        """Return true (t, x) coordinates of targets (validation labels)."""

        return self.events[self.target_indices]

    def provenance_row(self) -> dict[str, str | float]:
        """Return provenance fields attached to every output row."""

        return {
            "config_hash": self.config.config_hash(),
            "events_digest": self.events_digest,
            "causal_digest": self.causal_digest,
            "scene_seed": float(self.config.seed),
            "target_count": float(self.target_indices.size),
        }


def _digest(array: NDArray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()[:16]


def two_sided_bracket_counts(
    causal: NDArray[np.bool_],
    chain_index_arrays: tuple[NDArray[np.int_], ...],
    tick_ranks: NDArray[np.float64],
    candidate_indices: NDArray[np.int_],
) -> NDArray[np.int_]:
    """Count chains that two-sided bracket each candidate target."""

    counts = np.zeros(candidate_indices.size, dtype=int)
    for position, target in enumerate(candidate_indices):
        for chain in chain_index_arrays:
            ticks = find_radar_ticks_from_order(causal, chain, int(target), tick_ranks)
            if ticks is not None:
                counts[position] += 1
    return counts


def build_positive_control_scene(
    config: PositiveControlSceneConfig,
) -> PositiveControlScene:
    """Build a scene per PC-V1 Section 4, enforcing validity preconditions."""

    if config.min_bracketing_chains > len(config.chain_positions):
        raise ValueError("min_bracketing_chains exceeds chain count")

    bulk = sprinkle_1p1_causal_diamond(
        config.n_events,
        T=config.diamond_T,
        seed=config.seed,
    )

    chain_blocks: list[NDArray[np.float64]] = []
    for x0 in config.chain_positions:
        chain_events, _ = make_stationary_observer_chain_1p1(
            config.chain_span,
            config.ticks_per_chain,
            x=float(x0),
        )
        inside = inertial_chain_inside_diamond_mask(chain_events, config.diamond_T)
        if not bool(np.all(inside)):
            raise SceneValidityError(
                f"chain at x0={x0} leaves the diamond; adjust chain_span"
            )
        chain_blocks.append(chain_events)

    events = np.vstack([bulk, *chain_blocks])
    chain_index_arrays = []
    start = bulk.shape[0]
    for block in chain_blocks:
        chain_index_arrays.append(observer_chain_indices(start, block.shape[0]))
        start += block.shape[0]

    causal = causal_matrix_1p1(events)
    tick_ranks = np.arange(config.ticks_per_chain, dtype=np.float64)

    in_band = (np.abs(bulk[:, 0]) <= config.target_band_t) & (
        np.abs(bulk[:, 1]) <= config.target_band_x
    )
    candidates = np.flatnonzero(in_band)
    bracket_counts = two_sided_bracket_counts(
        causal,
        tuple(chain_index_arrays),
        tick_ranks,
        candidates,
    )
    eligible = candidates[bracket_counts >= config.min_bracketing_chains]

    if eligible.size < config.min_targets:
        raise SceneValidityError(
            f"only {eligible.size} eligible targets "
            f"(minimum {config.min_targets}); seed recorded as invalid"
        )
    if eligible.size > config.max_targets:
        subsample_rng = np.random.default_rng(config.seed + 1)
        chosen = subsample_rng.choice(
            eligible.size,
            size=config.max_targets,
            replace=False,
        )
        eligible = np.sort(eligible[chosen])

    return PositiveControlScene(
        config=config,
        events=events,
        causal=causal,
        bulk_count=bulk.shape[0],
        chain_index_arrays=tuple(chain_index_arrays),
        tick_ranks=tick_ranks,
        target_indices=eligible.astype(int),
        events_digest=_digest(events),
        causal_digest=_digest(causal),
    )
