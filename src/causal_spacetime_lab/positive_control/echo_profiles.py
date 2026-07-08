"""Order-level bracket-width echo profiles measured from a scene.

The observable is the passive bracket width in tick-index units
(PC-V1 Section 5): W[j, r] = min succeeding tick index - max preceding tick
index on chain r. Clock values are never used, only tick order.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.discrete_radar import find_radar_ticks_from_order
from causal_spacetime_lab.positive_control.scene import PositiveControlScene


@dataclass(frozen=True)
class EchoProfileMatrix:
    """Measured bracket-width echo profiles for a scene's targets."""

    delay_ranks: NDArray[np.float64]
    reachable: NDArray[np.bool_]
    target_indices: NDArray[np.int_]

    @property
    def target_count(self) -> int:
        return int(self.delay_ranks.shape[0])

    @property
    def reference_count(self) -> int:
        return int(self.delay_ranks.shape[1])


def measure_bracket_echo_profiles(scene: PositiveControlScene) -> EchoProfileMatrix:
    """Measure W[j, r] for every target j and reference chain r."""

    n_targets = scene.target_indices.size
    n_chains = len(scene.chain_index_arrays)
    delays = np.full((n_targets, n_chains), np.nan, dtype=np.float64)
    reachable = np.zeros((n_targets, n_chains), dtype=bool)

    for row, target in enumerate(scene.target_indices):
        for column, chain in enumerate(scene.chain_index_arrays):
            ticks = find_radar_ticks_from_order(
                scene.causal,
                chain,
                int(target),
                scene.tick_ranks,
            )
            if ticks is None:
                continue
            rank_minus, rank_plus = ticks
            delays[row, column] = float(rank_plus - rank_minus)
            reachable[row, column] = True

    return EchoProfileMatrix(
        delay_ranks=delays,
        reachable=reachable,
        target_indices=scene.target_indices.copy(),
    )


def shuffle_profile_columns(
    profiles: EchoProfileMatrix,
    seed: int,
) -> EchoProfileMatrix:
    """Destructive marginal-preserving null: permute each column independently.

    Cross-column consistency (the latent spatial signal) is destroyed while
    each column's marginal distribution is preserved (PC-V1 control C-NULL).
    """

    rng = np.random.default_rng(seed)
    delays = profiles.delay_ranks.copy()
    reachable = profiles.reachable.copy()
    for column in range(profiles.reference_count):
        permutation = rng.permutation(profiles.target_count)
        delays[:, column] = delays[permutation, column]
        reachable[:, column] = reachable[permutation, column]
    return EchoProfileMatrix(
        delay_ranks=delays,
        reachable=reachable,
        target_indices=profiles.target_indices.copy(),
    )


def relabel_targets(
    profiles: EchoProfileMatrix,
    truth_coords: NDArray[np.float64],
    seed: int,
) -> tuple[EchoProfileMatrix, NDArray[np.float64]]:
    """Symmetry control C-SYM: one consistent permutation of target labels.

    Profiles and truth labels are permuted together, so every metric must be
    unchanged up to sampling noise.
    """

    rng = np.random.default_rng(seed)
    permutation = rng.permutation(profiles.target_count)
    relabeled = EchoProfileMatrix(
        delay_ranks=profiles.delay_ranks[permutation],
        reachable=profiles.reachable[permutation],
        target_indices=profiles.target_indices[permutation],
    )
    return relabeled, np.asarray(truth_coords, dtype=float)[permutation]
