"""Null baselines for response-profile pairwise comparison diagnostics."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _copy_profile_arrays(
    profile: EchoResponseProfile,
) -> tuple[np.ndarray, np.ndarray]:
    return profile.delay_rank_matrix.copy(), profile.reachable_matrix.copy()


def shuffle_profile_delays_within_protocols(
    profile: EchoResponseProfile,
    seed: int | None = None,
) -> EchoResponseProfile:
    """Shuffle reachable delay ranks within each protocol column."""

    rng = np.random.default_rng(seed)
    delays, reachable = _copy_profile_arrays(profile)
    for column in range(delays.shape[1]):
        mask = reachable[:, column]
        values = delays[mask, column].copy()
        rng.shuffle(values)
        delays[mask, column] = values
    return EchoResponseProfile(
        profile.target_event_ids.copy(),
        list(profile.protocol_labels),
        delays,
        reachable,
    )


def shuffle_profile_reachability_masks(
    profile: EchoResponseProfile,
    seed: int | None = None,
) -> EchoResponseProfile:
    """Shuffle reachability masks while preserving per-column counts."""

    rng = np.random.default_rng(seed)
    original_delays = profile.delay_rank_matrix
    delays = np.full_like(original_delays, -1)
    reachable = profile.reachable_matrix.copy()
    for column in range(reachable.shape[1]):
        rng.shuffle(reachable[:, column])
        source_values = original_delays[profile.reachable_matrix[:, column], column]
        values = source_values.copy()
        rng.shuffle(values)
        delays[reachable[:, column], column] = values[: np.sum(reachable[:, column])]
    return EchoResponseProfile(
        profile.target_event_ids.copy(),
        list(profile.protocol_labels),
        delays,
        reachable,
    )


def permute_target_profiles(
    profile: EchoResponseProfile,
    independently_per_protocol: bool = True,
    seed: int | None = None,
) -> EchoResponseProfile:
    """Permute target rows to break cross-protocol profile alignment."""

    rng = np.random.default_rng(seed)
    delays, reachable = _copy_profile_arrays(profile)
    if independently_per_protocol:
        for column in range(delays.shape[1]):
            permutation = rng.permutation(delays.shape[0])
            delays[:, column] = delays[permutation, column]
            reachable[:, column] = reachable[permutation, column]
    else:
        permutation = rng.permutation(delays.shape[0])
        delays = delays[permutation]
        reachable = reachable[permutation]
    return EchoResponseProfile(
        profile.target_event_ids.copy(),
        list(profile.protocol_labels),
        delays,
        reachable,
    )


def random_profile_with_same_marginals(
    profile: EchoResponseProfile,
    seed: int | None = None,
) -> EchoResponseProfile:
    """Build a random profile preserving per-column reachable delay multisets."""

    rng = np.random.default_rng(seed)
    n_targets, n_protocols = profile.delay_rank_matrix.shape
    delays = np.full((n_targets, n_protocols), -1, dtype=int)
    reachable = np.zeros((n_targets, n_protocols), dtype=bool)
    for column in range(n_protocols):
        source_mask = profile.reachable_matrix[:, column]
        values = profile.delay_rank_matrix[source_mask, column].copy()
        rng.shuffle(values)
        selected = rng.choice(n_targets, size=values.size, replace=False)
        reachable[selected, column] = True
        delays[selected, column] = values
    return EchoResponseProfile(
        profile.target_event_ids.copy(),
        list(profile.protocol_labels),
        delays,
        reachable,
    )
