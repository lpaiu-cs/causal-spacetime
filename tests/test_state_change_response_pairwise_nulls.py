from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_pairwise_nulls import (
    permute_target_profiles,
    random_profile_with_same_marginals,
    shuffle_profile_delays_within_protocols,
    shuffle_profile_reachability_masks,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _profile() -> EchoResponseProfile:
    delays = np.asarray(
        [[1, 2, 3], [2, -1, 3], [3, 2, -1], [1, 4, 5]],
        dtype=int,
    )
    return EchoResponseProfile(
        target_event_ids=np.arange(4, dtype=int),
        protocol_labels=["a", "b", "c"],
        delay_rank_matrix=delays,
        reachable_matrix=delays >= 0,
    )


def _column_multisets(profile: EchoResponseProfile) -> list[list[int]]:
    return [
        sorted(
            int(value)
            for value in profile.delay_rank_matrix[
                profile.reachable_matrix[:, column],
                column,
            ]
        )
        for column in range(profile.delay_rank_matrix.shape[1])
    ]


def _reachable_counts(profile: EchoResponseProfile) -> list[int]:
    return [
        int(np.sum(profile.reachable_matrix[:, column]))
        for column in range(profile.reachable_matrix.shape[1])
    ]


def test_shuffle_profile_delays_within_protocols_preserves_multisets() -> None:
    profile = _profile()
    shuffled = shuffle_profile_delays_within_protocols(profile, seed=0)

    assert _column_multisets(shuffled) == _column_multisets(profile)


def test_shuffle_profile_reachability_masks_preserves_counts() -> None:
    profile = _profile()
    shuffled = shuffle_profile_reachability_masks(profile, seed=0)

    assert _reachable_counts(shuffled) == _reachable_counts(profile)


def test_permute_target_profiles_preserves_column_marginals() -> None:
    profile = _profile()
    permuted = permute_target_profiles(profile, seed=0)

    assert _column_multisets(permuted) == _column_multisets(profile)
    assert _reachable_counts(permuted) == _reachable_counts(profile)


def test_random_profile_with_same_marginals_preserves_marginals() -> None:
    profile = _profile()
    random_profile = random_profile_with_same_marginals(profile, seed=0)

    assert _column_multisets(random_profile) == _column_multisets(profile)
    assert _reachable_counts(random_profile) == _reachable_counts(profile)
