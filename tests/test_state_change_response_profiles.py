from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_profiles import (
    compare_response_profiles,
    response_profile_equivalence_classes,
    response_profile_from_signatures,
    response_profile_separation_fraction,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)


def _signature(
    targets: list[int],
    delays: list[int],
    label: str,
) -> EchoResponseSignature:
    target_array = np.asarray(targets, dtype=int)
    delay_array = np.asarray(delays, dtype=int)
    reachable = delay_array >= 0
    return EchoResponseSignature(
        target_event_ids=target_array,
        delay_ranks=delay_array,
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delay_array, reachable),
        label=label,
    )


def test_response_profile_from_signatures_aligns_common_targets() -> None:
    profile = response_profile_from_signatures(
        [
            _signature([1, 2, 3], [2, 3, 4], "a"),
            _signature([2, 3, 4], [5, 6, 7], "b"),
        ]
    )

    assert np.array_equal(profile.target_event_ids, np.asarray([2, 3]))
    assert profile.delay_rank_matrix.shape == (2, 2)


def test_response_profile_separation_fraction_deterministic() -> None:
    profile = response_profile_from_signatures(
        [
            _signature([1, 2, 3], [2, 2, 2], "a"),
            _signature([1, 2, 3], [2, 3, 2], "b"),
        ]
    )

    assert response_profile_separation_fraction(profile) > 0.0


def test_response_profile_equivalence_classes_deterministic() -> None:
    profile = response_profile_from_signatures(
        [
            _signature([1, 2, 3], [2, 2, 4], "a"),
            _signature([1, 2, 3], [5, 5, 7], "b"),
        ]
    )

    classes = response_profile_equivalence_classes(profile)

    assert any(set(group.tolist()) == {1, 2} for group in classes)


def test_compare_response_profiles_identical_case() -> None:
    profile = response_profile_from_signatures(
        [
            _signature([1, 2, 3], [2, 3, 4], "a"),
            _signature([1, 2, 3], [5, 6, 7], "b"),
        ]
    )

    comparison = compare_response_profiles(profile, profile)

    assert comparison["profile_entry_agreement_fraction"] == 1.0
    assert comparison["reachable_agreement_fraction"] == 1.0
