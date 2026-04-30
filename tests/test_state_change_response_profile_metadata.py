from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_measurement_protocol import (
    default_earliest_full_protocol,
    default_gated_full_protocol,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity_checked,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_from_protocols,
)
from causal_spacetime_lab.state_change_response_profiles import (
    EchoResponseProfile,
    attach_profile_metadata,
    profile_admissible_for_pairwise_dissimilarity,
)


def _profile() -> EchoResponseProfile:
    return EchoResponseProfile(
        target_event_ids=np.asarray([1, 2, 3], dtype=int),
        protocol_labels=["r1", "r2"],
        delay_rank_matrix=np.asarray([[1, 2], [2, 3], [4, 5]], dtype=int),
        reachable_matrix=np.ones((3, 2), dtype=bool),
    )


def test_protocol_invariant_profile_metadata_is_admissible() -> None:
    protocol = default_earliest_full_protocol()
    metadata = profile_metadata_from_protocols(
        "family",
        [protocol, protocol],
        ["r1", "r2"],
        "refs",
    )
    wrapped = attach_profile_metadata(_profile(), metadata)

    assert profile_admissible_for_pairwise_dissimilarity(wrapped)
    result = pairwise_response_dissimilarity_checked(
        wrapped,
        PairwiseResponseComparisonProtocol(name="gap", method="rank_gap_mean"),
    )
    assert result.dissimilarity_values.size == 3


def test_protocol_mixed_profile_metadata_is_blocked() -> None:
    metadata = profile_metadata_from_protocols(
        "family",
        [default_earliest_full_protocol(), default_gated_full_protocol()],
        ["r1", "r2"],
        "refs",
    )
    wrapped = attach_profile_metadata(_profile(), metadata)

    assert not profile_admissible_for_pairwise_dissimilarity(wrapped)
    with pytest.raises(ValueError):
        pairwise_response_dissimilarity_checked(
            wrapped,
            PairwiseResponseComparisonProtocol(name="gap", method="rank_gap_mean"),
        )


def test_underspecified_metadata_is_blocked() -> None:
    metadata = profile_metadata_from_protocols("family", [], [], "")

    assert metadata.profile_invariance_status == "underspecified"
    assert not metadata.admissible_for_pairwise_dissimilarity


def test_exploratory_mixed_context_is_report_only() -> None:
    metadata = profile_metadata_from_protocols(
        "family",
        [default_earliest_full_protocol(), default_gated_full_protocol()],
        ["r1", "r2"],
        "refs",
        exploratory_mixed_context=True,
    )

    assert metadata.exploratory_mixed_context
    assert metadata.profile_invariance_status == "protocol_mixed"
    assert not metadata.admissible_for_pairwise_dissimilarity
