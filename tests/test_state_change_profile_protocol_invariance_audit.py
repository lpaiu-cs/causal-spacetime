from __future__ import annotations

from causal_spacetime_lab.state_change_measurement_protocol import (
    default_earliest_full_protocol,
)
from causal_spacetime_lab.state_change_profile_protocol_invariance_audit import (
    audit_manifest_protocol_metadata,
    audit_profile_metadata,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_from_protocols,
)


def test_audit_manifest_protocol_metadata_underspecified_for_old_manifest() -> None:
    row = audit_manifest_protocol_metadata(
        {"manifest_id": "old", "profile_label": "old_profile"}
    )

    assert row.profile_invariance_status == "underspecified"
    assert not row.admissible_for_pairwise_dissimilarity


def test_audit_profile_metadata_protocol_invariant() -> None:
    protocol = default_earliest_full_protocol()
    metadata = profile_metadata_from_protocols(
        "family",
        [protocol, protocol],
        ["r1", "r2"],
        "refs",
    )
    row = audit_profile_metadata(metadata)

    assert row.profile_invariance_status == "protocol_invariant"
    assert row.admissible_for_pairwise_dissimilarity
