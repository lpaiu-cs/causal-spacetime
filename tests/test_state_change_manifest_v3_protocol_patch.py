from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)


def test_patched_v3_structured_families_are_protocol_invariant() -> None:
    patches = default_v3_protocol_invariant_family_patches()
    structured = [item for item in patches if item.family_kind == "structured"]

    assert structured
    assert all(item.measurement_protocol_hash for item in structured)
    assert all(item.admissible_for_pairwise_dissimilarity for item in structured)
    assert all(item.execution_status == "planned_only" for item in structured)


def test_failed_controls_v3_is_report_only() -> None:
    patches = default_v3_protocol_invariant_family_patches()
    failed = [
        item for item in patches if item.patched_family_name == "failed_controls_v3"
    ]

    assert failed
    assert failed[0].family_kind == "failed_control"
    assert failed[0].admissible_for_pairwise_dissimilarity is False
