from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_failure_decomposition import (
    decompose_v4_protocol_family_failures,
)
from tests.test_state_change_manifest_v4_carry_forward import _precondition, _row


def test_decompose_v4_protocol_family_failures_separates_metric_failures() -> None:
    records = decompose_v4_protocol_family_failures(
        [_row(heldout=0.5)],
        default_cross_family_robustness_criteria(),
        [_precondition()],
    )

    assert any(record.root_cause_category == "heldout_failure" for record in records)


def test_decompose_v4_protocol_family_failures_separates_provenance_failures() -> None:
    records = decompose_v4_protocol_family_failures(
        [_row()],
        default_cross_family_robustness_criteria(),
        [_precondition(passed=False)],
    )

    assert any(record.root_cause_category == "provenance_failure" for record in records)
