from __future__ import annotations

from dataclasses import asdict

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_threshold_sensitivity import (
    v4_protocol_threshold_sensitivity_rows,
)
from tests.test_state_change_manifest_v4_carry_forward import _precondition, _row


def test_v4_protocol_threshold_sensitivity_rows_does_not_mutate_default() -> None:
    before = asdict(default_cross_family_robustness_criteria())
    rows = v4_protocol_threshold_sensitivity_rows(
        [_row()],
        [_precondition()],
        heldout_thresholds=[0.2],
        null_gap_thresholds=[0.1],
        stricter_pass_thresholds=[0.5],
    )
    after = asdict(default_cross_family_robustness_criteria())

    assert rows[0]["note"] == "sensitivity_only_not_retuning"
    assert before == after
