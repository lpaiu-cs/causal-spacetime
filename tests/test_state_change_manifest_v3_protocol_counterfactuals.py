from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    V3ProtocolBlockingRecord,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_counterfactuals import (
    family_would_remain_blocked_after_single_fix,
    heldout_threshold_counterfactual_report,
)


def test_report_only_counterfactuals_do_not_change_decisions() -> None:
    rows = heldout_threshold_counterfactual_report(
        [
            {
                "family_name": "family",
                "family_kind": "structured",
                "mean_heldout_violation": 0.3,
            }
        ],
        hypothetical_thresholds=[0.2],
    )

    assert rows[0]["output_note"] == "report_only_not_decision_changing"


def test_single_fix_counterfactual_detects_remaining_blocker() -> None:
    records = [
        V3ProtocolBlockingRecord(
            family_name="family",
            family_kind="structured",
            criterion_name="mean_heldout_violation",
            observed_value=0.3,
            threshold_value=0.2,
            criterion_direction="max_allowed",
            blocking_type="measured_blocking",
            status="fail",
            root_cause_category="heldout_failure",
            explanation="test",
        ),
        V3ProtocolBlockingRecord(
            family_name="family",
            family_kind="structured",
            criterion_name="latent_order_disagreement",
            observed_value=0.6,
            threshold_value=0.3,
            criterion_direction="max_allowed",
            blocking_type="measured_blocking",
            status="fail",
            root_cause_category="latent_order_instability",
            explanation="test",
        ),
    ]

    rows = family_would_remain_blocked_after_single_fix(
        records,
        ignored_root_cause="heldout_failure",
    )

    assert rows[0]["would_remain_blocked"] == 1.0
