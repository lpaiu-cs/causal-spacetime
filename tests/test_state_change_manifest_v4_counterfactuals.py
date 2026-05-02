from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_blocking_analysis import (
    V4ProtocolBlockingRecord,
)
from causal_spacetime_lab.state_change_manifest_v4_counterfactuals import (
    family_would_remain_blocked_after_single_fix,
    heldout_threshold_counterfactual_report,
)


def test_report_only_counterfactuals_do_not_change_decisions() -> None:
    rows = heldout_threshold_counterfactual_report(
        [
            {
                "family_name": "family",
                "family_kind": "structured",
                "mean_heldout_violation": 0.4,
            }
        ],
        hypothetical_thresholds=[0.5],
    )
    records = [
        V4ProtocolBlockingRecord(
            "family",
            "structured",
            "mean_heldout_violation",
            0.4,
            0.2,
            "max_allowed",
            "measured_blocking",
            "fail",
            "heldout_failure",
            "fails",
        ),
        V4ProtocolBlockingRecord(
            "family",
            "structured",
            "latent_order_disagreement",
            0.5,
            0.2,
            "max_allowed",
            "measured_blocking",
            "fail",
            "latent_order_instability",
            "fails",
        ),
    ]
    single_fix = family_would_remain_blocked_after_single_fix(
        records,
        ignored_root_cause="heldout_failure",
    )

    assert rows[0]["note"] == "report_only_not_decision_changing"
    assert single_fix[0]["would_remain_blocked"] == 1.0
