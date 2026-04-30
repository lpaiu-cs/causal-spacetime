from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    V3ProtocolBlockingRecord,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_criterion_margins import (
    criterion_margins_from_v3_protocol_blocking_records,
    summarize_v3_protocol_criterion_margins,
)


def _record(status: str, observed: float, threshold: float) -> V3ProtocolBlockingRecord:
    return V3ProtocolBlockingRecord(
        family_name="family",
        family_kind="structured",
        criterion_name="mean_heldout_violation",
        observed_value=observed,
        threshold_value=threshold,
        criterion_direction="max_allowed",
        blocking_type="measured_blocking" if status == "fail" else "not_blocking",
        status=status,
        root_cause_category="heldout_failure",
        explanation="test",
    )


def test_criterion_margins_use_correct_sign_convention() -> None:
    margins = criterion_margins_from_v3_protocol_blocking_records(
        [_record("pass", 0.1, 0.2), _record("fail", 0.3, 0.2)]
    )

    assert margins[0].margin > 0
    assert margins[1].margin < 0


def test_summarize_v3_protocol_criterion_margins_counts_failures() -> None:
    summary = summarize_v3_protocol_criterion_margins(
        criterion_margins_from_v3_protocol_blocking_records(
            [_record("fail", 0.3, 0.2)]
        )
    )

    assert summary[0]["failed_criterion_count"] == 1.0
    assert summary[0]["measured_failure_count"] == 1.0
