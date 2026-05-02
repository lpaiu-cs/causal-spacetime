from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_blocking_analysis import (
    V4ProtocolBlockingRecord,
)
from causal_spacetime_lab.state_change_manifest_v4_criterion_margins import (
    criterion_margins_from_v4_blocking_records,
)


def test_criterion_margins_use_correct_sign_convention() -> None:
    records = [
        V4ProtocolBlockingRecord(
            "family",
            "structured",
            "min_metric",
            0.8,
            1.0,
            "min_required",
            "measured_blocking",
            "fail",
            "coverage_failure",
            "fails",
        ),
        V4ProtocolBlockingRecord(
            "family",
            "structured",
            "max_metric",
            0.1,
            0.2,
            "max_allowed",
            "not_blocking",
            "pass",
            "heldout_failure",
            "passes",
        ),
        V4ProtocolBlockingRecord(
            "family",
            "structured",
            "boolean_metric",
            0.0,
            1.0,
            "boolean_required",
            "measured_blocking",
            "fail",
            "unknown",
            "fails",
        ),
    ]

    margins = criterion_margins_from_v4_blocking_records(records)

    assert margins[0].margin < 0
    assert margins[1].margin > 0
    assert margins[2].margin == -1.0
