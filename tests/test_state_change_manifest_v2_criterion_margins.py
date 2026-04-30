from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    V2BlockingCriterionRecord,
)
from causal_spacetime_lab.state_change_manifest_v2_criterion_margins import (
    criterion_margins_from_blocking_records,
)


def test_criterion_margins_sign_convention() -> None:
    records = [
        V2BlockingCriterionRecord(
            family_name="f",
            family_kind="structured",
            criterion_name="manifest_count",
            observed_value=5.0,
            threshold_value=3.0,
            criterion_direction="min_required",
            blocking_type="not_blocking",
            status="pass",
            explanation="pass",
        ),
        V2BlockingCriterionRecord(
            family_name="f",
            family_kind="structured",
            criterion_name="heldout",
            observed_value=0.5,
            threshold_value=0.2,
            criterion_direction="max_allowed",
            blocking_type="measured_blocking",
            status="fail",
            explanation="fail",
        ),
    ]

    margins = criterion_margins_from_blocking_records(records)

    assert margins[0].margin > 0
    assert margins[1].margin < 0

