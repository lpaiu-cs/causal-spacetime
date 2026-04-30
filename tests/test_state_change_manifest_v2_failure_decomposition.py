from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_failure_decomposition import (
    decompose_v2_family_failures,
    summarize_v2_failure_records,
    v2_missing_metric_impact_rows,
)


def test_decompose_v2_family_failures_separates_measured_and_missing() -> None:
    rows = [
        {
            "family_name": "family",
            "family_kind": "structured",
            "manifest_count": 1.0,
            "mean_heldout_violation": 0.5,
        }
    ]

    records = decompose_v2_family_failures(
        rows,
        default_cross_family_robustness_criteria(),
    )
    statuses = {record.status for record in records}

    assert "measured_failure" in statuses
    assert "missing_metric" in statuses
    assert summarize_v2_failure_records(records)
    assert v2_missing_metric_impact_rows(records)
