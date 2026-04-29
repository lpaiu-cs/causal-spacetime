from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_diagnostic_completeness import (
    diagnostic_completeness_for_family,
    diagnostic_completeness_table,
    required_cross_family_metrics,
)


def test_diagnostic_completeness_for_family_detects_missing_metrics() -> None:
    report = diagnostic_completeness_for_family(
        {"family_name": "partial", "manifest_count": 2.0}
    )

    assert report.family_name == "partial"
    assert report.missing_metric_count == len(required_cross_family_metrics()) - 1
    assert "mean_heldout_violation" in report.missing_metrics


def test_diagnostic_completeness_table_flattens_missing_metrics() -> None:
    rows = diagnostic_completeness_table(
        [{"family_name": "partial", "manifest_count": 2.0}]
    )

    assert rows[0]["family_name"] == "partial"
    assert "mean_heldout_violation" in str(rows[0]["missing_metrics"])
