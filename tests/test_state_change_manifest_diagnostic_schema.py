from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
    diagnostic_metric_requirement_table,
    missing_metric_remediation_priority,
)


def test_default_diagnostic_metric_requirements_returns_14_metrics() -> None:
    requirements = default_diagnostic_metric_requirements()

    assert len(requirements) == 14


def test_diagnostic_metric_requirement_table_includes_target_coverage() -> None:
    rows = diagnostic_metric_requirement_table()

    assert any(row["metric_name"] == "target_coverage_fraction" for row in rows)


def test_missing_metric_remediation_priority_hard_metric_high() -> None:
    assert missing_metric_remediation_priority("mean_heldout_violation") == "high"
