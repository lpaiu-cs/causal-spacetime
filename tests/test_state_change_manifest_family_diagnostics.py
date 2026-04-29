from __future__ import annotations

import pytest

from causal_spacetime_lab.state_change_manifest_family import ManifestFamilyAssignment
from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    failed_manifest_accounting_summary,
    manifest_fit_diagnostic_row,
    summarize_family_fit_diagnostics,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationFit,
)


def _fit() -> ManifestRepresentationFit:
    return ManifestRepresentationFit(
        manifest_id="m1",
        embedding_dim=1,
        eligible=True,
        fitted=True,
        reason_not_fit="",
        train_violation_rate=0.1,
        heldout_violation_rate=0.25,
        train_hinge_loss=0.01,
        heldout_hinge_loss=0.02,
        train_constraint_count=10,
        heldout_constraint_count=5,
        target_count=4,
        embedding=None,
    )


def _assignment() -> ManifestFamilyAssignment:
    return ManifestFamilyAssignment(
        manifest_id="m1",
        family_name="eligible_rank_gap",
        family_kind="structured",
        eligible=True,
        failed_reasons=[],
    )


def test_manifest_fit_diagnostic_row_computes_generalization_gap() -> None:
    row = manifest_fit_diagnostic_row(_fit(), _assignment())

    assert row.generalization_gap == pytest.approx(0.15)


def test_summarize_family_fit_diagnostics_groups_rows() -> None:
    row = manifest_fit_diagnostic_row(_fit(), _assignment())

    summary = summarize_family_fit_diagnostics([row])

    assert summary[0]["family_name"] == "eligible_rank_gap"
    assert summary[0]["fitted_count"] == 1.0


def test_failed_manifest_accounting_summary_counts_reasons() -> None:
    rows = failed_manifest_accounting_summary(
        [
            _assignment(),
            ManifestFamilyAssignment(
                manifest_id="m2",
                family_name="ineligible_reported",
                family_kind="ineligible_control",
                eligible=False,
                failed_reasons=["low_signal"],
            ),
        ]
    )

    assert any(row["reason"] == "low_signal" for row in rows)
