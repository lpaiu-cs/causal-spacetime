"""Exact sanity checks for v2 blocking root-cause analysis."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    decompose_v2_blocking_by_family,
)


def _strong_structural_row() -> dict[str, float | str]:
    return {
        "family_name": "structural_low_count",
        "family_kind": "structured",
        "manifest_count": 1.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.05,
        "mean_generalization_gap": 0.01,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.02,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 1.0,
        "restart_std": 0.02,
        "latent_order_disagreement": 0.05,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def _measured_row() -> dict[str, float | str]:
    row = dict(_strong_structural_row())
    row["family_name"] = "measured_high_heldout"
    row["manifest_count"] = 5.0
    row["mean_heldout_violation"] = 0.6
    return row


def _failed_control_row() -> dict[str, float | str]:
    row = dict(_strong_structural_row())
    row["family_name"] = "failed_control_case"
    row["family_kind"] = "failed_control"
    return row


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic blocking classification checks."""

    records = decompose_v2_blocking_by_family(
        [_strong_structural_row(), _measured_row(), _failed_control_row()],
        default_cross_family_robustness_criteria(),
    )
    structural = any(
        record.family_name == "structural_low_count"
        and record.criterion_name == "manifest_count"
        and record.blocking_type == "structural_blocking"
        and record.status == "fail"
        for record in records
    )
    measured = any(
        record.family_name == "measured_high_heldout"
        and record.criterion_name == "mean_heldout_violation"
        and record.blocking_type == "measured_blocking"
        and record.status == "fail"
        for record in records
    )
    control = any(
        record.family_name == "failed_control_case"
        and record.blocking_type == "control_family_blocking"
        for record in records
    )
    return [
        {"check": "structural_manifest_count_detected", "passed": float(structural)},
        {"check": "measured_heldout_detected", "passed": float(measured)},
        {"check": "control_family_detected", "passed": float(control)},
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v2_blocking_analysis_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v2 blocking analysis exact sanity: {path}")


if __name__ == "__main__":
    main()

