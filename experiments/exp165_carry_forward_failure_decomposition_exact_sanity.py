"""Exact sanity checks for carry-forward failure decomposition."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    decompose_family_metric_failures,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _strong_family() -> dict[str, float | str]:
    return {
        "family_name": "strong_structured",
        "family_kind": "structured",
        "manifest_count": 4.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.08,
        "mean_generalization_gap": 0.03,
        "stricter_threshold_pass_fraction": 0.75,
        "destructive_null_gap": 0.22,
        "symmetry_control_gap": 0.03,
        "target_coverage_fraction": 0.9,
        "pair_node_coverage_fraction": 0.7,
        "restart_std": 0.03,
        "latent_order_disagreement": 0.12,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic decomposition checks."""

    criteria = default_cross_family_robustness_criteria()
    heldout_failure = dict(_strong_family())
    heldout_failure["family_name"] = "heldout_failure_family"
    heldout_failure["mean_heldout_violation"] = 0.35
    missing_coverage = dict(_strong_family())
    missing_coverage["family_name"] = "missing_coverage_family"
    missing_coverage.pop("target_coverage_fraction")

    strong_records = decompose_family_metric_failures(_strong_family(), criteria)
    heldout_records = decompose_family_metric_failures(heldout_failure, criteria)
    missing_records = decompose_family_metric_failures(missing_coverage, criteria)

    heldout_status = next(
        record.status
        for record in heldout_records
        if record.criterion_name == "mean_heldout_violation"
    )
    coverage_status = next(
        record.status
        for record in missing_records
        if record.criterion_name == "target_coverage_fraction"
    )
    strong_hard_pass = all(
        record.status == "pass"
        for record in strong_records
        if record.criterion_type == "hard"
    )
    checks = [
        ("measured_heldout_failure", heldout_status == "measured_failure"),
        ("missing_coverage_is_missing_metric", coverage_status == "missing_metric"),
        ("strong_family_hard_criteria_pass", strong_hard_pass),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact sanity CSV."""

    path = output_dir / "data" / "carry_forward_failure_decomposition_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote carry-forward failure decomposition exact sanity: {path}")


if __name__ == "__main__":
    main()
