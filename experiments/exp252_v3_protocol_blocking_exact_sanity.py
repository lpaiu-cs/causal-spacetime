"""Exact sanity checks for v3 protocol blocking decomposition."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    decompose_v3_protocol_blocking_by_family,
)


def _metric_row(
    family_name: str,
    family_kind: str = "structured",
    *,
    heldout: float = 0.05,
    latent: float = 0.1,
) -> dict[str, float | str]:
    return {
        "family_name": family_name,
        "family_kind": family_kind,
        "manifest_count": 3.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": heldout,
        "mean_generalization_gap": 0.01,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.01,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 0.8,
        "restart_std": 0.01,
        "latent_order_disagreement": latent,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
        "preconditions_passed": 1.0,
    }


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic blocking sanity checks."""

    rows = [
        _metric_row("passing"),
        _metric_row("heldout_bad", heldout=0.5),
        _metric_row("latent_bad", latent=0.8),
        _metric_row("failed", family_kind="failed_control"),
    ]
    records = decompose_v3_protocol_blocking_by_family(
        rows,
        [],
        [],
        default_cross_family_robustness_criteria(),
    )
    by_family = {}
    for record in records:
        by_family.setdefault(record.family_name, []).append(record)
    return [
        {
            "check": "heldout_failure_is_measured_blocking",
            "passed": float(
                any(
                    record.root_cause_category == "heldout_failure"
                    and record.blocking_type == "measured_blocking"
                    and record.status == "fail"
                    for record in by_family["heldout_bad"]
                )
            ),
        },
        {
            "check": "latent_order_failure_is_measured_blocking",
            "passed": float(
                any(
                    record.root_cause_category == "latent_order_instability"
                    and record.blocking_type == "measured_blocking"
                    and record.status == "fail"
                    for record in by_family["latent_bad"]
                )
            ),
        },
        {
            "check": "failed_control_is_control_family_blocking",
            "passed": float(
                all(
                    record.blocking_type == "control_family_blocking"
                    for record in by_family["failed"]
                )
            ),
        },
        {
            "check": "passing_family_has_no_blocking_failures",
            "passed": float(
                not any(record.status == "fail" for record in by_family["passing"])
            ),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v3_protocol_blocking_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v3 protocol blocking exact sanity: {path}")


if __name__ == "__main__":
    main()

