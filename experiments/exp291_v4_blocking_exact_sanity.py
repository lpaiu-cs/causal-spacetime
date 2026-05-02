"""Exact sanity checks for v4 blocked-decision decomposition."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_blocked_v5_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_blocking_analysis import (
    decompose_v4_blocking_by_family,
)


def _base_row() -> dict[str, float | str]:
    return {
        "family_name": "passing_v4",
        "family_kind": "structured",
        "manifest_count": 5.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.01,
        "mean_generalization_gap": 0.0,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 1.0,
        "symmetry_control_gap": 0.0,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 1.0,
        "restart_std": 0.0,
        "latent_order_disagreement": 0.0,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
        "preconditions_passed": 1.0,
    }


def run_experiment(output_dir: Path = Path("outputs")) -> list[dict[str, float | str]]:
    criteria = default_cross_family_robustness_criteria()
    rows: list[dict[str, float | str]] = []
    for family_name, updates in [
        ("heldout_fail", {"mean_heldout_violation": 1.0}),
        ("stricter_fail", {"stricter_threshold_pass_fraction": 0.0}),
        ("null_fail", {"destructive_null_gap": -0.1}),
        ("latent_fail", {"latent_order_disagreement": 1.0}),
        ("coverage_fail", {"pair_node_coverage_fraction": 0.0}),
    ]:
        row = _base_row()
        row["family_name"] = family_name
        row.update(updates)
        records = decompose_v4_blocking_by_family([row], [], [], criteria)
        failed_roots = {
            record.root_cause_category
            for record in records
            if record.status == "fail"
        }
        expected = {
            "heldout_fail": "heldout_failure",
            "stricter_fail": "stricter_pass_failure",
            "null_fail": "null_separation_failure",
            "latent_fail": "latent_order_instability",
            "coverage_fail": "coverage_failure",
        }[family_name]
        rows.append(
            {
                "check": family_name,
                "passed": float(expected in failed_roots),
                "root_causes": ";".join(sorted(failed_roots)),
            }
        )
    control = _base_row()
    control["family_name"] = "failed_controls_v4"
    control["family_kind"] = "failed_control"
    control_records = decompose_v4_blocking_by_family([control], [], [], criteria)
    rows.append(
        {
            "check": "failed_control_family",
            "passed": float(
                all(
                    record.blocking_type == "control_family_blocking"
                    for record in control_records
                )
            ),
            "root_causes": ";".join(
                sorted({record.root_cause_category for record in control_records})
            ),
        }
    )
    return rows


def main() -> None:
    output_dir = Path("outputs")
    path = write_csv(
        run_experiment(output_dir),
        data_path(output_dir, "v4_blocking_exact_sanity.csv"),
        ["check", "passed", "root_causes"],
    )
    print(f"Wrote v4 blocking exact sanity: {path}")


if __name__ == "__main__":
    main()
