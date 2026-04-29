"""Exact sanity checks for cross-family robustness decisions."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    decide_family_robustness,
    default_cross_family_robustness_criteria,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _strong_family() -> dict[str, float | str]:
    return {
        "family_name": "strong_structured",
        "family_kind": "structured",
        "manifest_count": 4.0,
        "fitted_fraction": 0.9,
        "no_fit_fraction": 0.1,
        "mean_heldout_violation": 0.08,
        "mean_generalization_gap": 0.02,
        "stricter_threshold_pass_fraction": 0.75,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.02,
        "target_coverage_fraction": 0.9,
        "pair_node_coverage_fraction": 0.7,
        "restart_std": 0.03,
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run deterministic robustness decision checks."""

    criteria = default_cross_family_robustness_criteria()
    weak = dict(_strong_family())
    weak["family_name"] = "weak_structured"
    weak["mean_heldout_violation"] = 0.35
    ineligible = dict(_strong_family())
    ineligible["family_name"] = "ineligible_family"
    ineligible["family_kind"] = "ineligible_control"
    decisions = [
        decide_family_robustness(_strong_family(), criteria),
        decide_family_robustness(weak, criteria),
        decide_family_robustness(ineligible, criteria),
    ]
    checks = [
        ("strong_carry_forward", decisions[0].decision == "carry_forward"),
        (
            "weak_not_carried_forward",
            decisions[1].decision in {"blocked", "provisional"},
        ),
        ("ineligible_report_only", decisions[2].decision == "report_only"),
    ]
    return [
        {
            "check": name,
            "passed": float(passed),
            "decision": decisions[index].decision,
        }
        for index, (name, passed) in enumerate(checks)
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact sanity CSV."""

    path = output_dir / "data" / "cross_family_robustness_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed", "decision"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    path = write_outputs(rows)
    print(f"Wrote cross-family robustness exact sanity: {path}")


if __name__ == "__main__":
    main()
