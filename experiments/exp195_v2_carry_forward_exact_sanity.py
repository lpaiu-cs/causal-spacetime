"""Exact sanity checks for v2 carry-forward decision utilities."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_carry_forward import (
    decide_v2_family_robustness,
    v2_decision_to_row,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _metric_row(
    family_name: str,
    family_kind: str,
    *,
    strong: bool,
) -> dict[str, float | str]:
    if family_kind == "failed_control":
        heldout = 1.0
        return {
            "family_name": family_name,
            "family_kind": family_kind,
            "manifest_count": 1.0,
            "fitted_fraction": 0.0,
            "no_fit_fraction": 1.0,
            "mean_heldout_violation": heldout,
            "mean_generalization_gap": 0.0,
            "stricter_threshold_pass_fraction": 0.0,
            "destructive_null_gap": 0.0,
            "symmetry_control_gap": 0.0,
            "target_coverage_fraction": 0.0,
            "pair_node_coverage_fraction": 0.0,
            "restart_std": 0.0,
            "latent_order_disagreement": 0.0,
            "no_retuning_audit_pass": 1.0,
            "failed_accounting_present": 1.0,
        }
    return {
        "family_name": family_name,
        "family_kind": family_kind,
        "manifest_count": 3.0 if strong else 1.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.05 if strong else 0.5,
        "mean_generalization_gap": 0.01,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.01,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 0.8,
        "restart_std": 0.02,
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic v2 carry-forward sanity checks."""

    criteria = default_cross_family_robustness_criteria()
    decisions = decide_v2_family_robustness(
        [
            _metric_row("strong_v2", "structured", strong=True),
            _metric_row("weak_v2", "structured", strong=False),
            _metric_row("failed_v2", "failed_control", strong=False),
        ],
        criteria,
    )
    rows = [
        v2_decision_to_row(
            decision,
            diagnostic_complete=True,
        )
        for decision in decisions
    ]
    by_family = {row["family_name"]: row for row in rows}
    checks = [
        ("strong_carry_forward", by_family["strong_v2"]["decision"] == "carry_forward"),
        (
            "weak_not_carry_forward",
            by_family["weak_v2"]["decision"] in {"blocked", "provisional"},
        ),
        (
            "failed_control_reported",
            by_family["failed_v2"]["decision"] == "failed_control",
        ),
        ("decision_rows_include_diagnostic_complete", "diagnostic_complete" in rows[0]),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact sanity rows."""

    return write_csv(
        rows,
        output_dir / "data" / "v2_carry_forward_exact_sanity.csv",
        ["check", "passed"],
    )


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote v2 carry-forward exact sanity: {path}")


if __name__ == "__main__":
    main()
