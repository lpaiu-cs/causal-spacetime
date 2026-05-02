"""Exact sanity checks for v4 protocol carry-forward decisions."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_carry_forward import (
    decide_v4_protocol_family_robustness,
)
from causal_spacetime_lab.state_change_manifest_v4_preconditions import (
    V4ProtocolPreconditionReport,
)


def _metric_row(name: str, kind: str, heldout: float = 0.1) -> dict[str, float | str]:
    return {
        "family_name": name,
        "family_kind": kind,
        "manifest_count": 3.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": heldout,
        "mean_generalization_gap": 0.01,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.02,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 0.8,
        "restart_std": 0.01,
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def _precondition(name: str, passed: bool = True) -> V4ProtocolPreconditionReport:
    failed = [] if passed else ["invalid_provenance"]
    return V4ProtocolPreconditionReport(
        family_name=name,
        family_kind="structured",
        manifest_count=3,
        diagnostic_complete=True,
        all_manifests_have_measurement_protocol=True,
        all_manifests_have_profile_metadata=True,
        all_manifests_have_handoff_provenance=True,
        all_structured_protocol_invariant=True,
        all_structured_parameter_complete=True,
        all_structured_admissible_for_pairwise_dissimilarity=True,
        all_structured_valid_provenance=passed,
        report_only_controls_ineligible=True,
        failed_controls_ineligible=True,
        preconditions_passed=passed,
        failed_preconditions=failed,
        warning_preconditions=[],
    )


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic M45 decision sanity rows."""

    metrics = [
        _metric_row("strong", "structured"),
        _metric_row("weak", "structured", heldout=0.5),
        _metric_row("bad_provenance", "structured"),
        _metric_row("report", "report_only"),
        _metric_row("failed", "failed_control"),
    ]
    decisions = decide_v4_protocol_family_robustness(
        metrics,
        default_cross_family_robustness_criteria(),
        [
            _precondition("strong"),
            _precondition("weak"),
            _precondition("bad_provenance", passed=False),
        ],
    )
    by_family = {decision.family_name: decision for decision in decisions}
    return [
        {
            "check": "strong_carry_forward",
            "passed": float(by_family["strong"].decision == "carry_forward"),
        },
        {
            "check": "weak_not_carry_forward",
            "passed": float(by_family["weak"].decision in {"blocked", "provisional"}),
        },
        {
            "check": "provenance_failed_not_carry_forward",
            "passed": float(by_family["bad_provenance"].decision != "carry_forward"),
        },
        {
            "check": "report_only_stays_report_only",
            "passed": float(by_family["report"].decision == "report_only"),
        },
        {
            "check": "failed_control_stays_failed_control",
            "passed": float(by_family["failed"].decision == "failed_control"),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v4_protocol_carry_forward_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v4 protocol carry-forward exact sanity: {path}")


if __name__ == "__main__":
    main()
