from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_carry_forward import (
    decide_v4_protocol_family_robustness,
    v4_protocol_decision_to_row,
    v4_protocol_metrics_rows_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_v4_preconditions import (
    V4ProtocolPreconditionReport,
)


def _row(
    heldout: float = 0.05,
    *,
    family_kind: str = "structured",
) -> dict[str, float | str]:
    return {
        "family_name": "family",
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
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def _precondition(passed: bool = True) -> V4ProtocolPreconditionReport:
    return V4ProtocolPreconditionReport(
        family_name="family",
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
        failed_preconditions=[] if passed else ["invalid_provenance"],
        warning_preconditions=[],
    )


def test_v4_protocol_metrics_rows_from_bundle_reads_metric_rows() -> None:
    rows = v4_protocol_metrics_rows_from_bundle({"metrics": [{"manifest_count": "3"}]})

    assert rows == [{"manifest_count": 3.0}]


def test_decide_v4_protocol_family_robustness_strong_row() -> None:
    decisions = decide_v4_protocol_family_robustness(
        [_row()],
        default_cross_family_robustness_criteria(),
        [_precondition()],
    )

    assert decisions[0].decision == "carry_forward"


def test_decide_v4_protocol_family_robustness_weak_row() -> None:
    decisions = decide_v4_protocol_family_robustness(
        [_row(heldout=0.5)],
        default_cross_family_robustness_criteria(),
        [_precondition()],
    )

    assert decisions[0].decision in {"blocked", "provisional"}


def test_decide_v4_protocol_family_robustness_precondition_failed() -> None:
    decisions = decide_v4_protocol_family_robustness(
        [_row()],
        default_cross_family_robustness_criteria(),
        [_precondition(passed=False)],
    )

    assert decisions[0].decision != "carry_forward"


def test_decide_v4_protocol_family_robustness_controls_remain_controls() -> None:
    decisions = decide_v4_protocol_family_robustness(
        [
            {**_row(family_kind="report_only"), "family_name": "report"},
            {**_row(family_kind="failed_control"), "family_name": "failed"},
        ],
        default_cross_family_robustness_criteria(),
        [],
    )
    by_family = {decision.family_name: decision for decision in decisions}

    assert by_family["report"].decision == "report_only"
    assert by_family["failed"].decision == "failed_control"


def test_v4_protocol_decision_to_row_includes_precondition_flags() -> None:
    decision = decide_v4_protocol_family_robustness(
        [_row()],
        default_cross_family_robustness_criteria(),
        [_precondition()],
    )[0]

    row = v4_protocol_decision_to_row(
        decision,
        diagnostic_complete=True,
        preconditions_passed=True,
        failed_preconditions=[],
    )

    assert row["diagnostic_complete"] == 1.0
    assert row["preconditions_passed"] == 1.0
