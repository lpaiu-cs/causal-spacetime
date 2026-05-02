from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_remediation_iteration_risk import (
    audit_v5_remediation_iteration_risk,
)


def test_remediation_iteration_risk_flags_blocking_changes() -> None:
    records = audit_v5_remediation_iteration_risk(
        [{"root_cause_category": "heldout_failure", "status": "fail"}],
        [],
        [
            {
                "family_name": "threshold_change",
                "linked_v4_root_causes": "heldout_failure",
                "profile_resolution_policy": "threshold adjustment",
                "remediation_rationale": "",
            },
            {
                "family_name": "try more families until one passes",
                "linked_v4_root_causes": "heldout_failure",
                "profile_resolution_policy": "fixed_protocol_columns_only",
                "remediation_rationale": "",
            },
        ],
    )

    assert [record.risk_level for record in records] == ["blocking", "blocking"]


def test_remediation_iteration_risk_allows_no_post_hoc_hard_case_removal_policy(
) -> None:
    records = audit_v5_remediation_iteration_risk(
        [{"root_cause_category": "coverage_failure", "status": "fail"}],
        [],
        [
            {
                "family_name": "coverage_design",
                "linked_v4_root_causes": "coverage_failure",
                "target_inclusion_policy": (
                    "predeclare_target_and_pair_coverage_floor;"
                    "no_post_hoc_hard_case_removal"
                ),
                "profile_resolution_policy": "increase_strict_pair_resolution",
                "remediation_rationale": (
                    "Target pair-node coverage without removing hard cases."
                ),
            },
        ],
    )

    assert records[0].allowed_as_v5_design
    assert records[0].risk_category == "coverage_design"


def test_remediation_iteration_risk_honors_declared_nonblocking_category() -> None:
    records = audit_v5_remediation_iteration_risk(
        [{"root_cause_category": "coverage_failure", "status": "fail"}],
        [],
        [
            {
                "family_name": "stability_design",
                "linked_v4_root_causes": "latent_order_instability",
                "target_inclusion_policy": "predeclare_target_coverage_floor",
                "iteration_risk_category": "stability_design",
                "remediation_rationale": "Predeclare restart outputs.",
            },
        ],
    )

    assert records[0].allowed_as_v5_design
    assert records[0].risk_category == "stability_design"
