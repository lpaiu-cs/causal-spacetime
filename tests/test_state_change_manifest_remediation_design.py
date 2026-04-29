from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
)
from causal_spacetime_lab.state_change_manifest_remediation_design import (
    default_remediation_proposals,
    propose_remediations_for_failure_records,
)


def test_default_remediation_proposals_preregister_execution_changes() -> None:
    proposals = default_remediation_proposals()
    execution_changing = [
        proposal
        for proposal in proposals
        if proposal.requires_new_preregistration
    ]

    assert execution_changing
    assert all(
        not proposal.allowed_in_current_milestone
        for proposal in execution_changing
    )


def test_propose_remediations_for_failure_records_returns_relevant_proposals() -> None:
    record = CriterionFailureRecord(
        "family",
        "structured",
        "mean_heldout_violation",
        "hard",
        0.4,
        0.2,
        "measured_failure",
        "heldout_failure",
        "failed",
    )

    proposals = propose_remediations_for_failure_records([record])

    assert any(
        proposal.proposal_name == "add_protocol_columns"
        for proposal in proposals
    )
