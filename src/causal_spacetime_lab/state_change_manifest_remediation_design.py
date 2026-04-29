"""Report-only upstream remediation design proposals."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
)


@dataclass(frozen=True)
class RemediationProposal:
    """A report-only upstream remediation design proposal."""

    proposal_name: str
    target_root_cause: str
    description: str
    expected_effect: str
    requires_new_preregistration: bool
    allowed_in_current_milestone: bool


def default_remediation_proposals() -> list[RemediationProposal]:
    """Return default report-only remediation proposals."""

    return [
        RemediationProposal(
            "add_protocol_columns",
            "heldout_failure",
            "Design richer response profiles with additional protocol columns.",
            "May improve held-out agreement in new preregistered manifests.",
            True,
            False,
        ),
        RemediationProposal(
            "increase_target_count",
            "coverage_failure",
            "Design future manifests with more eligible target events.",
            "May improve target and pair-node coverage.",
            True,
            False,
        ),
        RemediationProposal(
            "improve_reachability_filter",
            "coverage_failure",
            "Predeclare stricter reachability filters before manifest creation.",
            "May reduce sparse or unevaluable target-pair coverage.",
            True,
            False,
        ),
        RemediationProposal(
            "increase_unique_rank_resolution",
            "generalization_gap_failure",
            "Design profiles with richer response-rank variation.",
            "May reduce tie-driven instability in future manifests.",
            True,
            False,
        ),
        RemediationProposal(
            "compute_missing_coverage_metrics",
            "missing_metric",
            "Add reporting for target and pair-node coverage metrics.",
            "Improves diagnostic completeness without changing current decisions.",
            False,
            True,
        ),
        RemediationProposal(
            "compute_restart_stability_outputs",
            "restart_instability",
            "Report restart-stability diagnostics at family level.",
            "Separates missing restart diagnostics from measured instability.",
            False,
            True,
        ),
        RemediationProposal(
            "compute_latent_order_stability_outputs",
            "latent_order_instability",
            "Report latent-order stability diagnostics at family level.",
            "Separates missing latent-order diagnostics from measured instability.",
            False,
            True,
        ),
        RemediationProposal(
            "strengthen_null_taxonomy_reporting",
            "null_separation_failure",
            "Report destructive, symmetry-control, and marginal-preserving nulls.",
            "Improves future null-taxonomy interpretability.",
            False,
            True,
        ),
        RemediationProposal(
            "reduce_constraint_sampling_bias",
            "unknown",
            "Predeclare balanced constraint sampling for future manifests.",
            "May reduce family-level diagnostic sensitivity to sampling skew.",
            True,
            False,
        ),
        RemediationProposal(
            "generate_new_preregistered_manifest_family",
            "missing_output",
            "Create a new family only through a fresh preregistered manifest plan.",
            "May supply complete outputs for future carry-forward decisions.",
            True,
            False,
        ),
    ]


def propose_remediations_for_failure_records(
    records: list[CriterionFailureRecord],
) -> list[RemediationProposal]:
    """Return proposals matching observed failure root causes."""

    root_causes = {
        record.root_cause_category
        for record in records
        if record.status in {"measured_failure", "missing_metric", "warning_only"}
    }
    proposals = default_remediation_proposals()
    proposal_names: set[str] = {
        proposal.proposal_name
        for proposal in proposals
        if proposal.target_root_cause in root_causes
    }
    missing_criteria = {
        record.criterion_name
        for record in records
        if record.status == "missing_metric"
    }
    if {
        "target_coverage_fraction",
        "pair_node_coverage_fraction",
    } & missing_criteria:
        proposal_names.add("compute_missing_coverage_metrics")
    if "restart_std" in missing_criteria:
        proposal_names.add("compute_restart_stability_outputs")
    if "latent_order_disagreement" in missing_criteria:
        proposal_names.add("compute_latent_order_stability_outputs")
    selected = [
        proposal for proposal in proposals if proposal.proposal_name in proposal_names
    ]
    if not selected and root_causes:
        selected = [
            proposal
            for proposal in proposals
            if proposal.target_root_cause == "unknown"
        ]
    return selected


def remediation_proposal_table(
    proposals: list[RemediationProposal],
) -> list[dict[str, float | str]]:
    """Convert remediation proposals to CSV rows."""

    rows: list[dict[str, float | str]] = []
    for proposal in proposals:
        row = asdict(proposal)
        row["requires_new_preregistration"] = float(
            proposal.requires_new_preregistration
        )
        row["allowed_in_current_milestone"] = float(
            proposal.allowed_in_current_milestone
        )
        rows.append(row)
    return rows
