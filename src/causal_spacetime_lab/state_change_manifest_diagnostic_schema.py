"""Diagnostic-complete metric requirements for manifest family outputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass

CRITERION_TYPES = {"hard", "non_hard", "warning", "accounting"}
MISSING_STATUSES = {
    "blocking_missing_metric",
    "warning_missing_metric",
    "not_applicable",
}


@dataclass(frozen=True)
class DiagnosticMetricRequirement:
    """One required family-level diagnostic metric."""

    metric_name: str
    criterion_type: str
    required_for_carry_forward: bool
    source_output_file: str
    source_milestone: str
    missing_status: str
    description: str

    def __post_init__(self) -> None:
        if self.criterion_type not in CRITERION_TYPES:
            allowed = ", ".join(sorted(CRITERION_TYPES))
            raise ValueError(f"criterion_type must be one of: {allowed}")
        if self.missing_status not in MISSING_STATUSES:
            allowed = ", ".join(sorted(MISSING_STATUSES))
            raise ValueError(f"missing_status must be one of: {allowed}")


def _missing_status_for_type(criterion_type: str) -> str:
    if criterion_type == "hard":
        return "blocking_missing_metric"
    return "warning_missing_metric"


def _requirement(
    metric_name: str,
    criterion_type: str,
    source_output_file: str,
    source_milestone: str,
    description: str,
) -> DiagnosticMetricRequirement:
    return DiagnosticMetricRequirement(
        metric_name=metric_name,
        criterion_type=criterion_type,
        required_for_carry_forward=True,
        source_output_file=source_output_file,
        source_milestone=source_milestone,
        missing_status=_missing_status_for_type(criterion_type),
        description=description,
    )


def default_diagnostic_metric_requirements() -> list[DiagnosticMetricRequirement]:
    """Return the diagnostic-complete cross-family metric requirements."""

    fit_summary = "outputs/data/manifest_family_fit_summary.csv"
    null_taxonomy = "outputs/data/manifest_family_null_taxonomy.csv"
    stricter = "outputs/data/manifest_family_stricter_criteria.csv"
    failed_accounting = "outputs/data/manifest_family_failed_manifest_accounting.csv"
    no_retuning = "outputs/data/manifest_family_no_retuning_audit.csv"
    coverage = "outputs/data/response_constraint_pool_coverage.csv"
    stability = "outputs/data/frozen_manifest_fit_stability.csv"
    return [
        _requirement(
            "manifest_count",
            "hard",
            fit_summary,
            "Milestone 33",
            "Number of manifests assigned to the family.",
        ),
        _requirement(
            "fitted_fraction",
            "non_hard",
            fit_summary,
            "Milestone 33",
            "Fraction of family manifests that produced fit rows.",
        ),
        _requirement(
            "no_fit_fraction",
            "non_hard",
            fit_summary,
            "Milestone 33",
            "Fraction of family manifests retained as no-fit rows.",
        ),
        _requirement(
            "mean_heldout_violation",
            "hard",
            fit_summary,
            "Milestone 33",
            "Mean held-out response-comparison violation rate.",
        ),
        _requirement(
            "mean_generalization_gap",
            "non_hard",
            fit_summary,
            "Milestone 33",
            "Mean held-out minus train violation rate.",
        ),
        _requirement(
            "stricter_threshold_pass_fraction",
            "non_hard",
            stricter,
            "Milestone 33",
            "Fraction passing fixed stricter family-level diagnostics.",
        ),
        _requirement(
            "destructive_null_gap",
            "non_hard",
            null_taxonomy,
            "Milestone 33",
            "Destructive-null held-out violation minus structured violation.",
        ),
        _requirement(
            "symmetry_control_gap",
            "non_hard",
            null_taxonomy,
            "Milestone 33",
            "Absolute gap between symmetry-control and structured violation.",
        ),
        _requirement(
            "target_coverage_fraction",
            "warning",
            coverage,
            "Milestone 30 aggregation",
            "Fraction of targets touched by the constraint pool.",
        ),
        _requirement(
            "pair_node_coverage_fraction",
            "warning",
            coverage,
            "Milestone 30 aggregation",
            "Fraction of unordered target-pair nodes touched by constraints.",
        ),
        _requirement(
            "restart_std",
            "warning",
            stability,
            "Milestone 32 aggregation",
            "Standard deviation of held-out violation across restarts.",
        ),
        _requirement(
            "latent_order_disagreement",
            "warning",
            stability,
            "Milestone 32 aggregation",
            "Mean disagreement among latent-order diagnostics across restarts.",
        ),
        _requirement(
            "no_retuning_audit_pass",
            "accounting",
            no_retuning,
            "Milestone 33",
            "Indicator that fixed settings were not changed after fit results.",
        ),
        _requirement(
            "failed_accounting_present",
            "accounting",
            failed_accounting,
            "Milestone 33",
            "Indicator that failed and ineligible families remain reported.",
        ),
    ]


def diagnostic_metric_requirement_table() -> list[dict[str, float | str]]:
    """Return diagnostic requirements as flat table rows."""

    rows: list[dict[str, float | str]] = []
    for requirement in default_diagnostic_metric_requirements():
        row = asdict(requirement)
        row["required_for_carry_forward"] = float(
            requirement.required_for_carry_forward
        )
        rows.append(row)
    return rows


def missing_metric_remediation_priority(metric_name: str) -> str:
    """Return remediation priority for a missing diagnostic metric."""

    requirements = {
        requirement.metric_name: requirement
        for requirement in default_diagnostic_metric_requirements()
    }
    criterion_type = requirements.get(
        metric_name,
        DiagnosticMetricRequirement(
            metric_name=metric_name,
            criterion_type="warning",
            required_for_carry_forward=False,
            source_output_file="",
            source_milestone="",
            missing_status="warning_missing_metric",
            description="Unknown metric.",
        ),
    ).criterion_type
    if criterion_type == "hard":
        return "high"
    if criterion_type == "non_hard":
        return "medium"
    return "low"
