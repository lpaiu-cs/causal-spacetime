"""Planned-only v5 protocol remediation design utilities."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class V5ProtocolFamilyDesign:
    """One planned-only v5 protocol family design row."""

    family_name: str
    family_kind: str
    planned_manifest_count: int
    measurement_protocol_family: str
    handoff_provenance_type: str
    profile_resolution_policy: str
    target_inclusion_policy: str
    null_taxonomy_policy: str
    stability_policy: str
    comparison_method: str
    margin_policy: str
    margin_value: float
    execution_status: str
    linked_v4_root_causes: list[str]
    remediation_rationale: str
    iteration_risk_category: str

    def __post_init__(self) -> None:
        if self.execution_status != "planned_only":
            raise ValueError("v5 designs must be planned_only in M46")
        if self.family_kind == "structured" and self.planned_manifest_count < 3:
            raise ValueError("structured v5 designs need at least 3 manifests")
        if not self.linked_v4_root_causes:
            raise ValueError("v5 designs must link to v4 root causes")


def _design(
    family_name: str,
    measurement_protocol_family: str,
    profile_resolution_policy: str,
    linked: list[str],
    rationale: str,
    risk_category: str,
    *,
    family_kind: str = "structured",
    comparison_method: str = "rank_gap_mean",
) -> V5ProtocolFamilyDesign:
    return V5ProtocolFamilyDesign(
        family_name=family_name,
        family_kind=family_kind,
        planned_manifest_count=5 if family_kind == "structured" else 3,
        measurement_protocol_family=measurement_protocol_family,
        handoff_provenance_type=(
            "hybrid_template_instantiated_from_profile"
            if family_kind == "structured"
            else "report_only_control"
        ),
        profile_resolution_policy=profile_resolution_policy,
        target_inclusion_policy=(
            "predeclare_target_and_pair_coverage_floor;"
            "no_post_hoc_hard_case_removal"
        ),
        null_taxonomy_policy=(
            "predeclare_destructive_nulls;"
            "separate_symmetry_controls;"
            "no_post_hoc_null_relabeling"
        ),
        stability_policy=(
            "predeclare_optimizer_restarts;"
            "require_manifest_transfer_split;"
            "require_latent_order_stability_outputs"
        ),
        comparison_method=comparison_method,
        margin_policy="fixed_margin",
        margin_value=0.05,
        execution_status="planned_only",
        linked_v4_root_causes=linked,
        remediation_rationale=rationale,
        iteration_risk_category=risk_category,
    )


def default_v5_protocol_family_designs() -> list[V5ProtocolFamilyDesign]:
    """Return planned-only v5 protocol family designs."""

    return [
        _design(
            "rank_gap_earliest_full_manifest_transfer_v5",
            "earliest_full_reference",
            "fixed_protocol_columns_only;manifest_transfer_split",
            ["heldout_failure", "latent_order_instability"],
            (
                "Test held-out stability through preregistered manifest-transfer "
                "splits under the same fixed criteria."
            ),
            "stability_design",
        ),
        _design(
            "rank_gap_high_coverage_strict_pair_v5",
            "earliest_full_reference",
            "increase_strict_pair_resolution;coverage_floor_before_fit",
            ["coverage_failure", "heldout_failure"],
            (
                "Target pair-node coverage and strict-pair resolution without "
                "removing hard cases after evaluation."
            ),
            "coverage_design",
        ),
        _design(
            "rank_gap_latent_stability_replicated_v5",
            "gated_full_reference",
            "low_tie_reference_grid;replicated_reference_sets",
            ["latent_order_instability", "heldout_failure"],
            (
                "Assess whether preregistered reference replication improves "
                "latent-order stability."
            ),
            "stability_design",
        ),
        _design(
            "combined_null_separation_v5",
            "earliest_full_reference",
            "fixed_protocol_columns_only",
            ["null_separation_failure", "stricter_pass_failure"],
            (
                "Separate destructive-null evaluation from symmetry controls "
                "with fixed taxonomy and fixed thresholds."
            ),
            "justified_protocol_design",
            comparison_method="combined_gap_and_mismatch",
        ),
        _design(
            "rank_gap_low_tie_reference_diverse_v5",
            "earliest_retained_reference",
            "reduce_tie_fraction;increase_reference_set_diversity",
            ["latent_order_instability", "stricter_pass_failure"],
            (
                "Increase preregistered reference diversity to reduce tied "
                "profile comparisons."
            ),
            "profile_resolution",
        ),
        _design(
            "failed_controls_v5",
            "earliest_full_reference",
            "fixed_protocol_columns_only",
            ["control_family_blocking"],
            "Preserve failed-control accounting as ineligible controls.",
            "diagnostic_completeness",
            family_kind="failed_control",
        ),
        _design(
            "report_only_immediate_edge_v5",
            "immediate_edge_reference",
            "fixed_protocol_columns_only",
            ["control_family_blocking"],
            "Preserve immediate-edge report-only controls as ineligible rows.",
            "diagnostic_completeness",
            family_kind="report_only",
        ),
    ]


def v5_protocol_family_design_table(
    designs: list[V5ProtocolFamilyDesign],
) -> list[dict[str, float | str]]:
    """Convert planned v5 designs to CSV-safe rows."""

    rows: list[dict[str, float | str]] = []
    for design in designs:
        row = asdict(design)
        row["planned_manifest_count"] = float(design.planned_manifest_count)
        row["linked_v4_root_causes"] = ";".join(design.linked_v4_root_causes)
        rows.append(row)  # type: ignore[arg-type]
    return rows
