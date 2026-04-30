"""Planned-only v4 protocol remediation design utilities."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class V4ProtocolFamilyDesign:
    """One planned-only v4 protocol family design row."""

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
    linked_v3_root_causes: list[str]

    def __post_init__(self) -> None:
        if self.execution_status != "planned_only":
            raise ValueError("v4 designs must be planned_only in M43")


def _design(
    family_name: str,
    measurement_protocol_family: str,
    profile_resolution_policy: str,
    linked: list[str],
    *,
    family_kind: str = "structured",
    comparison_method: str = "rank_gap_mean",
) -> V4ProtocolFamilyDesign:
    return V4ProtocolFamilyDesign(
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
        target_inclusion_policy="predeclare_target_and_pair_coverage_floor",
        null_taxonomy_policy="separate_destructive_nulls_from_symmetry_controls",
        stability_policy=(
            "require_restart_stability_outputs;"
            "require_latent_order_stability_outputs;"
            "predeclare_optimizer_restarts"
        ),
        comparison_method=comparison_method,
        margin_policy="fixed_margin",
        margin_value=0.05,
        execution_status="planned_only",
        linked_v3_root_causes=linked,
    )


def default_v4_protocol_family_designs() -> list[V4ProtocolFamilyDesign]:
    """Return planned-only v4 protocol family designs."""

    return [
        _design(
            "rank_gap_earliest_full_stability_v4",
            "earliest_full_reference",
            "increase_unique_delay_ranks",
            ["heldout_failure", "latent_order_instability"],
        ),
        _design(
            "rank_gap_earliest_retained_resolution_v4",
            "earliest_retained_reference",
            "reduce_tie_fraction",
            ["heldout_failure", "restart_instability"],
        ),
        _design(
            "rank_gap_gated_full_stability_v4",
            "gated_full_reference",
            "increase_reference_set_diversity",
            ["heldout_failure", "latent_order_instability"],
        ),
        _design(
            "combined_earliest_full_stability_v4",
            "earliest_full_reference",
            "fixed_protocol_columns_only",
            ["generalization_gap_failure", "null_separation_failure"],
            comparison_method="combined_gap_and_mismatch",
        ),
        _design(
            "rank_gap_high_resolution_reference_set_v4",
            "earliest_full_reference",
            "increase_unique_delay_ranks;increase_reference_set_diversity",
            ["coverage_failure", "latent_order_instability"],
        ),
        _design(
            "failed_controls_v4",
            "earliest_full_reference",
            "fixed_protocol_columns_only",
            ["control_family_blocking"],
            family_kind="failed_control",
        ),
        _design(
            "report_only_immediate_edge_v4",
            "immediate_edge_reference",
            "fixed_protocol_columns_only",
            ["control_family_blocking"],
            family_kind="report_only",
        ),
    ]


def v4_protocol_family_design_table(
    designs: list[V4ProtocolFamilyDesign],
) -> list[dict[str, float | str]]:
    """Convert planned v4 designs to CSV-safe rows."""

    rows: list[dict[str, float | str]] = []
    for design in designs:
        row = asdict(design)
        row["planned_manifest_count"] = float(design.planned_manifest_count)
        row["linked_v3_root_causes"] = ";".join(design.linked_v3_root_causes)
        rows.append(row)  # type: ignore[arg-type]
    return rows

