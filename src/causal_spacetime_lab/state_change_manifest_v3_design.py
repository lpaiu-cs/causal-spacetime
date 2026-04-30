"""Preregistered v3 manifest-family design tables."""

from __future__ import annotations

from dataclasses import asdict, dataclass

EXECUTION_STATUSES = {"planned_only", "executed"}


@dataclass(frozen=True)
class V3ManifestFamilyDesign:
    """Planned v3 manifest-family design, not an executed result."""

    family_name: str
    family_kind: str
    planned_manifest_count: int
    profile_column_policy: str
    target_inclusion_policy: str
    rank_resolution_policy: str
    coverage_policy: str
    null_taxonomy_policy: str
    restart_stability_policy: str
    latent_order_stability_policy: str
    comparison_method: str
    min_margin: float
    execution_status: str

    def __post_init__(self) -> None:
        if self.execution_status not in EXECUTION_STATUSES:
            allowed = ", ".join(sorted(EXECUTION_STATUSES))
            raise ValueError(f"execution_status must be one of: {allowed}")


def default_v3_manifest_family_designs() -> list[V3ManifestFamilyDesign]:
    """Return preregistered planned-only v3 family designs."""

    common = {
        "family_kind": "structured",
        "planned_manifest_count": 5,
        "target_inclusion_policy": "predeclare target inclusion before generation",
        "coverage_policy": "produce target and pair-node coverage per family",
        "null_taxonomy_policy": "report destructive symmetry and marginal nulls",
        "restart_stability_policy": "produce restart_std metric",
        "latent_order_stability_policy": "produce latent_order_disagreement metric",
        "min_margin": 0.05,
        "execution_status": "planned_only",
    }
    coverage_common = {
        key: value
        for key, value in common.items()
        if key != "target_inclusion_policy"
    }
    return [
        V3ManifestFamilyDesign(
            family_name="rank_gap_multi_manifest_v3",
            profile_column_policy="multi-seed replicated rank-gap profiles",
            rank_resolution_policy="same as v2 baseline",
            comparison_method="rank_gap_mean",
            **common,
        ),
        V3ManifestFamilyDesign(
            family_name="rank_gap_more_protocol_columns_v3",
            profile_column_policy="more protocol columns than v2",
            rank_resolution_policy="same as v2 baseline",
            comparison_method="rank_gap_mean",
            **common,
        ),
        V3ManifestFamilyDesign(
            family_name="rank_gap_rank_resolution_enriched_v3",
            profile_column_policy="multi-seed replicated rank-gap profiles",
            rank_resolution_policy="wider and more delay-rank levels",
            comparison_method="rank_gap_mean",
            **common,
        ),
        V3ManifestFamilyDesign(
            family_name="rank_gap_coverage_enriched_v3",
            profile_column_policy="multi-seed replicated rank-gap profiles",
            rank_resolution_policy="same as v2 baseline",
            target_inclusion_policy="predeclare enriched reachable target inclusion",
            comparison_method="rank_gap_mean",
            **coverage_common,
        ),
        V3ManifestFamilyDesign(
            family_name="combined_diagnostic_complete_v3",
            profile_column_policy="multi-seed replicated combined profiles",
            rank_resolution_policy="same as v2 combined baseline",
            comparison_method="combined_gap_and_mismatch",
            **common,
        ),
        V3ManifestFamilyDesign(
            family_name="failed_controls_v3",
            family_kind="failed_control",
            planned_manifest_count=3,
            profile_column_policy="report-only failed controls",
            target_inclusion_policy="intentionally failed control inclusion",
            rank_resolution_policy="control-only",
            coverage_policy="control coverage reported separately",
            null_taxonomy_policy="control null taxonomy reported separately",
            restart_stability_policy="control no-fit rows",
            latent_order_stability_policy="control no-fit rows",
            comparison_method="rank_gap_mean",
            min_margin=0.05,
            execution_status="planned_only",
        ),
    ]


def v3_manifest_family_design_table(
    designs: list[V3ManifestFamilyDesign],
) -> list[dict[str, float | str]]:
    """Convert v3 family designs to CSV rows."""

    return [asdict(design) for design in designs]
