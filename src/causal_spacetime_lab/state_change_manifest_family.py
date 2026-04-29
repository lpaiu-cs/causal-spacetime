"""Manifest-family assignment for frozen response-comparison manifests."""

from __future__ import annotations

from dataclasses import dataclass

from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
)

FAMILY_KINDS = {
    "structured",
    "ineligible_control",
    "destructive_null",
    "symmetry_control",
    "marginal_preserving_null",
    "failed_control",
}


@dataclass(frozen=True)
class ManifestFamilySpec:
    """Rule for assigning frozen manifests to a preregistered family."""

    family_name: str
    family_kind: str
    description: str
    include_eligible: bool = True
    include_ineligible: bool = True
    manifest_id_patterns: list[str] | None = None
    comparison_methods: list[str] | None = None
    min_margin_values: list[float] | None = None

    def __post_init__(self) -> None:
        if self.family_kind not in FAMILY_KINDS:
            allowed = ", ".join(sorted(FAMILY_KINDS))
            raise ValueError(f"family_kind must be one of: {allowed}")


@dataclass(frozen=True)
class ManifestFamilyAssignment:
    """Assignment of one frozen manifest to one comparison family."""

    manifest_id: str
    family_name: str
    family_kind: str
    eligible: bool
    failed_reasons: list[str]


def default_manifest_family_specs() -> list[ManifestFamilySpec]:
    """Return preregistered default manifest-family specifications."""

    return [
        ManifestFamilySpec(
            family_name="eligible_rank_gap",
            family_kind="structured",
            description="Eligible rank-gap response-comparison manifests.",
            include_eligible=True,
            include_ineligible=False,
            comparison_methods=["rank_gap_mean", "rank_gap_median"],
        ),
        ManifestFamilySpec(
            family_name="eligible_combined",
            family_kind="structured",
            description="Eligible combined gap-and-reachability manifests.",
            include_eligible=True,
            include_ineligible=False,
            comparison_methods=["combined_gap_and_mismatch"],
        ),
        ManifestFamilySpec(
            family_name="eligible_other",
            family_kind="structured",
            description="Other eligible response-comparison manifests.",
            include_eligible=True,
            include_ineligible=False,
        ),
        ManifestFamilySpec(
            family_name="ineligible_reported",
            family_kind="ineligible_control",
            description="Ineligible manifests that remain in accounting.",
            include_eligible=False,
            include_ineligible=True,
        ),
        ManifestFamilySpec(
            family_name="failed_synthetic_controls",
            family_kind="failed_control",
            description="Synthetic failed controls used for reporting checks.",
            include_eligible=False,
            include_ineligible=True,
            manifest_id_patterns=["failed_control", "ineligible_exact"],
        ),
    ]


def _matches_spec(
    dataset: ManifestConstraintDataset,
    spec: ManifestFamilySpec,
) -> bool:
    if dataset.eligible and not spec.include_eligible:
        return False
    if not dataset.eligible and not spec.include_ineligible:
        return False
    manifest_json = dataset.manifest_json
    method = str(manifest_json.get("comparison_method", ""))
    if spec.comparison_methods is not None and method not in spec.comparison_methods:
        return False
    if spec.min_margin_values is not None:
        margin = float(manifest_json.get("min_margin", float("nan")))
        if not any(abs(margin - value) < 1e-12 for value in spec.min_margin_values):
            return False
    if spec.manifest_id_patterns is not None:
        haystack = " ".join(
            [
                dataset.manifest_id,
                str(manifest_json.get("profile_label", "")),
                " ".join(dataset.failed_reasons),
            ]
        ).lower()
        return any(pattern.lower() in haystack for pattern in spec.manifest_id_patterns)
    return True


def assign_manifest_to_family(
    dataset: ManifestConstraintDataset,
    specs: list[ManifestFamilySpec],
) -> ManifestFamilyAssignment:
    """Assign one frozen manifest to the first matching preregistered family."""

    ordered_specs = sorted(
        specs,
        key=lambda spec: 0 if spec.family_kind == "failed_control" else 1,
    )
    for spec in ordered_specs:
        if _matches_spec(dataset, spec):
            return ManifestFamilyAssignment(
                manifest_id=dataset.manifest_id,
                family_name=spec.family_name,
                family_kind=spec.family_kind,
                eligible=dataset.eligible,
                failed_reasons=list(dataset.failed_reasons),
            )
    fallback = "eligible_other" if dataset.eligible else "ineligible_reported"
    fallback_kind = "structured" if dataset.eligible else "ineligible_control"
    return ManifestFamilyAssignment(
        manifest_id=dataset.manifest_id,
        family_name=fallback,
        family_kind=fallback_kind,
        eligible=dataset.eligible,
        failed_reasons=list(dataset.failed_reasons),
    )


def assign_manifest_families(
    datasets: list[ManifestConstraintDataset],
    specs: list[ManifestFamilySpec],
) -> list[ManifestFamilyAssignment]:
    """Assign all loaded manifests to preregistered families."""

    return [assign_manifest_to_family(dataset, specs) for dataset in datasets]
