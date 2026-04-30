"""V2 handoff manifest generation pipeline."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v2_profiles import (
    V2ProfileGenerationConfig,
    build_v2_response_profile,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import V2ManifestFamilySpec
from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    HandoffDecision,
    ResponseConstraintHandoffManifest,
    manifest_digest,
    manifest_to_jsonable,
    write_handoff_manifest,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)


@dataclass(frozen=True)
class V2ManifestGenerationConfig:
    """Configuration for v2 manifest generation."""

    max_constraints: int = 5000
    min_margin: float = 0.05
    train_fraction: float = 0.7
    bootstrap_count: int = 20
    null_repetitions: int = 5
    constraint_seed: int = 0


def comparison_protocol_for_v2_family(
    family_name: str,
) -> PairwiseResponseComparisonProtocol:
    """Return the predeclared pairwise protocol for one v2 family."""

    if family_name.startswith("combined_"):
        return PairwiseResponseComparisonProtocol(
            name=f"{family_name}_combined_gap_and_mismatch",
            method="combined_gap_and_mismatch",
            missing_policy="penalize_mismatch",
            min_common_protocols=1,
            reachability_mismatch_penalty=0.5,
        )
    return PairwiseResponseComparisonProtocol(
        name=f"{family_name}_rank_gap_mean",
        method="rank_gap_mean",
        missing_policy="common_reachable",
        min_common_protocols=1,
    )


def _with_manifest_id(
    manifest: ResponseConstraintHandoffManifest,
) -> ResponseConstraintHandoffManifest:
    digest = manifest_digest(manifest_to_jsonable(manifest))
    return replace(manifest, manifest_id=digest)


def build_v2_handoff_manifest_for_family(
    family_spec: V2ManifestFamilySpec,
    profile_config: V2ProfileGenerationConfig,
    generation_config: V2ManifestGenerationConfig,
) -> ResponseConstraintHandoffManifest:
    """Build one production v2 handoff manifest."""

    profile = build_v2_response_profile(profile_config)
    protocol = comparison_protocol_for_v2_family(family_spec.family_name)
    gate = ConstraintValidationGate(name="v2_manifest_generation_gate")
    manifest = build_candidate_handoff_manifest(
        profile,
        protocol,
        gate,
        max_constraints=generation_config.max_constraints,
        min_margin=generation_config.min_margin,
        train_fraction=generation_config.train_fraction,
        constraint_seed=generation_config.constraint_seed,
        bootstrap_count=generation_config.bootstrap_count,
        null_repetitions=generation_config.null_repetitions,
        source_label=family_spec.family_name,
    )
    manifest = replace(manifest, created_by_milestone="37")
    if family_spec.family_kind == "failed_control":
        manifest = replace(
            manifest,
            handoff_decision=HandoffDecision(
                eligible=False,
                failed_reasons=["v2_report_only_failed_control"],
                warning_reasons=["not_carry_forward_evaluated"],
            ),
        )
    return _with_manifest_id(manifest)


def build_v2_handoff_manifests(
    family_specs: list[V2ManifestFamilySpec],
    profile_configs: list[V2ProfileGenerationConfig],
    generation_config: V2ManifestGenerationConfig,
) -> list[ResponseConstraintHandoffManifest]:
    """Build all v2 handoff manifests."""

    config_by_family = {config.family_name: config for config in profile_configs}
    manifests: list[ResponseConstraintHandoffManifest] = []
    for spec in family_specs:
        manifests.append(
            build_v2_handoff_manifest_for_family(
                spec,
                config_by_family[spec.family_name],
                generation_config,
            )
        )
    return manifests


def write_v2_handoff_manifests(
    manifests: list[ResponseConstraintHandoffManifest],
    output_dir: Path,
) -> list[Path]:
    """Write production v2 handoff manifests under ``manifests_v2``."""

    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for manifest in manifests:
        safe_family = manifest.profile_label.replace("/", "_")
        path = output_dir / f"{safe_family}_{manifest.manifest_id[:12]}.json"
        paths.append(write_handoff_manifest(manifest, path))
    return paths
