"""Generate protocol-invariant, provenance-aware patched v3 manifests."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from pathlib import Path

from causal_spacetime_lab.state_change_handoff_provenance import (
    HandoffProvenanceMetadata,
    default_hybrid_handoff_provenance,
    default_report_only_handoff_provenance,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    V3ProtocolExecutionFamilySpec,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_profiles import (
    V3ProtocolProfileGenerationConfig,
    build_v3_protocol_response_profile,
)
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
    build_candidate_handoff_manifest_checked,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_jsonable,
)
from causal_spacetime_lab.state_change_top_down_handoff_template import (
    TopDownHandoffTemplate,
    template_from_v3_protocol_family_spec,
    top_down_template_digest,
    top_down_template_jsonable,
)


@dataclass(frozen=True)
class V3ProtocolManifestGenerationConfig:
    """Generation config for patched v3 handoff manifests."""

    max_constraints: int = 5000
    train_fraction: float = 0.7
    bootstrap_count: int = 20
    null_repetitions: int = 5
    constraint_seed: int = 0


def comparison_protocol_for_v3_protocol_family(
    family_spec: V3ProtocolExecutionFamilySpec,
) -> PairwiseResponseComparisonProtocol:
    """Return pairwise comparison protocol for a patched v3 family."""

    method = family_spec.comparison_method
    penalty = 0.5 if method == "combined_gap_and_mismatch" else 1.0
    return PairwiseResponseComparisonProtocol(
        name=f"{family_spec.family_name}_{method}",
        method=method,
        missing_policy=family_spec.measurement_protocol.missing_policy,
        min_common_protocols=1,
        reachability_mismatch_penalty=penalty,
    )


def top_down_template_for_v3_protocol_family(
    family_spec: V3ProtocolExecutionFamilySpec,
) -> TopDownHandoffTemplate:
    """Build a top-down handoff template for a patched v3 family."""

    patch_like = _patch_like_from_execution_spec(family_spec)
    template = template_from_v3_protocol_family_spec(patch_like)
    if template.comparison_method == family_spec.comparison_method:
        return template
    return replace(template, comparison_method=family_spec.comparison_method)


def handoff_provenance_for_v3_protocol_family(
    family_spec: V3ProtocolExecutionFamilySpec,
    template: TopDownHandoffTemplate,
    profile_config: V3ProtocolProfileGenerationConfig,
    *,
    design_source_path: str = (
        "outputs/remediation/v3_protocol_patched_preregistration_m40.json"
    ),
    design_digest: str = "",
) -> HandoffProvenanceMetadata:
    """Return provenance metadata for one patched v3 handoff manifest."""

    template_hash = top_down_template_digest(template)
    if family_spec.family_kind in {"failed_control", "report_only"}:
        return default_report_only_handoff_provenance(
            design_source_label=family_spec.family_name,
            design_source_path=design_source_path,
            design_digest=design_digest or template_hash,
            reason="report-only or failed-control patched v3 family",
        )
    return default_hybrid_handoff_provenance(
        design_source_label=family_spec.family_name,
        design_source_path=design_source_path,
        design_digest=design_digest or template_hash,
        constraint_selection_policy=template.constraint_selection_policy,
        template_id=template.template_id,
        template_hash=template_hash,
    )


def build_v3_protocol_handoff_manifest_for_config(
    family_spec: V3ProtocolExecutionFamilySpec,
    profile_config: V3ProtocolProfileGenerationConfig,
    generation_config: V3ProtocolManifestGenerationConfig,
) -> ResponseConstraintHandoffManifest:
    """Build one protocol-invariant, provenance-aware v3 handoff manifest."""

    profile_with_metadata = build_v3_protocol_response_profile(profile_config)
    protocol = comparison_protocol_for_v3_protocol_family(family_spec)
    if protocol.missing_policy != family_spec.measurement_protocol.missing_policy:
        raise ValueError("comparison missing policy must match measurement protocol")
    template = top_down_template_for_v3_protocol_family(family_spec)
    design_digest = _design_digest_from_path(
        Path("outputs/remediation/v3_protocol_patched_preregistration_m40.json")
    )
    provenance = handoff_provenance_for_v3_protocol_family(
        family_spec,
        template,
        profile_config,
        design_digest=design_digest,
    )
    gate = ConstraintValidationGate(name="v3_protocol_manifest_generation_gate")
    source_label = f"{family_spec.family_name}_m{profile_config.manifest_index}"
    if family_spec.family_kind in {"failed_control", "report_only"}:
        manifest = build_candidate_handoff_manifest(
            profile_with_metadata.profile,
            protocol,
            gate,
            max_constraints=generation_config.max_constraints,
            min_margin=family_spec.min_margin,
            train_fraction=generation_config.train_fraction,
            constraint_seed=(
                generation_config.constraint_seed + profile_config.manifest_index
            ),
            bootstrap_count=generation_config.bootstrap_count,
            null_repetitions=generation_config.null_repetitions,
            source_label=source_label,
        )
    else:
        manifest = build_candidate_handoff_manifest_checked(
            profile_with_metadata,
            protocol,
            gate,
            provenance,
            max_constraints=generation_config.max_constraints,
            train_fraction=generation_config.train_fraction,
            constraint_seed=(
                generation_config.constraint_seed + profile_config.manifest_index
            ),
            bootstrap_count=generation_config.bootstrap_count,
            null_repetitions=generation_config.null_repetitions,
            source_label=source_label,
        )
    profile_json = profile_metadata_jsonable(profile_with_metadata.metadata)
    template_json = top_down_template_jsonable(template)
    provenance_json = provenance.__dict__.copy()
    decision = manifest.handoff_decision
    if family_spec.family_kind in {"failed_control", "report_only"}:
        decision = HandoffDecision(
            eligible=False,
            failed_reasons=["report_only_or_protocol_inadmissible"],
            warning_reasons=["not_carry_forward_evaluated"],
        )
    manifest = replace(
        manifest,
        created_by_milestone="41",
        handoff_decision=decision,
        profile_metadata=profile_json,
        measurement_protocol=profile_json["measurement_protocol"],
        measurement_protocol_id=str(profile_json["measurement_protocol_id"]),
        measurement_protocol_hash=str(profile_json["measurement_protocol_hash"]),
        reference_set_id=str(profile_json["reference_set_id"]),
        reference_chain_ids=list(profile_with_metadata.metadata.reference_chain_ids),
        profile_invariance_status=profile_with_metadata.metadata.profile_invariance_status,
        admissible_for_pairwise_dissimilarity=(
            family_spec.admissible_for_pairwise_dissimilarity
        ),
        handoff_provenance=provenance_json,
        handoff_provenance_type=provenance.provenance_type,
        handoff_design_digest=provenance.design_digest,
        top_down_template=template_json,
        top_down_template_id=template.template_id,
        top_down_template_hash=top_down_template_digest(template),
    )
    digest = manifest_digest(manifest_to_jsonable(manifest))
    return replace(manifest, manifest_id=digest)


def build_v3_protocol_handoff_manifests(
    family_specs: list[V3ProtocolExecutionFamilySpec],
    profile_configs: list[V3ProtocolProfileGenerationConfig],
    generation_config: V3ProtocolManifestGenerationConfig,
) -> list[ResponseConstraintHandoffManifest]:
    """Build all protocol-invariant patched v3 manifests."""

    spec_by_family = {spec.family_name: spec for spec in family_specs}
    manifests: list[ResponseConstraintHandoffManifest] = []
    for config in profile_configs:
        manifests.append(
            build_v3_protocol_handoff_manifest_for_config(
                spec_by_family[config.family_name],
                config,
                generation_config,
            )
        )
    return manifests


def write_v3_protocol_handoff_manifests(
    manifests: list[ResponseConstraintHandoffManifest],
    output_dir: Path,
    *,
    overwrite: bool = True,
) -> list[Path]:
    """Write patched v3 manifests under ``outputs/manifests_v3``."""

    output_dir.mkdir(parents=True, exist_ok=True)
    existing = list(output_dir.glob("*.json"))
    incompatible = [
        path
        for path in existing
        if not path.name.startswith("m41_") and not _manifest_has_m41_metadata(path)
    ]
    if incompatible:
        raise FileExistsError(
            "outputs/manifests_v3 contains incompatible manifest files: "
            + ", ".join(path.name for path in incompatible[:5])
        )
    if overwrite:
        for path in output_dir.glob("m41_*.json"):
            path.unlink()
    paths: list[Path] = []
    for manifest in manifests:
        safe_family = manifest.profile_label.replace("/", "_")
        path = output_dir / f"m41_{safe_family}_{manifest.manifest_id[:12]}.json"
        paths.append(write_handoff_manifest(manifest, path))
    return paths


def _design_digest_from_path(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest_has_m41_metadata(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    return (
        "measurement_protocol_id" in text
        and "profile_metadata" in text
        and "handoff_provenance" in text
    )


def _patch_like_from_execution_spec(family_spec: V3ProtocolExecutionFamilySpec):
    from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
        V3ProtocolInvariantFamilyPatch,
    )

    protocol = family_spec.measurement_protocol
    return V3ProtocolInvariantFamilyPatch(
        original_family_name=family_spec.family_name,
        patched_family_name=family_spec.family_name,
        family_kind=family_spec.family_kind,
        measurement_protocol_id=family_spec.measurement_protocol_id,
        measurement_protocol_hash=family_spec.measurement_protocol_hash,
        echo_rule=protocol.echo_rule,
        emission_rule=protocol.emission_rule,
        gate_rule=protocol.gate_rule,
        reference_subsampling_rule=protocol.reference_subsampling_rule,
        spectrum_type=protocol.spectrum_type,
        normalization_rule=protocol.normalization_rule,
        missing_policy=protocol.missing_policy,
        tie_policy=protocol.tie_policy,
        margin_policy=protocol.margin_policy,
        planned_manifest_count=family_spec.planned_manifest_count,
        profile_family_semantics=family_spec.profile_family_semantics,
        admissible_for_pairwise_dissimilarity=(
            family_spec.admissible_for_pairwise_dissimilarity
        ),
        execution_status="planned_only",
    )
