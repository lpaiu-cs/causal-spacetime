"""Generate preregistered protocol-invariant v4 handoff manifests."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    V4ProtocolExecutionFamilySpec,
)
from causal_spacetime_lab.state_change_manifest_v4_profiles import (
    V4ProtocolProfileGenerationConfig,
    build_v4_protocol_response_profile,
)
from causal_spacetime_lab.state_change_manifest_v4_protocol_mapping import (
    comparison_protocol_for_v4_family,
)
from causal_spacetime_lab.state_change_manifest_v4_provenance import (
    handoff_provenance_for_v4_family,
    top_down_template_for_v4_family,
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
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_jsonable,
)
from causal_spacetime_lab.state_change_top_down_handoff_template import (
    top_down_template_digest,
    top_down_template_jsonable,
)


@dataclass(frozen=True)
class V4ProtocolManifestGenerationConfig:
    """Generation config for preregistered v4 handoff manifests."""

    max_constraints: int = 5000
    train_fraction: float = 0.7
    bootstrap_count: int = 20
    null_repetitions: int = 5
    constraint_seed: int = 0


def build_v4_protocol_handoff_manifest_for_config(
    family_spec: V4ProtocolExecutionFamilySpec,
    profile_config: V4ProtocolProfileGenerationConfig,
    generation_config: V4ProtocolManifestGenerationConfig,
) -> ResponseConstraintHandoffManifest:
    """Build one production v4 handoff manifest."""

    profile_with_metadata = build_v4_protocol_response_profile(profile_config)
    comparison_protocol = comparison_protocol_for_v4_family(
        family_spec,
        profile_config.measurement_protocol,
    )
    template = top_down_template_for_v4_family(
        family_spec,
        profile_config.measurement_protocol,
    )
    design_digest = _design_digest_from_path(
        Path("outputs/remediation/v4_protocol_preregistration_spec_m43.json")
    )
    provenance = handoff_provenance_for_v4_family(
        family_spec,
        template,
        profile_config,
        design_digest=design_digest,
    )
    gate = ConstraintValidationGate(name="v4_protocol_manifest_generation_gate")
    source_label = f"{family_spec.family_name}_m{profile_config.manifest_index}"
    if family_spec.family_kind in {"failed_control", "report_only"}:
        manifest = build_candidate_handoff_manifest(
            profile_with_metadata.profile,
            comparison_protocol,
            gate,
            max_constraints=generation_config.max_constraints,
            min_margin=family_spec.margin_value,
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
            comparison_protocol,
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
    decision = manifest.handoff_decision
    admissible = family_spec.family_kind == "structured"
    if not admissible:
        decision = HandoffDecision(
            eligible=False,
            failed_reasons=["report_only_or_protocol_inadmissible"],
            warning_reasons=["not_carry_forward_evaluated"],
        )
    manifest = replace(
        manifest,
        created_by_milestone="44",
        handoff_decision=decision,
        profile_metadata=profile_json,
        measurement_protocol=profile_json["measurement_protocol"],
        measurement_protocol_id=str(profile_json["measurement_protocol_id"]),
        measurement_protocol_hash=str(profile_json["measurement_protocol_hash"]),
        reference_set_id=str(profile_json["reference_set_id"]),
        reference_chain_ids=list(profile_with_metadata.metadata.reference_chain_ids),
        profile_invariance_status=profile_with_metadata.metadata.profile_invariance_status,
        admissible_for_pairwise_dissimilarity=admissible,
        handoff_provenance=provenance.__dict__.copy(),
        handoff_provenance_type=provenance.provenance_type,
        handoff_design_digest=provenance.design_digest,
        top_down_template=template_json,
        top_down_template_id=template.template_id,
        top_down_template_hash=top_down_template_digest(template),
    )
    digest = manifest_digest(manifest_to_jsonable(manifest))
    return replace(manifest, manifest_id=digest)


def build_v4_protocol_handoff_manifests(
    family_specs: list[V4ProtocolExecutionFamilySpec],
    profile_configs: list[V4ProtocolProfileGenerationConfig],
    generation_config: V4ProtocolManifestGenerationConfig,
) -> list[ResponseConstraintHandoffManifest]:
    """Build all preregistered v4 manifests."""

    spec_by_family = {spec.family_name: spec for spec in family_specs}
    return [
        build_v4_protocol_handoff_manifest_for_config(
            spec_by_family[config.family_name],
            config,
            generation_config,
        )
        for config in profile_configs
    ]


def write_v4_protocol_handoff_manifests(
    manifests: list[ResponseConstraintHandoffManifest],
    output_dir: Path,
    *,
    overwrite: bool = True,
) -> list[Path]:
    """Write production v4 handoff manifests under ``outputs/manifests_v4``."""

    if output_dir.name != "manifests_v4":
        raise ValueError("v4 manifests must be written under outputs/manifests_v4")
    output_dir.mkdir(parents=True, exist_ok=True)
    incompatible = [
        path
        for path in output_dir.glob("*.json")
        if not path.name.startswith("m44_") and not _manifest_has_m44_metadata(path)
    ]
    if incompatible:
        raise FileExistsError(
            "outputs/manifests_v4 contains incompatible manifest files: "
            + ", ".join(path.name for path in incompatible[:5])
        )
    if overwrite:
        for path in output_dir.glob("m44_*.json"):
            path.unlink()
    paths: list[Path] = []
    for manifest in manifests:
        safe_family = manifest.profile_label.replace("/", "_")
        path = output_dir / f"m44_{safe_family}_{manifest.manifest_id[:12]}.json"
        paths.append(write_handoff_manifest(manifest, path))
    return paths


def _design_digest_from_path(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest_has_m44_metadata(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    return (
        "measurement_protocol_id" in text
        and "profile_metadata" in text
        and "handoff_provenance" in text
        and '"created_by_milestone": "44"' in text
    )
