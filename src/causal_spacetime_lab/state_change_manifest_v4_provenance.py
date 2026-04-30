"""Top-down template and provenance helpers for v4 handoff manifests."""

from __future__ import annotations

from dataclasses import replace

from causal_spacetime_lab.state_change_handoff_provenance import (
    HandoffProvenanceMetadata,
    default_hybrid_handoff_provenance,
    default_report_only_handoff_provenance,
)
from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    V4ProtocolExecutionFamilySpec,
)
from causal_spacetime_lab.state_change_manifest_v4_profiles import (
    V4ProtocolProfileGenerationConfig,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    measurement_protocol_hash,
    measurement_protocol_id,
)
from causal_spacetime_lab.state_change_top_down_handoff_template import (
    TopDownHandoffTemplate,
    top_down_template_digest,
)


def top_down_template_for_v4_family(
    spec: V4ProtocolExecutionFamilySpec,
    measurement_protocol: MeasurementProtocolSpec,
) -> TopDownHandoffTemplate:
    """Build a preregistered v4 handoff template."""

    if spec.family_kind in {"failed_control", "report_only"}:
        policy = "report_only_no_constraints"
    elif spec.comparison_method == "combined_gap_and_mismatch":
        policy = "profile_instantiated_combined_gap"
    else:
        policy = "profile_instantiated_rank_gap"
    template = TopDownHandoffTemplate(
        template_id="pending",
        template_label=f"{spec.family_name}_template",
        family_name=spec.family_name,
        measurement_protocol_id=measurement_protocol_id(measurement_protocol),
        measurement_protocol_hash=measurement_protocol_hash(measurement_protocol),
        comparison_method=spec.comparison_method,
        missing_policy=measurement_protocol.missing_policy,
        min_common_protocols=1,
        margin_policy=spec.margin_policy,
        margin_value=float(spec.margin_value),
        constraint_selection_policy=policy,
        target_inclusion_policy=spec.target_inclusion_policy,
        reference_set_policy=spec.profile_resolution_policy,
        forbidden_interpretations=[
            "metric reconstruction",
            "geometry recovery",
            "diagnostic-complete output is not carry-forward success",
        ],
    )
    digest = top_down_template_digest(template)
    return TopDownHandoffTemplate(
        template_id=digest,
        template_label=template.template_label,
        family_name=template.family_name,
        measurement_protocol_id=template.measurement_protocol_id,
        measurement_protocol_hash=template.measurement_protocol_hash,
        comparison_method=template.comparison_method,
        missing_policy=template.missing_policy,
        min_common_protocols=template.min_common_protocols,
        margin_policy=template.margin_policy,
        margin_value=template.margin_value,
        constraint_selection_policy=template.constraint_selection_policy,
        target_inclusion_policy=template.target_inclusion_policy,
        reference_set_policy=template.reference_set_policy,
        forbidden_interpretations=template.forbidden_interpretations,
    )


def handoff_provenance_for_v4_family(
    spec: V4ProtocolExecutionFamilySpec,
    template: TopDownHandoffTemplate,
    profile_config: V4ProtocolProfileGenerationConfig,
    *,
    design_source_path: str = (
        "outputs/remediation/v4_protocol_preregistration_spec_m43.json"
    ),
    design_digest: str = "",
) -> HandoffProvenanceMetadata:
    """Return provenance for one v4 handoff manifest."""

    template_hash = top_down_template_digest(template)
    digest = design_digest or template_hash
    if spec.family_kind in {"failed_control", "report_only"}:
        metadata = default_report_only_handoff_provenance(
            design_source_label=spec.family_name,
            design_source_path=design_source_path,
            design_digest=digest,
            reason=(
                "failed-control accounting"
                if spec.family_kind == "failed_control"
                else "report-only control"
            ),
        )
        return replace(metadata, created_by_milestone="Milestone 44")
    metadata = default_hybrid_handoff_provenance(
        design_source_label=spec.family_name,
        design_source_path=design_source_path,
        design_digest=digest,
        constraint_selection_policy=template.constraint_selection_policy,
        template_id=template.template_id,
        template_hash=template_hash,
    )
    return replace(metadata, created_by_milestone="Milestone 44")
