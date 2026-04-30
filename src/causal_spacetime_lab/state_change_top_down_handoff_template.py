"""Top-down handoff templates for preregistered response-comparison handoffs."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    V3ProtocolInvariantFamilyPatch,
)

CONSTRAINT_SELECTION_POLICIES = {
    "profile_instantiated_rank_gap",
    "profile_instantiated_combined_gap",
    "declared_constraint_schema_only",
    "report_only_no_constraints",
}


@dataclass(frozen=True)
class TopDownHandoffTemplate:
    """Preregistered source of handoff family and constraint schema."""

    template_id: str
    template_label: str
    family_name: str
    measurement_protocol_id: str
    measurement_protocol_hash: str
    comparison_method: str
    missing_policy: str
    min_common_protocols: int
    margin_policy: str
    margin_value: float
    constraint_selection_policy: str
    target_inclusion_policy: str
    reference_set_policy: str
    forbidden_interpretations: list[str]

    def __post_init__(self) -> None:
        if self.constraint_selection_policy not in CONSTRAINT_SELECTION_POLICIES:
            allowed = ", ".join(sorted(CONSTRAINT_SELECTION_POLICIES))
            raise ValueError(
                f"constraint_selection_policy must be one of: {allowed}"
            )


def top_down_template_jsonable(
    template: TopDownHandoffTemplate,
) -> dict[str, object]:
    """Return JSON-compatible top-down template metadata."""

    return asdict(template)


def top_down_template_digest(template: TopDownHandoffTemplate) -> str:
    """Return stable digest for a top-down handoff template."""

    payload = dict(top_down_template_jsonable(template))
    payload.pop("template_id", None)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def template_from_v3_protocol_family_spec(
    family_spec: V3ProtocolInvariantFamilyPatch,
) -> TopDownHandoffTemplate:
    """Build a top-down template from a patched v3 protocol family."""

    if family_spec.family_kind in {"failed_control", "report_only"}:
        policy = "report_only_no_constraints"
    elif family_spec.patched_family_name.startswith("combined_"):
        policy = "profile_instantiated_combined_gap"
    else:
        policy = "profile_instantiated_rank_gap"
    comparison_method = (
        "combined_gap_and_mismatch"
        if family_spec.patched_family_name.startswith("combined_")
        else "rank_gap_mean"
    )
    template = TopDownHandoffTemplate(
        template_id="pending",
        template_label=f"{family_spec.patched_family_name}_template",
        family_name=family_spec.patched_family_name,
        measurement_protocol_id=family_spec.measurement_protocol_id,
        measurement_protocol_hash=family_spec.measurement_protocol_hash,
        comparison_method=comparison_method,
        missing_policy=family_spec.missing_policy,
        min_common_protocols=1,
        margin_policy=family_spec.margin_policy,
        margin_value=0.05,
        constraint_selection_policy=policy,
        target_inclusion_policy="predeclared_profile_targets",
        reference_set_policy="reference-chain variation inside fixed protocol",
        forbidden_interpretations=[
            "metric reconstruction",
            "geometry recovery",
            "top-down handoff is evidence of geometry",
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
