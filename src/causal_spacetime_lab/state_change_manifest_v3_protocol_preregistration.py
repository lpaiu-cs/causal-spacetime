"""Patched v3 preregistration after protocol-invariance audit."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    V3ProtocolInvariantFamilyPatch,
    default_v3_protocol_invariant_family_patches,
    v3_protocol_patch_table,
)


@dataclass(frozen=True)
class V3ProtocolPatchedPreregistrationSpec:
    """Planned-only v3 protocol-patched preregistration spec."""

    spec_id: str
    created_by_milestone: str
    based_on_spec: str
    patch_reason: str
    patched_families: list[V3ProtocolInvariantFamilyPatch]
    execution_allowed_in_current_milestone: bool
    forbidden_actions: list[str]
    forbidden_interpretations: list[str]


def forbidden_protocol_patch_interpretations() -> list[str]:
    """Return forbidden interpretations for the protocol patch."""

    return [
        "protocol patch proves geometry",
        "protocol-invariant profile is distance",
        "patched v3 family will pass",
        "measurement-protocol metadata recovers metric",
    ]


def _forbidden_actions() -> list[str]:
    return [
        "execute_v3_generation_in_m40",
        "run_stress_tests",
        "retune_thresholds",
        "change_v2_decisions",
        "change_m34_criteria",
    ]


def build_v3_protocol_patched_preregistration(
    base_spec_path: Path,
    patches: list[V3ProtocolInvariantFamilyPatch],
) -> V3ProtocolPatchedPreregistrationSpec:
    """Build a planned-only protocol-patched v3 preregistration spec."""

    if not patches:
        patches = default_v3_protocol_invariant_family_patches()
    based_on_spec = str(base_spec_path)
    if base_spec_path.exists():
        try:
            payload = json.loads(base_spec_path.read_text(encoding="utf-8"))
            based_on_spec = str(payload.get("spec_id") or base_spec_path)
        except json.JSONDecodeError:
            based_on_spec = str(base_spec_path)
    spec = V3ProtocolPatchedPreregistrationSpec(
        spec_id="pending",
        created_by_milestone="Milestone 40",
        based_on_spec=based_on_spec,
        patch_reason=(
            "pre-execution design correction: one response profile may vary "
            "reference chains only inside a fixed measurement protocol"
        ),
        patched_families=patches,
        execution_allowed_in_current_milestone=False,
        forbidden_actions=_forbidden_actions(),
        forbidden_interpretations=forbidden_protocol_patch_interpretations(),
    )
    jsonable = v3_protocol_patched_preregistration_to_jsonable(spec)
    digest = _digest({**jsonable, "spec_id": "pending"})
    return V3ProtocolPatchedPreregistrationSpec(
        spec_id=digest,
        created_by_milestone=spec.created_by_milestone,
        based_on_spec=spec.based_on_spec,
        patch_reason=spec.patch_reason,
        patched_families=spec.patched_families,
        execution_allowed_in_current_milestone=(
            spec.execution_allowed_in_current_milestone
        ),
        forbidden_actions=spec.forbidden_actions,
        forbidden_interpretations=spec.forbidden_interpretations,
    )


def v3_protocol_patched_preregistration_to_jsonable(
    spec: V3ProtocolPatchedPreregistrationSpec,
) -> dict[str, object]:
    """Convert a protocol-patched preregistration spec to JSON data."""

    return {
        "spec_id": spec.spec_id,
        "created_by_milestone": spec.created_by_milestone,
        "based_on_spec": spec.based_on_spec,
        "patch_reason": spec.patch_reason,
        "patched_families": v3_protocol_patch_table(spec.patched_families),
        "execution_allowed_in_current_milestone": bool(
            spec.execution_allowed_in_current_milestone
        ),
        "forbidden_actions": list(spec.forbidden_actions),
        "forbidden_interpretations": list(spec.forbidden_interpretations),
    }


def write_v3_protocol_patched_preregistration(
    spec: V3ProtocolPatchedPreregistrationSpec,
    output_path: Path,
) -> Path:
    """Write planned-only protocol-patched v3 preregistration JSON."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(
            v3_protocol_patched_preregistration_to_jsonable(spec),
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return output_path


def _digest(jsonable: dict[str, object]) -> str:
    payload = json.dumps(jsonable, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
