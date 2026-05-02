"""No-execution audits for planned-only v5 protocol preregistration."""

from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v5_preregistration import (
    V5ProtocolPreregistrationSpec,
)


def audit_v5_preregistration_only(
    spec: V5ProtocolPreregistrationSpec,
) -> dict[str, float | str]:
    """Audit that v5 remains planned-only in Milestone 46."""

    structured = [
        family
        for family in spec.planned_families
        if family.family_kind == "structured"
    ]
    return {
        "execution_allowed_false": float(
            not spec.execution_allowed_in_current_milestone
        ),
        "all_families_planned_only": float(
            all(
                family.execution_status == "planned_only"
                for family in spec.planned_families
            )
        ),
        "structured_manifest_count_at_least_3": float(
            all(family.planned_manifest_count >= 3 for family in structured)
        ),
        "forbids_threshold_retuning": float(
            "threshold retuning" in spec.forbidden_actions
        ),
        "forbidden_interpretations_nonempty": float(
            bool(spec.forbidden_interpretations)
        ),
        "remediation_iteration_risk_audit_required": float(
            spec.remediation_iteration_risk_audit_required
        ),
    }


def check_no_v5_execution_outputs(output_dir: Path) -> dict[str, float | str]:
    """Check that M46 did not create v5 production outputs."""

    manifest_dir = output_dir / "manifests_v5"
    data_dir = output_dir / "data"
    registry_dir = output_dir / "carry_forward_v5"
    registry_protocol_dir = output_dir / "carry_forward_v5_protocol"
    return {
        "no_v5_production_manifests": float(
            not manifest_dir.exists() or not list(manifest_dir.glob("*.json"))
        ),
        "no_v5_carry_forward_registry": float(
            not registry_dir.exists() and not registry_protocol_dir.exists()
        ),
        "no_v5_stress_test_outputs": float(
            not any(data_dir.glob("v5_*stress*.csv"))
        ),
        "no_v5_fit_outputs": float(
            not any(data_dir.glob("v5_*fit*.csv"))
            and not any(data_dir.glob("v5_*representation*.csv"))
        ),
    }
