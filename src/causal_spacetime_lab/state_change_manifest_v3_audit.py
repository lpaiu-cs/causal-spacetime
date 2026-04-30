"""Audits for planned-only v3 preregistration."""

from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_preregistration import (
    V3PreregistrationSpec,
)


def audit_v3_preregistration_only(
    spec: V3PreregistrationSpec,
) -> dict[str, float | str]:
    """Audit that v3 is preregistered but not executed in Milestone 39."""

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
    }


def check_no_v3_execution_outputs(
    output_dir: Path,
) -> dict[str, float | str]:
    """Check that M39 did not create v3 production outputs."""

    manifest_dir = output_dir / "manifests_v3"
    data_dir = output_dir / "data"
    registry_path = output_dir / "carry_forward_v3" / "carry_forward_registry_v3.json"
    return {
        "no_v3_production_manifests": float(
            not manifest_dir.exists() or not list(manifest_dir.glob("*.json"))
        ),
        "no_v3_carry_forward_registry": float(
            not registry_path.exists()
        ),
        "no_v3_stress_test_outputs": float(
            not any(data_dir.glob("v3_*stress*.csv"))
        ),
        "no_v3_fit_outputs": float(
            not any(data_dir.glob("v3_*fit*.csv"))
            and not any(data_dir.glob("v3_*representation*.csv"))
        ),
    }
