from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v4_audit import (
    audit_v4_preregistration_only,
    check_no_v4_execution_outputs,
)
from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v4_preregistration import (
    build_v4_protocol_preregistration_spec,
)


def test_v4_audit_passes_report_only_spec() -> None:
    spec = build_v4_protocol_preregistration_spec(
        default_v4_protocol_family_designs()
    )
    row = audit_v4_preregistration_only(spec)

    assert row["execution_allowed_false"] == 1.0
    assert row["all_families_planned_only"] == 1.0


def test_check_no_v4_execution_outputs_detects_absent_outputs(tmp_path: Path) -> None:
    row = check_no_v4_execution_outputs(tmp_path)

    assert row["no_v4_production_manifests"] == 1.0
    assert row["no_v4_fit_outputs"] == 1.0
