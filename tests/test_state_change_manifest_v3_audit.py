from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_audit import (
    audit_v3_preregistration_only,
    check_no_v3_execution_outputs,
)
from causal_spacetime_lab.state_change_manifest_v3_design import (
    default_v3_manifest_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v3_preregistration import (
    build_v3_preregistration_spec,
)


def test_v3_audit_passes_report_only_spec(tmp_path: Path) -> None:
    spec = build_v3_preregistration_spec(default_v3_manifest_family_designs())
    audit = audit_v3_preregistration_only(spec)
    execution = check_no_v3_execution_outputs(tmp_path)

    assert all(float(value) == 1.0 for value in audit.values())
    assert all(float(value) == 1.0 for value in execution.values())

