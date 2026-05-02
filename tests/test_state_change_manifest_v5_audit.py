from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v5_audit import (
    audit_v5_preregistration_only,
    check_no_v5_execution_outputs,
)
from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v5_preregistration import (
    build_v5_protocol_preregistration_spec,
)


def test_v5_audit_detects_no_v5_production_outputs(tmp_path: Path) -> None:
    spec = build_v5_protocol_preregistration_spec(
        default_v5_protocol_family_designs()
    )
    audit = audit_v5_preregistration_only(spec)
    no_outputs = check_no_v5_execution_outputs(tmp_path)

    assert audit["execution_allowed_false"] == 1.0
    assert no_outputs["no_v5_production_manifests"] == 1.0
