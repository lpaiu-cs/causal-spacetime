from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v5_preregistration import (
    build_v5_protocol_preregistration_spec,
    v5_protocol_preregistration_to_jsonable,
)


def test_v5_preregistration_execution_allowed_false() -> None:
    spec = build_v5_protocol_preregistration_spec(
        default_v5_protocol_family_designs()
    )
    jsonable = v5_protocol_preregistration_to_jsonable(spec)

    assert spec.execution_allowed_in_current_milestone is False
    assert jsonable["remediation_iteration_risk_audit_required"] is True
