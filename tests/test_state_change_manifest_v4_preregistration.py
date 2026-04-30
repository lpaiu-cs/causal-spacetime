from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v4_preregistration import (
    build_v4_protocol_preregistration_spec,
    v4_protocol_preregistration_digest,
    v4_protocol_preregistration_to_jsonable,
)


def test_v4_preregistration_execution_allowed_false() -> None:
    spec = build_v4_protocol_preregistration_spec(
        default_v4_protocol_family_designs()
    )

    assert not spec.execution_allowed_in_current_milestone
    assert spec.created_by_milestone == "Milestone 43"


def test_v4_preregistration_digest_is_stable() -> None:
    spec = build_v4_protocol_preregistration_spec(
        default_v4_protocol_family_designs()
    )
    jsonable = v4_protocol_preregistration_to_jsonable(spec)

    assert v4_protocol_preregistration_digest(jsonable) == (
        v4_protocol_preregistration_digest(jsonable)
    )
