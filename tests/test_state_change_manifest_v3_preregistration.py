from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_design import (
    default_v3_manifest_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v3_preregistration import (
    build_v3_preregistration_spec,
    v3_preregistration_digest,
    v3_preregistration_to_jsonable,
)


def test_v3_preregistration_spec_execution_disallowed_and_digest_stable() -> None:
    spec = build_v3_preregistration_spec(default_v3_manifest_family_designs())
    jsonable = v3_preregistration_to_jsonable(spec)

    assert spec.execution_allowed_in_current_milestone is False
    assert v3_preregistration_digest(jsonable) == v3_preregistration_digest(jsonable)

