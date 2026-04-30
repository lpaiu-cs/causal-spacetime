from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preregistration import (
    build_v3_protocol_patched_preregistration,
    v3_protocol_patched_preregistration_to_jsonable,
)


def test_patched_preregistration_execution_disallowed() -> None:
    spec = build_v3_protocol_patched_preregistration(
        Path("missing.json"),
        default_v3_protocol_invariant_family_patches(),
    )
    jsonable = v3_protocol_patched_preregistration_to_jsonable(spec)

    assert spec.execution_allowed_in_current_milestone is False
    assert jsonable["execution_allowed_in_current_milestone"] is False
    assert spec.forbidden_interpretations
