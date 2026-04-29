from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    build_carry_forward_registry,
    forbidden_carry_forward_interpretations,
    registry_digest,
    registry_to_jsonable,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)


def _decision() -> FamilyRobustnessDecision:
    return FamilyRobustnessDecision(
        family_name="eligible_rank_gap",
        family_kind="structured",
        decision="carry_forward",
        passed=True,
        failed_reasons=[],
        warning_reasons=[],
        key_metrics={"mean_heldout_violation": 0.1},
    )


def test_forbidden_carry_forward_interpretations_nonempty() -> None:
    assert forbidden_carry_forward_interpretations()


def test_build_carry_forward_registry_serializes() -> None:
    registry = build_carry_forward_registry([_decision()])
    payload = registry_to_jsonable(registry)

    assert payload["records"]
    assert registry_digest(payload) == registry.registry_id
