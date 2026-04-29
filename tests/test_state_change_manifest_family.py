from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.abspath("experiments"))

from manifest_family_experiment_helpers import ensure_failed_control_manifest
from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_dataset
from causal_spacetime_lab.state_change_manifest_family import (
    assign_manifest_to_family,
    default_manifest_family_specs,
)


def test_default_manifest_family_specs_returns_expected_families() -> None:
    specs = default_manifest_family_specs()

    names = {spec.family_name for spec in specs}
    assert "eligible_rank_gap" in names
    assert "ineligible_reported" in names


def test_assign_manifest_to_family_deterministic_eligible_gap(tmp_path) -> None:  # type: ignore[no-untyped-def]
    dataset = load_manifest_dataset(build_exact_manifest(tmp_path, "manifest.json"))

    assignment = assign_manifest_to_family(dataset, default_manifest_family_specs())

    assert assignment.family_name == "eligible_rank_gap"
    assert assignment.family_kind == "structured"


def test_assign_manifest_to_family_deterministic_ineligible(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["manifest_id"] = "ordinary_ineligible"
    payload["handoff_decision"]["eligible"] = False
    payload["handoff_decision"]["failed_reasons"] = ["low_heldout_agreement"]
    path.write_text(json.dumps(payload), encoding="utf-8")
    dataset = load_manifest_dataset(path)

    assignment = assign_manifest_to_family(dataset, default_manifest_family_specs())

    assert assignment.family_name == "ineligible_reported"


def test_assign_manifest_to_family_failed_control(tmp_path) -> None:  # type: ignore[no-untyped-def]
    dataset = ensure_failed_control_manifest(tmp_path)

    assignment = assign_manifest_to_family(dataset, default_manifest_family_specs())

    assert assignment.family_name == "failed_synthetic_controls"
