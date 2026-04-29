from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import (
    forbidden_manifest_fields,
    load_manifest_dataset,
    manifest_integrity_report,
)


def test_load_manifest_dataset_validates_shape_and_splits(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)

    assert dataset.constraints.ndim == 2
    assert dataset.constraints.shape[1] == 4
    assert set(dataset.train_constraint_indices).isdisjoint(
        set(dataset.heldout_constraint_indices)
    )
    assert dataset.train_constraints.shape[0] == dataset.train_constraint_indices.size


def test_manifest_integrity_report_detects_no_embedding_fields(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)
    report = manifest_integrity_report(dataset)

    assert report["has_embedding_fields"] == 0.0
    assert report["has_metric_fields"] == 0.0
    assert report["passed_integrity"] == 1.0
    assert "embedding" in forbidden_manifest_fields()


def test_load_manifest_dataset_rejects_embedding_field(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["embedding"] = [[0.0]]
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="forbidden"):
        load_manifest_dataset(path)
