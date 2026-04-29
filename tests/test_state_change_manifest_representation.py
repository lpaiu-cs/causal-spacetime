from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.abspath("experiments"))

from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_dataset
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_ordinal_representation,
    representation_fit_to_row,
)


def test_fit_skips_ineligible_by_default(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["handoff_decision"]["eligible"] = False
    payload["handoff_decision"]["failed_reasons"] = ["test_failure"]
    path.write_text(json.dumps(payload), encoding="utf-8")
    dataset = load_manifest_dataset(path)

    fit = fit_manifest_ordinal_representation(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=20, restarts=1),
    )

    assert not fit.fitted
    assert fit.reason_not_fit == "manifest_ineligible"


def test_fit_manifest_ordinal_representation_uses_train_constraints(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)

    fit = fit_manifest_ordinal_representation(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=30, restarts=1),
    )

    assert fit.fitted
    assert fit.train_constraint_count == dataset.train_constraints.shape[0]
    assert fit.heldout_constraint_count == dataset.heldout_constraints.shape[0]
    assert fit.embedding is not None


def test_representation_fit_to_row_excludes_embedding_by_default(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)
    fit = fit_manifest_ordinal_representation(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=20, restarts=1),
    )

    row = representation_fit_to_row(fit)

    assert "embedding_values" not in row
    assert row["representation_kind"] == "latent_ordinal_representation"
