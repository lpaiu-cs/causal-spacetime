from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath("experiments"))

from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_dataset
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
)
from causal_spacetime_lab.state_change_manifest_representation_nulls import (
    evaluate_manifest_representation_nulls,
    permute_manifest_target_labels,
    random_manifest_constraints,
    shuffle_constraint_sides_for_manifest,
)


def test_shuffle_constraint_sides_for_manifest_preserves_shape() -> None:
    constraints = np.asarray([[0, 1, 2, 3], [1, 2, 3, 4]], dtype=int)

    shuffled = shuffle_constraint_sides_for_manifest(constraints, seed=0)

    assert shuffled.shape == constraints.shape


def test_random_manifest_constraints_generates_valid_quadruplets() -> None:
    constraints = random_manifest_constraints(5, 20, seed=1)

    assert constraints.shape == (20, 4)
    assert np.all(constraints >= 0)
    assert np.all(constraints < 5)
    assert np.all(constraints[:, 0] != constraints[:, 1])
    assert np.all(constraints[:, 2] != constraints[:, 3])


def test_permute_manifest_target_labels_preserves_valid_ids() -> None:
    constraints = np.asarray([[0, 1, 2, 3], [1, 2, 3, 4]], dtype=int)

    permuted = permute_manifest_target_labels(constraints, 5, seed=2)

    assert permuted.shape == constraints.shape
    assert np.all(permuted >= 0)
    assert np.all(permuted < 5)


def test_evaluate_manifest_nulls_returns_expected_summaries(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)

    summaries = evaluate_manifest_representation_nulls(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=20, restarts=1),
        null_repetitions=1,
        seed=0,
    )

    assert {summary.null_type for summary in summaries} == {
        "shuffled_sides",
        "random_constraints",
        "permuted_targets",
    }
