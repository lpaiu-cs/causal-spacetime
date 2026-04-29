from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.abspath("experiments"))

from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_dataset
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
)
from causal_spacetime_lab.state_change_manifest_representation_stability import (
    fit_manifest_restarts,
    heldout_violation_stability_summary,
    pairwise_latent_order_stability,
)


def test_fit_manifest_restarts_returns_requested_count(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)

    fits = fit_manifest_restarts(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=20, restarts=1),
        restart_count=3,
    )

    assert len(fits) == 3
    assert all(fit.fitted for fit in fits)


def test_heldout_violation_stability_summary_fields(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)
    fits = fit_manifest_restarts(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=20, restarts=1),
        restart_count=2,
    )

    report = heldout_violation_stability_summary(fits)

    assert report["fit_count"] == 2.0
    assert "mean_heldout_violation_rate" in report


def test_pairwise_latent_order_stability_fields(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = build_exact_manifest(tmp_path, "manifest.json")
    dataset = load_manifest_dataset(path)
    fits = fit_manifest_restarts(
        dataset,
        ManifestRepresentationConfig(embedding_dim=1, steps=20, restarts=1),
        restart_count=2,
    )

    report = pairwise_latent_order_stability(fits, sample_pair_count=20, seed=0)

    assert "mean_pair_order_disagreement" in report
