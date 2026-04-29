"""Latent ordinal representation fits for frozen handoff manifests.

The fitted arrays in this module are mathematical latent representation
variables for response-comparison constraints. They are not physical
coordinates and are not calibrated metric geometry.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.ordinal_embedding import (
    fit_ordinal_embedding_gradient_descent,
    quadruplet_hinge_loss,
    quadruplet_violation_rate,
)
from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
)


@dataclass(frozen=True)
class ManifestRepresentationConfig:
    """Configuration for a latent ordinal representation diagnostic."""

    embedding_dim: int
    steps: int = 800
    restarts: int = 3
    learning_rate: float = 0.05
    seed: int = 0


@dataclass(frozen=True)
class ManifestRepresentationFit:
    """Result of fitting one frozen manifest split.

    ``embedding`` is a latent ordinal representation, not physical coordinates.
    """

    manifest_id: str
    embedding_dim: int
    eligible: bool
    fitted: bool
    reason_not_fit: str
    train_violation_rate: float
    heldout_violation_rate: float
    train_hinge_loss: float
    heldout_hinge_loss: float
    train_constraint_count: int
    heldout_constraint_count: int
    target_count: int
    embedding: NDArray[np.float64] | None


def _empty_fit(
    dataset: ManifestConstraintDataset,
    config: ManifestRepresentationConfig,
    reason: str,
) -> ManifestRepresentationFit:
    return ManifestRepresentationFit(
        manifest_id=dataset.manifest_id,
        embedding_dim=int(config.embedding_dim),
        eligible=dataset.eligible,
        fitted=False,
        reason_not_fit=reason,
        train_violation_rate=float("nan"),
        heldout_violation_rate=float("nan"),
        train_hinge_loss=float("nan"),
        heldout_hinge_loss=float("nan"),
        train_constraint_count=int(dataset.train_constraints.shape[0]),
        heldout_constraint_count=int(dataset.heldout_constraints.shape[0]),
        target_count=int(dataset.target_event_ids.size),
        embedding=None,
    )


def fit_manifest_ordinal_representation(
    dataset: ManifestConstraintDataset,
    config: ManifestRepresentationConfig,
    *,
    fit_ineligible: bool = False,
) -> ManifestRepresentationFit:
    """Fit a latent ordinal representation using only manifest train constraints."""

    if config.embedding_dim <= 0:
        raise ValueError("embedding_dim must be positive")
    if config.steps <= 0:
        raise ValueError("steps must be positive")
    if config.restarts <= 0:
        raise ValueError("restarts must be positive")
    if config.learning_rate <= 0:
        raise ValueError("learning_rate must be positive")
    if not dataset.eligible and not fit_ineligible:
        return _empty_fit(dataset, config, "manifest_ineligible")
    if dataset.train_constraints.shape[0] == 0:
        return _empty_fit(dataset, config, "no_train_constraints")
    if dataset.target_event_ids.size < 2:
        return _empty_fit(dataset, config, "too_few_targets")

    embedding, _diagnostics = fit_ordinal_embedding_gradient_descent(
        int(dataset.target_event_ids.size),
        int(config.embedding_dim),
        dataset.train_constraints,
        steps=int(config.steps),
        learning_rate=float(config.learning_rate),
        seed=int(config.seed),
        restarts=int(config.restarts),
    )
    return ManifestRepresentationFit(
        manifest_id=dataset.manifest_id,
        embedding_dim=int(config.embedding_dim),
        eligible=dataset.eligible,
        fitted=True,
        reason_not_fit="",
        train_violation_rate=quadruplet_violation_rate(
            embedding,
            dataset.train_constraints,
        ),
        heldout_violation_rate=quadruplet_violation_rate(
            embedding,
            dataset.heldout_constraints,
        ),
        train_hinge_loss=quadruplet_hinge_loss(embedding, dataset.train_constraints),
        heldout_hinge_loss=quadruplet_hinge_loss(
            embedding,
            dataset.heldout_constraints,
        ),
        train_constraint_count=int(dataset.train_constraints.shape[0]),
        heldout_constraint_count=int(dataset.heldout_constraints.shape[0]),
        target_count=int(dataset.target_event_ids.size),
        embedding=embedding.astype(float, copy=True),
    )


def fit_manifest_dimension_curve(
    dataset: ManifestConstraintDataset,
    dims: list[int],
    base_config: ManifestRepresentationConfig,
    *,
    fit_ineligible: bool = False,
) -> list[ManifestRepresentationFit]:
    """Fit a manifest over candidate latent representation dimensions."""

    fits: list[ManifestRepresentationFit] = []
    for index, dim in enumerate(dims):
        config = ManifestRepresentationConfig(
            embedding_dim=int(dim),
            steps=base_config.steps,
            restarts=base_config.restarts,
            learning_rate=base_config.learning_rate,
            seed=base_config.seed + 10_000 * index,
        )
        fits.append(
            fit_manifest_ordinal_representation(
                dataset,
                config,
                fit_ineligible=fit_ineligible,
            )
        )
    return fits


def representation_fit_to_row(
    fit: ManifestRepresentationFit,
    *,
    include_embedding: bool = False,
) -> dict[str, float | str]:
    """Convert a representation fit to a CSV-safe row.

    Raw latent representation variables are excluded by default.
    """

    row: dict[str, float | str] = {
        "manifest_id": fit.manifest_id,
        "embedding_dim": float(fit.embedding_dim),
        "eligible": float(fit.eligible),
        "fitted": float(fit.fitted),
        "reason_not_fit": fit.reason_not_fit,
        "train_violation_rate": float(fit.train_violation_rate),
        "heldout_violation_rate": float(fit.heldout_violation_rate),
        "train_hinge_loss": float(fit.train_hinge_loss),
        "heldout_hinge_loss": float(fit.heldout_hinge_loss),
        "train_constraint_count": float(fit.train_constraint_count),
        "heldout_constraint_count": float(fit.heldout_constraint_count),
        "target_count": float(fit.target_count),
        "representation_kind": "latent_ordinal_representation",
    }
    if include_embedding:
        row["embedding_values"] = (
            "" if fit.embedding is None else np.array2string(fit.embedding)
        )
    return row
