"""Restart stability diagnostics for latent manifest representations."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.ordinal_embedding import squared_distance_matrix
from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    ManifestRepresentationFit,
    fit_manifest_ordinal_representation,
)


def fit_manifest_restarts(
    dataset: ManifestConstraintDataset,
    config: ManifestRepresentationConfig,
    restart_count: int = 10,
) -> list[ManifestRepresentationFit]:
    """Fit a frozen manifest repeatedly with different random seeds."""

    if restart_count <= 0:
        raise ValueError("restart_count must be positive")
    fits: list[ManifestRepresentationFit] = []
    for restart in range(restart_count):
        restart_config = ManifestRepresentationConfig(
            embedding_dim=config.embedding_dim,
            steps=config.steps,
            restarts=config.restarts,
            learning_rate=config.learning_rate,
            seed=config.seed + 10_000 * restart,
        )
        fits.append(fit_manifest_ordinal_representation(dataset, restart_config))
    return fits


def heldout_violation_stability_summary(
    fits: list[ManifestRepresentationFit],
) -> dict[str, float]:
    """Summarize held-out violation variation across representation restarts."""

    values = np.asarray(
        [
            fit.heldout_violation_rate
            for fit in fits
            if fit.fitted and np.isfinite(fit.heldout_violation_rate)
        ],
        dtype=float,
    )
    if values.size == 0:
        return {
            "fit_count": 0.0,
            "mean_heldout_violation_rate": float("nan"),
            "std_heldout_violation_rate": float("nan"),
            "min_heldout_violation_rate": float("nan"),
            "max_heldout_violation_rate": float("nan"),
        }
    return {
        "fit_count": float(values.size),
        "mean_heldout_violation_rate": float(np.mean(values)),
        "std_heldout_violation_rate": float(np.std(values)),
        "min_heldout_violation_rate": float(np.min(values)),
        "max_heldout_violation_rate": float(np.max(values)),
    }


def _sample_pair_pair_indices(
    target_count: int,
    sample_pair_count: int,
    seed: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    if target_count < 3:
        return np.empty((0, 2), dtype=int), np.empty((0, 2), dtype=int)
    pairs = np.asarray(
        [(i, j) for i in range(target_count - 1) for j in range(i + 1, target_count)],
        dtype=int,
    )
    if pairs.shape[0] < 2:
        return np.empty((0, 2), dtype=int), np.empty((0, 2), dtype=int)
    rng = np.random.default_rng(seed)
    left = rng.integers(0, pairs.shape[0], size=sample_pair_count)
    right = rng.integers(0, pairs.shape[0], size=sample_pair_count)
    valid = left != right
    return pairs[left[valid]], pairs[right[valid]]


def _latent_pair_order_signs(
    embedding: np.ndarray,
    left_pairs: np.ndarray,
    right_pairs: np.ndarray,
) -> np.ndarray:
    distances = squared_distance_matrix(embedding)
    left = distances[left_pairs[:, 0], left_pairs[:, 1]]
    right = distances[right_pairs[:, 0], right_pairs[:, 1]]
    return np.sign(left - right)


def pairwise_latent_order_stability(
    fits: list[ManifestRepresentationFit],
    *,
    sample_pair_count: int = 1000,
    seed: int | None = None,
) -> dict[str, float]:
    """Compare latent pair-order signs across fitted representations.

    This is stability of fitted latent representation variables, not physical
    geometry.
    """

    fitted = [fit for fit in fits if fit.fitted and fit.embedding is not None]
    if len(fitted) < 2:
        return {
            "fit_pair_count": 0.0,
            "mean_pair_order_disagreement": float("nan"),
            "max_pair_order_disagreement": float("nan"),
        }
    if sample_pair_count <= 0:
        raise ValueError("sample_pair_count must be positive")
    left_pairs, right_pairs = _sample_pair_pair_indices(
        fitted[0].target_count,
        int(sample_pair_count),
        seed,
    )
    if left_pairs.size == 0:
        return {
            "fit_pair_count": 0.0,
            "mean_pair_order_disagreement": float("nan"),
            "max_pair_order_disagreement": float("nan"),
        }
    signs = [
        _latent_pair_order_signs(fit.embedding, left_pairs, right_pairs)
        for fit in fitted
        if fit.embedding is not None
    ]
    disagreements: list[float] = []
    for i in range(len(signs) - 1):
        for j in range(i + 1, len(signs)):
            comparable = (signs[i] != 0.0) | (signs[j] != 0.0)
            if np.any(comparable):
                disagreements.append(
                    float(np.mean(signs[i][comparable] != signs[j][comparable]))
                )
    values = np.asarray(disagreements, dtype=float)
    return {
        "fit_pair_count": float(values.size),
        "mean_pair_order_disagreement": (
            float(np.mean(values)) if values.size else float("nan")
        ),
        "max_pair_order_disagreement": (
            float(np.max(values)) if values.size else float("nan")
        ),
    }
