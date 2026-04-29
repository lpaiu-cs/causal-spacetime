"""Null baselines for frozen-manifest latent representation diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_ordinal_representation,
)


def shuffle_constraint_sides_for_manifest(
    constraints: NDArray[np.int_],
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Randomly swap left/right response-comparison pairs in constraints."""

    shuffled = np.asarray(constraints, dtype=int).copy()
    if shuffled.ndim != 2 or shuffled.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    rng = np.random.default_rng(seed)
    swap = rng.random(shuffled.shape[0]) < 0.5
    left = shuffled[swap, 0:2].copy()
    shuffled[swap, 0:2] = shuffled[swap, 2:4]
    shuffled[swap, 2:4] = left
    return shuffled


def random_manifest_constraints(
    target_count: int,
    constraint_count: int,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Generate random valid response-comparison quadruplets."""

    n = int(target_count)
    count = int(constraint_count)
    if n < 2:
        raise ValueError("target_count must be at least 2")
    if count < 0:
        raise ValueError("constraint_count must be nonnegative")
    rng = np.random.default_rng(seed)
    rows = np.empty((count, 4), dtype=int)
    accepted = 0
    while accepted < count:
        batch = max(256, 2 * (count - accepted))
        candidate = rng.integers(0, n, size=(batch, 4))
        left_ok = candidate[:, 0] != candidate[:, 1]
        right_ok = candidate[:, 2] != candidate[:, 3]
        left_pair = np.sort(candidate[:, 0:2], axis=1)
        right_pair = np.sort(candidate[:, 2:4], axis=1)
        different_pairs = np.any(left_pair != right_pair, axis=1)
        valid = candidate[left_ok & right_ok & different_pairs]
        take = min(valid.shape[0], count - accepted)
        if take:
            rows[accepted : accepted + take] = valid[:take]
            accepted += take
    return rows


def permute_manifest_target_labels(
    constraints: NDArray[np.int_],
    target_count: int,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Apply a random target-row permutation to all constraints."""

    constraint_array = np.asarray(constraints, dtype=int)
    if constraint_array.ndim != 2 or constraint_array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    n = int(target_count)
    if n <= 0:
        raise ValueError("target_count must be positive")
    if constraint_array.size and int(np.max(constraint_array)) >= n:
        raise IndexError("constraints contain an out-of-range target index")
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(n)
    return permutation[constraint_array].astype(int, copy=False)


@dataclass(frozen=True)
class ManifestNullFitSummary:
    """Held-out fit summary for one manifest representation null type."""

    manifest_id: str
    null_type: str
    embedding_dim: int
    repetitions: int
    mean_heldout_violation_rate: float
    std_heldout_violation_rate: float
    best_heldout_violation_rate: float
    structured_heldout_violation_rate: float
    structured_minus_null_mean: float


def _dataset_with_constraints(
    dataset: ManifestConstraintDataset,
    constraints: NDArray[np.int_],
) -> ManifestConstraintDataset:
    return ManifestConstraintDataset(
        manifest_id=dataset.manifest_id,
        eligible=dataset.eligible,
        failed_reasons=dataset.failed_reasons,
        target_event_ids=dataset.target_event_ids,
        constraints=np.asarray(constraints, dtype=int),
        margins=dataset.margins,
        train_constraint_indices=dataset.train_constraint_indices,
        heldout_constraint_indices=dataset.heldout_constraint_indices,
        train_constraints=np.asarray(constraints, dtype=int)[
            dataset.train_constraint_indices
        ],
        heldout_constraints=np.asarray(constraints, dtype=int)[
            dataset.heldout_constraint_indices
        ],
        forbidden_interpretations=dataset.forbidden_interpretations,
        manifest_json=dataset.manifest_json,
    )


def _null_constraints(
    dataset: ManifestConstraintDataset,
    null_type: str,
    seed: int,
) -> NDArray[np.int_]:
    if null_type == "shuffled_sides":
        return shuffle_constraint_sides_for_manifest(dataset.constraints, seed=seed)
    if null_type == "random_constraints":
        return random_manifest_constraints(
            int(dataset.target_event_ids.size),
            int(dataset.constraints.shape[0]),
            seed=seed,
        )
    if null_type == "permuted_targets":
        return permute_manifest_target_labels(
            dataset.constraints,
            int(dataset.target_event_ids.size),
            seed=seed,
        )
    raise ValueError(f"unknown null_type: {null_type}")


def evaluate_manifest_representation_nulls(
    dataset: ManifestConstraintDataset,
    config: ManifestRepresentationConfig,
    *,
    null_repetitions: int = 10,
    seed: int = 0,
) -> list[ManifestNullFitSummary]:
    """Evaluate representation nulls for one frozen manifest.

    Nulls test representation artifacts. Better-than-null fit is not evidence
    of physical geometry.
    """

    if null_repetitions <= 0:
        raise ValueError("null_repetitions must be positive")
    structured = fit_manifest_ordinal_representation(dataset, config)
    structured_heldout = float(structured.heldout_violation_rate)
    rows: list[ManifestNullFitSummary] = []
    for type_index, null_type in enumerate(
        ["shuffled_sides", "random_constraints", "permuted_targets"]
    ):
        values: list[float] = []
        for repetition in range(null_repetitions):
            null_constraints = _null_constraints(
                dataset,
                null_type,
                seed + 10_000 * type_index + repetition,
            )
            null_dataset = _dataset_with_constraints(dataset, null_constraints)
            null_config = ManifestRepresentationConfig(
                embedding_dim=config.embedding_dim,
                steps=config.steps,
                restarts=config.restarts,
                learning_rate=config.learning_rate,
                seed=config.seed + 100_000 * (type_index + 1) + repetition,
            )
            fit = fit_manifest_ordinal_representation(null_dataset, null_config)
            values.append(float(fit.heldout_violation_rate))
        finite = np.asarray([value for value in values if np.isfinite(value)])
        mean = float(np.mean(finite)) if finite.size else float("nan")
        rows.append(
            ManifestNullFitSummary(
                manifest_id=dataset.manifest_id,
                null_type=null_type,
                embedding_dim=int(config.embedding_dim),
                repetitions=int(null_repetitions),
                mean_heldout_violation_rate=mean,
                std_heldout_violation_rate=(
                    float(np.std(finite)) if finite.size else float("nan")
                ),
                best_heldout_violation_rate=(
                    float(np.min(finite)) if finite.size else float("nan")
                ),
                structured_heldout_violation_rate=structured_heldout,
                structured_minus_null_mean=structured_heldout - mean,
            )
        )
    return rows
