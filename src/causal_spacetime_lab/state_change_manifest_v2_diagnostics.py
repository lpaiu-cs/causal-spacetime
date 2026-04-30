"""Diagnostic-complete v2 manifest-family outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
    load_manifest_datasets,
)
from causal_spacetime_lab.state_change_manifest_null_taxonomy import classify_null_type
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationConfig,
    fit_manifest_dimension_curve,
    fit_manifest_ordinal_representation,
    representation_fit_to_row,
)
from causal_spacetime_lab.state_change_manifest_representation_nulls import (
    evaluate_manifest_representation_nulls,
)
from causal_spacetime_lab.state_change_manifest_representation_stability import (
    fit_manifest_restarts,
    heldout_violation_stability_summary,
    pairwise_latent_order_stability,
)


def _family_name(dataset: ManifestConstraintDataset) -> str:
    return str(dataset.manifest_json.get("profile_label", dataset.manifest_id))


def _family_kind(dataset: ManifestConstraintDataset) -> str:
    family = _family_name(dataset)
    if "failed" in family or not dataset.eligible:
        return "failed_control" if "failed" in family else "structured"
    return "structured"


def _finite_values(values: list[float]) -> list[float]:
    return [value for value in values if np.isfinite(value)]


def _mean(values: list[float]) -> float:
    finite = _finite_values(values)
    return float(np.mean(finite)) if finite else float("nan")


def _dataset_with_constraints(
    dataset: ManifestConstraintDataset,
    constraints: np.ndarray,
) -> ManifestConstraintDataset:
    array = np.asarray(constraints, dtype=int)
    return ManifestConstraintDataset(
        manifest_id=dataset.manifest_id,
        eligible=dataset.eligible,
        failed_reasons=dataset.failed_reasons,
        target_event_ids=dataset.target_event_ids,
        constraints=array,
        margins=dataset.margins,
        train_constraint_indices=dataset.train_constraint_indices,
        heldout_constraint_indices=dataset.heldout_constraint_indices,
        train_constraints=array[dataset.train_constraint_indices],
        heldout_constraints=array[dataset.heldout_constraint_indices],
        forbidden_interpretations=dataset.forbidden_interpretations,
        manifest_json=dataset.manifest_json,
    )


def _same_pair_marginal_constraint_null(
    constraints: np.ndarray,
    seed: int,
) -> np.ndarray:
    """Preserve left/right pair-node marginals while re-pairing constraints."""

    array = np.asarray(constraints, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("constraints must have shape (m, 4)")
    if array.shape[0] == 0:
        return array.copy()
    rng = np.random.default_rng(seed)
    left = array[:, 0:2].copy()
    right = array[:, 2:4].copy()
    right = right[rng.permutation(right.shape[0])]
    return np.hstack([left, right]).astype(int, copy=False)


def _same_marginal_null_row(
    dataset: ManifestConstraintDataset,
    config: ManifestRepresentationConfig,
    *,
    repetitions: int,
    seed: int,
) -> dict[str, float | str]:
    structured = fit_manifest_ordinal_representation(dataset, config)
    values: list[float] = []
    for repetition in range(repetitions):
        null_constraints = _same_pair_marginal_constraint_null(
            dataset.constraints,
            seed + repetition,
        )
        null_dataset = _dataset_with_constraints(dataset, null_constraints)
        null_config = ManifestRepresentationConfig(
            embedding_dim=config.embedding_dim,
            steps=config.steps,
            restarts=config.restarts,
            learning_rate=config.learning_rate,
            seed=config.seed + 50_000 + repetition,
        )
        fit = fit_manifest_ordinal_representation(null_dataset, null_config)
        values.append(float(fit.heldout_violation_rate))
    mean = _mean(values)
    finite = _finite_values(values)
    return {
        "manifest_id": dataset.manifest_id,
        "family_name": _family_name(dataset),
        "family_kind": _family_kind(dataset),
        "null_type": "random_same_marginals",
        "taxonomy_class": classify_null_type("random_same_marginals"),
        "embedding_dim": float(config.embedding_dim),
        "repetitions": float(repetitions),
        "mean_heldout_violation_rate": mean,
        "std_heldout_violation_rate": (
            float(np.std(finite)) if finite else float("nan")
        ),
        "best_heldout_violation_rate": (
            float(np.min(finite)) if finite else float("nan")
        ),
        "structured_heldout_violation_rate": float(
            structured.heldout_violation_rate
        ),
        "structured_minus_null_mean": float(structured.heldout_violation_rate) - mean,
        "taxonomy_note": "same_pair_marginal_constraint_null",
    }


def compute_v2_manifest_family_fit_rows(
    manifest_dir: Path,
    *,
    dims: list[int],
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Fit v2 manifests and return per-manifest dimension rows."""

    rows: list[dict[str, float | str]] = []
    datasets = load_manifest_datasets(manifest_dir, include_ineligible=True)
    for dataset_index, dataset in enumerate(datasets):
        base_config = ManifestRepresentationConfig(
            embedding_dim=int(dims[0]),
            steps=int(steps),
            restarts=int(restarts),
            seed=int(seed) + 1000 * dataset_index,
        )
        fits = fit_manifest_dimension_curve(dataset, dims, base_config)
        for fit in fits:
            row = representation_fit_to_row(fit)
            row["family_name"] = _family_name(dataset)
            row["family_kind"] = _family_kind(dataset)
            rows.append(row)
    return rows


def summarize_v2_fit_rows(
    fit_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Summarize v2 fit rows by family and latent dimension."""

    keys = sorted(
        {
            (
                str(row["family_name"]),
                str(row["family_kind"]),
                int(row["embedding_dim"]),
            )
            for row in fit_rows
        }
    )
    summaries: list[dict[str, float | str]] = []
    for family_name, family_kind, dim in keys:
        rows = [
            row
            for row in fit_rows
            if row["family_name"] == family_name
            and row["family_kind"] == family_kind
            and int(row["embedding_dim"]) == dim
        ]
        fitted = [row for row in rows if float(row["fitted"]) == 1.0]
        no_fit = len(rows) - len(fitted)
        train_values = [float(row["train_violation_rate"]) for row in fitted]
        heldout_values = [float(row["heldout_violation_rate"]) for row in fitted]
        summaries.append(
            {
                "family_name": family_name,
                "family_kind": family_kind,
                "embedding_dim": float(dim),
                "manifest_count": float(len(rows)),
                "fitted_count": float(len(fitted)),
                "no_fit_count": float(no_fit),
                "mean_train_violation": _mean(train_values),
                "mean_heldout_violation": _mean(heldout_values),
                "mean_generalization_gap": _mean(
                    [h - t for h, t in zip(heldout_values, train_values, strict=True)]
                ),
                "median_heldout_violation": (
                    float(np.median(_finite_values(heldout_values)))
                    if _finite_values(heldout_values)
                    else float("nan")
                ),
                "best_heldout_violation": min(_finite_values(heldout_values))
                if _finite_values(heldout_values)
                else float("nan"),
                "worst_heldout_violation": max(_finite_values(heldout_values))
                if _finite_values(heldout_values)
                else float("nan"),
            }
        )
    return summaries


def compute_v2_null_taxonomy_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    null_repetitions: int,
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute representation-null rows for v2 manifests."""

    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=False)
    ):
        config = ManifestRepresentationConfig(
            embedding_dim=int(embedding_dim),
            steps=int(steps),
            restarts=int(restarts),
            seed=int(seed) + 1000 * index,
        )
        for summary in evaluate_manifest_representation_nulls(
            dataset,
            config,
            null_repetitions=int(null_repetitions),
            seed=int(seed) + 10_000 * index,
        ):
            rows.append(
                {
                    "manifest_id": dataset.manifest_id,
                    "family_name": _family_name(dataset),
                    "family_kind": _family_kind(dataset),
                    "null_type": summary.null_type,
                    "taxonomy_class": classify_null_type(summary.null_type),
                    "embedding_dim": float(summary.embedding_dim),
                    "repetitions": float(summary.repetitions),
                    "mean_heldout_violation_rate": float(
                        summary.mean_heldout_violation_rate
                    ),
                    "std_heldout_violation_rate": float(
                        summary.std_heldout_violation_rate
                    ),
                    "best_heldout_violation_rate": float(
                        summary.best_heldout_violation_rate
                    ),
                    "structured_heldout_violation_rate": float(
                        summary.structured_heldout_violation_rate
                    ),
                    "structured_minus_null_mean": float(
                        summary.structured_minus_null_mean
                    ),
                    "taxonomy_note": (
                        "target_label_permutation_is_symmetry_control"
                        if summary.null_type == "permuted_targets"
                        else ""
                    ),
                }
            )
        rows.append(
            _same_marginal_null_row(
                dataset,
                config,
                repetitions=int(null_repetitions),
                seed=int(seed) + 20_000 * (index + 1),
            )
        )
    return rows


def compute_v2_stricter_criteria_rows(
    manifest_dir: Path,
    *,
    dims: list[int],
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Apply fixed diagnostic thresholds to v2 fit rows."""

    fit_rows = summarize_v2_fit_rows(
        compute_v2_manifest_family_fit_rows(
            manifest_dir,
            dims=dims,
            steps=steps,
            restarts=restarts,
            seed=seed,
        )
    )
    rows: list[dict[str, float | str]] = []
    for row in fit_rows:
        heldout = float(row["mean_heldout_violation"])
        gap = float(row["mean_generalization_gap"])
        passed = np.isfinite(heldout) and heldout <= 0.35
        passed = bool(passed and (not np.isfinite(gap) or gap <= 0.20))
        rows.append(
            {
                "family_name": row["family_name"],
                "family_kind": row["family_kind"],
                "embedding_dim": row["embedding_dim"],
                "threshold_pass": float(passed),
                "heldout_violation_threshold": 0.35,
                "generalization_gap_threshold": 0.20,
            }
        )
    return rows


def compute_v2_failed_accounting_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Report failed, ineligible, and no-fit v2 manifest accounting."""

    rows: list[dict[str, float | str]] = []
    for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True):
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name(dataset),
                "family_kind": _family_kind(dataset),
                "eligible": float(dataset.eligible),
                "failed_reasons": ";".join(dataset.failed_reasons),
                "row_type": "failed_accounting"
                if not dataset.eligible
                else "eligible_accounting",
            }
        )
    return rows


def compute_v2_restart_stability_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    restart_count: int,
    steps: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute v2 restart-stability diagnostics."""

    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=False)
    ):
        config = ManifestRepresentationConfig(
            embedding_dim=int(embedding_dim),
            steps=int(steps),
            restarts=1,
            seed=int(seed) + 1000 * index,
        )
        fits = fit_manifest_restarts(dataset, config, int(restart_count))
        summary = heldout_violation_stability_summary(fits)
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name(dataset),
                "family_kind": _family_kind(dataset),
                "embedding_dim": float(embedding_dim),
                "restart_count": float(restart_count),
                "restart_std": float(summary["std_heldout_violation_rate"]),
                **summary,
            }
        )
    return rows


def compute_v2_latent_order_stability_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    restart_count: int,
    steps: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute v2 latent-order stability diagnostics."""

    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=False)
    ):
        config = ManifestRepresentationConfig(
            embedding_dim=int(embedding_dim),
            steps=int(steps),
            restarts=1,
            seed=int(seed) + 1000 * index,
        )
        fits = fit_manifest_restarts(dataset, config, int(restart_count))
        summary = pairwise_latent_order_stability(
            fits,
            sample_pair_count=500,
            seed=int(seed) + 50_000 * index,
        )
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name(dataset),
                "family_kind": _family_kind(dataset),
                "embedding_dim": float(embedding_dim),
                "restart_count": float(restart_count),
                "latent_order_disagreement": float(
                    summary["mean_pair_order_disagreement"]
                ),
                **summary,
            }
        )
    return rows


def compute_v2_coverage_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Extract target and pair-node coverage metrics from v2 manifests."""

    rows: list[dict[str, float | str]] = []
    for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True):
        summary = dataset.manifest_json.get("validation_summary", {})
        if not isinstance(summary, dict):
            summary = {}
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name(dataset),
                "family_kind": _family_kind(dataset),
                "target_coverage_fraction": float(
                    summary.get("target_coverage_fraction", float("nan"))
                ),
                "pair_node_coverage_fraction": float(
                    summary.get("pair_node_coverage_fraction", float("nan"))
                ),
            }
        )
    return rows
