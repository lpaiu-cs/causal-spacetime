"""Diagnostics for protocol-invariant, provenance-aware v3 manifests."""

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
from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_coverage_rows,
    compute_v2_failed_accounting_rows,
    summarize_v2_fit_rows,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    MeasurementProtocolSpec,
    parameter_complete_for_protocol,
)
from causal_spacetime_lab.state_change_response_handoff import (
    manifest_provenance_integrity_report,
)


def compute_v3_protocol_manifest_family_fit_rows(
    manifest_dir: Path,
    *,
    dims: list[int],
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Fit protocol-bearing v3 manifests and return per-manifest rows."""

    _validate_manifest_directory(manifest_dir)
    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=True)
    ):
        family_kind = _family_kind_from_manifest(dataset.manifest_json)
        fit_dataset = (
            _eligible_copy(dataset) if family_kind == "structured" else dataset
        )
        base_config = ManifestRepresentationConfig(
            embedding_dim=int(dims[0]),
            steps=int(steps),
            restarts=int(restarts),
            seed=int(seed) + 1000 * index,
        )
        fits = fit_manifest_dimension_curve(
            fit_dataset,
            dims,
            base_config,
            fit_ineligible=family_kind == "structured",
        )
        for fit in fits:
            row = representation_fit_to_row(fit)
            row["family_name"] = _family_name_from_manifest(dataset.manifest_json)
            row["family_kind"] = family_kind
            rows.append(row)
    return rows


def summarize_v3_protocol_fit_rows(
    fit_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Summarize protocol-bearing v3 fit rows by family and dimension."""

    return summarize_v2_fit_rows(fit_rows)


def compute_v3_protocol_null_taxonomy_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    null_repetitions: int,
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute null taxonomy rows for protocol-bearing v3 manifests."""

    _validate_manifest_directory(manifest_dir)
    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=True)
    ):
        if _family_kind_from_manifest(dataset.manifest_json) != "structured":
            continue
        fit_dataset = _eligible_copy(dataset)
        config = ManifestRepresentationConfig(
            embedding_dim=int(embedding_dim),
            steps=int(steps),
            restarts=int(restarts),
            seed=int(seed) + 1000 * index,
        )
        for summary in evaluate_manifest_representation_nulls(
            fit_dataset,
            config,
            null_repetitions=int(null_repetitions),
            seed=int(seed) + 10_000 * index,
        ):
            rows.append(
                {
                    "manifest_id": dataset.manifest_id,
                    "family_name": _family_name_from_manifest(dataset.manifest_json),
                    "family_kind": "structured",
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
                fit_dataset, config, int(null_repetitions), int(seed) + index
            )
        )
    return rows


def compute_v3_protocol_stricter_criteria_rows(
    manifest_dir: Path,
    *,
    dims: list[int],
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute stricter diagnostic rows without carry-forward decisions."""

    _validate_manifest_directory(manifest_dir)
    fit_rows = summarize_v3_protocol_fit_rows(
        compute_v3_protocol_manifest_family_fit_rows(
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


def compute_v3_protocol_failed_accounting_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Report failed/report-only accounting for protocol-bearing v3 manifests."""

    _validate_manifest_directory(manifest_dir)
    rows = _rewrite_family_rows(
        manifest_dir,
        compute_v2_failed_accounting_rows(manifest_dir),
    )
    provenance = {
        dataset.manifest_id: str(
            dataset.manifest_json.get("handoff_provenance_type", "")
        )
        for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True)
    }
    for row in rows:
        row["handoff_provenance_type"] = provenance.get(str(row["manifest_id"]), "")
    return rows


def compute_v3_protocol_coverage_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Extract v3 target and pair-node coverage rows."""

    _validate_manifest_directory(manifest_dir)
    return _rewrite_family_rows(manifest_dir, compute_v2_coverage_rows(manifest_dir))


def compute_v3_protocol_restart_stability_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    restart_count: int,
    steps: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute v3 restart-stability rows."""

    _validate_manifest_directory(manifest_dir)
    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=True)
    ):
        if _family_kind_from_manifest(dataset.manifest_json) != "structured":
            continue
        config = ManifestRepresentationConfig(
            embedding_dim=int(embedding_dim),
            steps=int(steps),
            restarts=1,
            seed=int(seed) + 1000 * index,
        )
        fits = fit_manifest_restarts(
            _eligible_copy(dataset), config, int(restart_count)
        )
        summary = heldout_violation_stability_summary(fits)
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name_from_manifest(dataset.manifest_json),
                "family_kind": "structured",
                "embedding_dim": float(embedding_dim),
                "restart_count": float(restart_count),
                "restart_std": float(summary["std_heldout_violation_rate"]),
                **summary,
            }
        )
    return rows


def compute_v3_protocol_latent_order_stability_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    restart_count: int,
    steps: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute v3 latent-order stability rows."""

    _validate_manifest_directory(manifest_dir)
    rows: list[dict[str, float | str]] = []
    for index, dataset in enumerate(
        load_manifest_datasets(manifest_dir, include_ineligible=True)
    ):
        if _family_kind_from_manifest(dataset.manifest_json) != "structured":
            continue
        config = ManifestRepresentationConfig(
            embedding_dim=int(embedding_dim),
            steps=int(steps),
            restarts=1,
            seed=int(seed) + 1000 * index,
        )
        fits = fit_manifest_restarts(
            _eligible_copy(dataset), config, int(restart_count)
        )
        summary = pairwise_latent_order_stability(
            fits,
            sample_pair_count=500,
            seed=int(seed) + 50_000 * index,
        )
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name_from_manifest(dataset.manifest_json),
                "family_kind": "structured",
                "embedding_dim": float(embedding_dim),
                "restart_count": float(restart_count),
                "latent_order_disagreement": float(
                    summary["mean_pair_order_disagreement"]
                ),
                **summary,
            }
        )
    return rows


def manifest_metadata_audit_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Audit protocol/profile/provenance metadata on v3 manifests."""

    rows: list[dict[str, float | str]] = []
    for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True):
        manifest = dataset.manifest_json
        measurement = manifest.get("measurement_protocol")
        protocol_complete = 0.0
        if isinstance(measurement, dict):
            protocol = _protocol_from_manifest(measurement)
            protocol_complete = float(parameter_complete_for_protocol(protocol))
        provenance_report = manifest_provenance_integrity_report(manifest)
        profile_status = str(manifest.get("profile_invariance_status", ""))
        family_kind = _family_kind_from_manifest(manifest)
        structured = family_kind == "structured"
        rows.append(
            {
                "manifest_id": dataset.manifest_id,
                "family_name": _family_name_from_manifest(manifest),
                "family_kind": family_kind,
                "has_profile_metadata": float(
                    isinstance(manifest.get("profile_metadata"), dict)
                ),
                "has_measurement_protocol_id": float(
                    bool(manifest.get("measurement_protocol_id"))
                ),
                "has_measurement_protocol_hash": float(
                    bool(manifest.get("measurement_protocol_hash"))
                ),
                "has_handoff_provenance": float(
                    isinstance(manifest.get("handoff_provenance"), dict)
                ),
                "has_handoff_provenance_type": float(
                    bool(manifest.get("handoff_provenance_type"))
                ),
                "profile_invariance_status": profile_status,
                "parameter_complete": protocol_complete,
                "admissible_for_pairwise_dissimilarity": float(
                    bool(manifest.get("admissible_for_pairwise_dissimilarity", False))
                ),
                "valid_provenance": provenance_report["provenance_integrity_passed"],
                "structured_metadata_pass": float(
                    (not structured)
                    or (
                        profile_status == "protocol_invariant"
                        and protocol_complete == 1.0
                        and bool(
                            manifest.get(
                                "admissible_for_pairwise_dissimilarity",
                                False,
                            )
                        )
                        and float(provenance_report["provenance_integrity_passed"])
                        == 1.0
                    )
                ),
            }
        )
    return rows


def _validate_manifest_directory(manifest_dir: Path) -> None:
    for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True):
        manifest = dataset.manifest_json
        if not isinstance(manifest.get("measurement_protocol"), dict):
            raise ValueError("v3 protocol manifest missing measurement_protocol")
        if not isinstance(manifest.get("profile_metadata"), dict):
            raise ValueError("v3 protocol manifest missing profile_metadata")
        if not isinstance(manifest.get("handoff_provenance"), dict):
            raise ValueError("v3 protocol manifest missing handoff_provenance")
        if _family_kind_from_manifest(manifest) == "structured":
            if manifest.get("profile_invariance_status") != "protocol_invariant":
                raise ValueError("structured v3 manifest is not protocol_invariant")
            if not bool(manifest.get("admissible_for_pairwise_dissimilarity", False)):
                raise ValueError("structured v3 manifest is not admissible")
            protocol = _protocol_from_manifest(manifest["measurement_protocol"])
            if not parameter_complete_for_protocol(protocol):
                raise ValueError("structured v3 manifest has incomplete parameters")
            report = manifest_provenance_integrity_report(manifest)
            if float(report["provenance_integrity_passed"]) != 1.0:
                raise ValueError("structured v3 manifest has invalid provenance")


def _rewrite_family_rows(
    manifest_dir: Path,
    rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    lookup = {
        dataset.manifest_id: (
            _family_name_from_manifest(dataset.manifest_json),
            _family_kind_from_manifest(dataset.manifest_json),
        )
        for dataset in load_manifest_datasets(manifest_dir, include_ineligible=True)
    }
    rewritten: list[dict[str, float | str]] = []
    for row in rows:
        manifest_id = str(row.get("manifest_id", ""))
        family_name, family_kind = lookup.get(
            manifest_id,
            (str(row.get("family_name", "")), str(row.get("family_kind", ""))),
        )
        updated = dict(row)
        updated["family_name"] = family_name
        updated["family_kind"] = family_kind
        rewritten.append(updated)
    return rewritten


def _eligible_copy(dataset: ManifestConstraintDataset) -> ManifestConstraintDataset:
    return ManifestConstraintDataset(
        manifest_id=dataset.manifest_id,
        eligible=True,
        failed_reasons=dataset.failed_reasons,
        target_event_ids=dataset.target_event_ids,
        constraints=dataset.constraints,
        margins=dataset.margins,
        train_constraint_indices=dataset.train_constraint_indices,
        heldout_constraint_indices=dataset.heldout_constraint_indices,
        train_constraints=dataset.train_constraints,
        heldout_constraints=dataset.heldout_constraints,
        forbidden_interpretations=dataset.forbidden_interpretations,
        manifest_json=dataset.manifest_json,
    )


def _dataset_with_constraints(
    dataset: ManifestConstraintDataset,
    constraints: np.ndarray,
) -> ManifestConstraintDataset:
    array = np.asarray(constraints, dtype=int)
    return ManifestConstraintDataset(
        manifest_id=dataset.manifest_id,
        eligible=True,
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


def _same_marginal_constraints(constraints: np.ndarray, seed: int) -> np.ndarray:
    array = np.asarray(constraints, dtype=int)
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
    repetitions: int,
    seed: int,
) -> dict[str, float | str]:
    structured = fit_manifest_ordinal_representation(
        dataset,
        config,
        fit_ineligible=True,
    )
    values: list[float] = []
    for repetition in range(repetitions):
        null_dataset = _dataset_with_constraints(
            dataset,
            _same_marginal_constraints(dataset.constraints, seed + repetition),
        )
        fit = fit_manifest_ordinal_representation(
            null_dataset,
            ManifestRepresentationConfig(
                embedding_dim=config.embedding_dim,
                steps=config.steps,
                restarts=config.restarts,
                learning_rate=config.learning_rate,
                seed=config.seed + 80_000 + repetition,
            ),
            fit_ineligible=True,
        )
        values.append(float(fit.heldout_violation_rate))
    finite = [value for value in values if np.isfinite(value)]
    mean = float(np.mean(finite)) if finite else float("nan")
    return {
        "manifest_id": dataset.manifest_id,
        "family_name": _family_name_from_manifest(dataset.manifest_json),
        "family_kind": "structured",
        "null_type": "random_same_marginals",
        "taxonomy_class": classify_null_type("random_same_marginals"),
        "embedding_dim": float(config.embedding_dim),
        "repetitions": float(repetitions),
        "mean_heldout_violation_rate": mean,
        "std_heldout_violation_rate": float(np.std(finite)) if finite else float("nan"),
        "best_heldout_violation_rate": min(finite) if finite else float("nan"),
        "structured_heldout_violation_rate": float(structured.heldout_violation_rate),
        "structured_minus_null_mean": float(structured.heldout_violation_rate) - mean,
        "taxonomy_note": "same_pair_marginal_constraint_null",
    }


def _family_name_from_manifest(manifest: dict[str, object]) -> str:
    template = manifest.get("top_down_template")
    if isinstance(template, dict) and template.get("family_name"):
        return str(template["family_name"])
    profile_label = str(manifest.get("profile_label", ""))
    if "_m" in profile_label:
        return profile_label.rsplit("_m", 1)[0]
    return profile_label


def _family_kind_from_manifest(manifest: dict[str, object]) -> str:
    provenance_type = str(manifest.get("handoff_provenance_type", ""))
    family = _family_name_from_manifest(manifest)
    if provenance_type == "report_only_control":
        return "failed_control" if "failed" in family else "report_only"
    return "structured"


def _protocol_from_manifest(payload: object) -> MeasurementProtocolSpec:
    data = dict(payload) if isinstance(payload, dict) else {}
    params = data.get("protocol_parameters")
    return MeasurementProtocolSpec(
        echo_rule=str(data["echo_rule"]),
        emission_rule=str(data["emission_rule"]),
        gate_rule=str(data["gate_rule"]),
        reference_subsampling_rule=str(data["reference_subsampling_rule"]),
        spectrum_type=str(data["spectrum_type"]),
        normalization_rule=str(data["normalization_rule"]),
        missing_policy=str(data["missing_policy"]),
        tie_policy=str(data["tie_policy"]),
        margin_policy=str(data["margin_policy"]),
        protocol_label=str(data.get("protocol_label", "")),
        protocol_parameters=dict(params) if isinstance(params, dict) else {},
    )
