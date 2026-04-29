"""Shared helpers for frozen-manifest representation experiments."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from response_handoff_experiment_helpers import deterministic_handoff_profile

from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
    load_manifest_datasets,
)
from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    write_handoff_manifest,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)


def write_csv(
    rows: list[dict[str, float | str]],
    path: Path,
    fallback_fields: list[str],
) -> Path:
    """Write rows to CSV with stable fallback fields."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else fallback_fields
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def ensure_manifests(
    manifest_dir: Path,
    *,
    output_dir: Path = Path("outputs"),
) -> list[ManifestConstraintDataset]:
    """Load manifests, exporting lightweight eligible manifests if needed."""

    datasets = load_manifest_datasets(manifest_dir, include_ineligible=True)
    if datasets:
        return datasets
    from exp136_response_handoff_manifest_export import (
        ExperimentConfig as ManifestExportConfig,
    )
    from exp136_response_handoff_manifest_export import (
        run_experiment as run_manifest_export,
    )

    run_manifest_export(
        ManifestExportConfig(
            reference_length=48,
            emission_positions=(6, 12, 18),
            layer_delay_ranks=(3, 5, 8),
            targets_per_layer=5,
            max_manifests=3,
            output_dir=output_dir,
        )
    )
    return load_manifest_datasets(manifest_dir, include_ineligible=True)


def build_exact_manifest(output_dir: Path, filename: str) -> Path:
    """Build a deterministic small handoff manifest for exact checks."""

    manifest = build_candidate_handoff_manifest(
        deterministic_handoff_profile(),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
        ConstraintValidationGate(
            "manifest_representation_exact",
            min_constraint_count=1,
            min_evaluable_fraction=0.0,
            min_agreement_fraction=0.0,
            max_inversion_fraction=1.0,
            max_tie_or_unresolved_fraction=1.0,
            min_null_z_score=-999.0,
            min_bootstrap_confidence=0.0,
        ),
        max_constraints=40,
        min_margin=0.0,
        train_fraction=0.6,
        constraint_seed=0,
        bootstrap_count=2,
        null_repetitions=1,
        source_label="manifest_representation_exact",
    )
    path = output_dir / "manifests" / filename
    write_handoff_manifest(manifest, path)
    return path


def save_simple_bar(
    labels: list[str],
    values: list[float],
    path: Path,
    *,
    ylabel: str,
) -> Path:
    """Save a compact bar chart."""

    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, values)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", labelrotation=35)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def save_grouped_line(
    rows: list[dict[str, float | str]],
    *,
    x_key: str,
    y_key: str,
    group_key: str,
    path: Path,
    ylabel: str,
) -> Path:
    """Save a grouped line plot for finite numeric values."""

    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    groups = sorted({str(row[group_key]) for row in rows})
    fig, ax = plt.subplots(figsize=(8, 4))
    for group in groups:
        group_rows = [row for row in rows if str(row[group_key]) == group]
        xs = np.asarray([float(row[x_key]) for row in group_rows], dtype=float)
        ys = np.asarray([float(row[y_key]) for row in group_rows], dtype=float)
        finite = np.isfinite(xs) & np.isfinite(ys)
        if not np.any(finite):
            continue
        order = np.argsort(xs[finite])
        ax.plot(xs[finite][order], ys[finite][order], marker="o", label=group[:12])
    ax.set_xlabel(x_key)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    if groups:
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path
