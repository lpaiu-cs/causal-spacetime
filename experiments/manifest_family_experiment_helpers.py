"""Shared helpers for frozen-manifest family comparison experiments."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from manifest_representation_experiment_helpers import (
    build_exact_manifest,
    ensure_manifests,
)

from causal_spacetime_lab.state_change_manifest_dataset import (
    ManifestConstraintDataset,
    load_manifest_dataset,
    load_manifest_datasets,
)
from causal_spacetime_lab.state_change_manifest_family import (
    ManifestFamilyAssignment,
    assign_manifest_families,
    default_manifest_family_specs,
)


def write_csv(
    rows: list[dict[str, float | str]],
    path: Path,
    fallback_fields: list[str],
) -> Path:
    """Write CSV rows with fallback fields for empty outputs."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else fallback_fields
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def load_or_create_family_datasets(
    manifest_dir: Path,
    output_dir: Path,
) -> list[ManifestConstraintDataset]:
    """Load manifests, creating lightweight eligible manifests if missing."""

    datasets = ensure_manifests(manifest_dir, output_dir=output_dir)
    if datasets:
        return datasets
    return load_manifest_datasets(manifest_dir, include_ineligible=True)


def ensure_failed_control_manifest(output_dir: Path) -> ManifestConstraintDataset:
    """Create and load a deterministic failed manifest control."""

    base_path = build_exact_manifest(
        output_dir,
        "response_handoff_family_failed_control_base.json",
    )
    payload = json.loads(base_path.read_text(encoding="utf-8"))
    payload["manifest_id"] = "family_failed_control_manifest"
    payload["profile_label"] = "family_failed_control"
    payload["handoff_decision"]["eligible"] = False
    payload["handoff_decision"]["failed_reasons"] = [
        "synthetic_family_failed_control",
    ]
    failed_path = output_dir / "manifests" / "response_handoff_family_failed.json"
    failed_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return load_manifest_dataset(failed_path)


def assignments_for_datasets(
    datasets: list[ManifestConstraintDataset],
) -> list[ManifestFamilyAssignment]:
    """Assign datasets to default preregistered manifest families."""

    return assign_manifest_families(datasets, default_manifest_family_specs())


def assignment_lookup(
    assignments: list[ManifestFamilyAssignment],
) -> dict[str, ManifestFamilyAssignment]:
    """Return assignment lookup by manifest id."""

    return {assignment.manifest_id: assignment for assignment in assignments}


def save_bar_figure(
    labels: list[str],
    values: list[float],
    path: Path,
    *,
    ylabel: str,
) -> Path:
    """Save a simple bar figure."""

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


def save_family_metric_figure(
    rows: list[dict[str, float | str]],
    *,
    metric: str,
    path: Path,
    ylabel: str,
) -> Path:
    """Save mean family metric by family name."""

    grouped: dict[str, list[float]] = {}
    for row in rows:
        value = float(row[metric])
        if not np.isfinite(value):
            continue
        grouped.setdefault(str(row["family_name"]), []).append(value)
    labels = sorted(grouped)
    values = [float(np.mean(grouped[label])) for label in labels]
    return save_bar_figure(labels, values, path, ylabel=ylabel)
