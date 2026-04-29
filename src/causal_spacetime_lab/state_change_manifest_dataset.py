"""Frozen handoff manifest loading for representation diagnostics.

Milestone 32 consumes handoff manifests as frozen input specifications. A
loaded dataset exposes response-comparison constraints and the manifest-defined
train/held-out split. It must not contain fitted representation coordinates or
metric interpretation fields.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class ManifestConstraintDataset:
    """Validated response-comparison constraints from one frozen manifest."""

    manifest_id: str
    eligible: bool
    failed_reasons: list[str]
    target_event_ids: NDArray[np.int_]
    constraints: NDArray[np.int_]
    margins: NDArray[np.float64]
    train_constraint_indices: NDArray[np.int_]
    heldout_constraint_indices: NDArray[np.int_]
    train_constraints: NDArray[np.int_]
    heldout_constraints: NDArray[np.int_]
    forbidden_interpretations: list[str]
    manifest_json: dict[str, object]


def forbidden_manifest_fields() -> list[str]:
    """Return manifest field-name fragments forbidden at fit time."""

    return [
        "embedding",
        "coordinates",
        "metric",
        "spatial_distance",
        "radar_distance",
        "geometry",
    ]


def _field_names_containing(
    value: Any,
    forbidden_fragments: list[str],
    *,
    prefix: str = "",
) -> list[str]:
    hits: list[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            key_text = str(key)
            path = f"{prefix}.{key_text}" if prefix else key_text
            lower = key_text.lower()
            if any(fragment in lower for fragment in forbidden_fragments):
                hits.append(path)
            hits.extend(
                _field_names_containing(
                    item,
                    forbidden_fragments,
                    prefix=path,
                )
            )
    elif isinstance(value, list):
        for index, item in enumerate(value):
            hits.extend(
                _field_names_containing(
                    item,
                    forbidden_fragments,
                    prefix=f"{prefix}[{index}]",
                )
            )
    return hits


def _as_constraint_array(value: object) -> NDArray[np.int_]:
    array = np.asarray(value, dtype=int)
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError("manifest constraints must have shape (m, 4)")
    if np.any(array < 0):
        raise IndexError("manifest constraints contain negative target indices")
    return array


def _as_index_array(value: object, name: str) -> NDArray[np.int_]:
    array = np.asarray(value, dtype=int)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if np.any(array < 0):
        raise IndexError(f"{name} contains negative indices")
    return array


def _has_embedding_fields(manifest_json: dict[str, object]) -> bool:
    return bool(_field_names_containing(manifest_json, ["embedding", "coordinates"]))


def _has_metric_fields(manifest_json: dict[str, object]) -> bool:
    return bool(
        _field_names_containing(
            manifest_json,
            ["metric", "spatial_distance", "radar_distance", "geometry"],
        )
    )


def load_manifest_dataset(path: Path) -> ManifestConstraintDataset:
    """Load and validate a frozen response-comparison handoff manifest."""

    manifest_json = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(manifest_json, dict):
        raise ValueError("manifest JSON must be an object")
    forbidden_field_hits = _field_names_containing(
        manifest_json,
        forbidden_manifest_fields(),
    )
    if forbidden_field_hits:
        raise ValueError(
            "manifest contains forbidden fit or interpretation fields: "
            + ", ".join(forbidden_field_hits)
        )

    target_event_ids = _as_index_array(
        manifest_json.get("target_event_ids", []),
        "target_event_ids",
    )
    if target_event_ids.size == 0:
        raise ValueError("target_event_ids must be nonempty")
    constraints = _as_constraint_array(manifest_json.get("constraints", []))
    margins = np.asarray(manifest_json.get("margins", []), dtype=float)
    if margins.ndim != 1 or margins.shape[0] != constraints.shape[0]:
        raise ValueError("margins must be one-dimensional and match constraints")

    train_indices = _as_index_array(
        manifest_json.get("train_constraint_indices", []),
        "train_constraint_indices",
    )
    heldout_indices = _as_index_array(
        manifest_json.get("heldout_constraint_indices", []),
        "heldout_constraint_indices",
    )
    constraint_count = constraints.shape[0]
    if train_indices.size and int(np.max(train_indices)) >= constraint_count:
        raise IndexError("train_constraint_indices contain an out-of-range index")
    if heldout_indices.size and int(np.max(heldout_indices)) >= constraint_count:
        raise IndexError("heldout_constraint_indices contain an out-of-range index")
    if set(train_indices.tolist()) & set(heldout_indices.tolist()):
        raise ValueError("train and heldout constraint indices must be disjoint")
    if constraints.size and int(np.max(constraints)) >= target_event_ids.size:
        raise IndexError("constraints reference an invalid target row index")

    forbidden_interpretations = list(
        manifest_json.get("forbidden_interpretations", [])
    )
    if not forbidden_interpretations:
        raise ValueError("forbidden_interpretations must be nonempty")

    decision = manifest_json.get("handoff_decision", {})
    if not isinstance(decision, dict):
        raise ValueError("handoff_decision must be an object")
    eligible = bool(decision.get("eligible", False))
    failed_reasons = [str(reason) for reason in decision.get("failed_reasons", [])]

    return ManifestConstraintDataset(
        manifest_id=str(manifest_json.get("manifest_id", path.stem)),
        eligible=eligible,
        failed_reasons=failed_reasons,
        target_event_ids=target_event_ids.astype(int, copy=False),
        constraints=constraints.astype(int, copy=False),
        margins=margins.astype(float, copy=False),
        train_constraint_indices=train_indices.astype(int, copy=False),
        heldout_constraint_indices=heldout_indices.astype(int, copy=False),
        train_constraints=constraints[train_indices].astype(int, copy=True),
        heldout_constraints=constraints[heldout_indices].astype(int, copy=True),
        forbidden_interpretations=forbidden_interpretations,
        manifest_json=manifest_json,
    )


def load_manifest_datasets(
    manifest_dir: Path,
    *,
    include_ineligible: bool = True,
) -> list[ManifestConstraintDataset]:
    """Load all manifest JSON files from a directory."""

    if not manifest_dir.exists():
        return []
    datasets: list[ManifestConstraintDataset] = []
    for path in sorted(manifest_dir.glob("*.json")):
        dataset = load_manifest_dataset(path)
        if include_ineligible or dataset.eligible:
            datasets.append(dataset)
    return datasets


def manifest_integrity_report(
    dataset: ManifestConstraintDataset,
) -> dict[str, float | str]:
    """Return integrity diagnostics for a loaded frozen manifest."""

    train = set(dataset.train_constraint_indices.tolist())
    heldout = set(dataset.heldout_constraint_indices.tolist())
    train_heldout_disjoint = not bool(train & heldout)
    has_embedding_fields = _has_embedding_fields(dataset.manifest_json)
    has_metric_fields = _has_metric_fields(dataset.manifest_json)
    passed = (
        train_heldout_disjoint
        and not has_embedding_fields
        and not has_metric_fields
        and bool(dataset.forbidden_interpretations)
    )
    return {
        "manifest_id": dataset.manifest_id,
        "eligible": float(dataset.eligible),
        "constraint_count": float(dataset.constraints.shape[0]),
        "train_constraint_count": float(dataset.train_constraints.shape[0]),
        "heldout_constraint_count": float(dataset.heldout_constraints.shape[0]),
        "target_count": float(dataset.target_event_ids.size),
        "train_heldout_disjoint": float(train_heldout_disjoint),
        "forbidden_interpretation_count": float(
            len(dataset.forbidden_interpretations)
        ),
        "has_embedding_fields": float(has_embedding_fields),
        "has_metric_fields": float(has_metric_fields),
        "passed_integrity": float(passed),
    }
