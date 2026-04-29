"""Predeclared handoff manifests for response-comparison constraints."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)


@dataclass(frozen=True)
class HandoffValidationSummary:
    """Summary metrics used to decide handoff eligibility."""

    heldout_agreement_fraction: float
    heldout_inversion_fraction: float
    heldout_evaluable_fraction: float
    bootstrap_mean_agreement_fraction: float
    bootstrap_stable_constraint_fraction: float
    null_z_score: float
    best_null_agreement_fraction: float
    target_coverage_fraction: float
    pair_node_coverage_fraction: float
    constraint_count: int


@dataclass(frozen=True)
class HandoffDecision:
    """Pass/fail decision for a predeclared representability handoff."""

    eligible: bool
    failed_reasons: list[str]
    warning_reasons: list[str]


@dataclass(frozen=True)
class ResponseConstraintHandoffManifest:
    """Frozen input specification for future representability experiments.

    The constraints are response-comparison constraints, not distance-order
    constraints. This manifest does not contain fitted embeddings.
    """

    manifest_id: str
    created_by_milestone: str
    profile_label: str
    comparison_protocol_name: str
    comparison_method: str
    missing_policy: str
    min_common_protocols: int
    min_margin: float
    max_constraints: int
    constraint_seed: int
    train_fraction: float
    validation_gate_name: str
    validation_summary: HandoffValidationSummary
    handoff_decision: HandoffDecision
    target_event_ids: NDArray[np.int_]
    constraints: NDArray[np.int_]
    margins: NDArray[np.float64]
    train_constraint_indices: NDArray[np.int_]
    heldout_constraint_indices: NDArray[np.int_]
    null_baseline_labels: list[str]
    forbidden_interpretations: list[str]


def forbidden_interpretations_default() -> list[str]:
    """Return forbidden interpretations for handoff manifests."""

    return [
        "distance-order constraints",
        "spatial distance",
        "metric distance",
        "metric reconstruction",
        "geometry recovery",
        "calibrated radar distance",
    ]


def build_handoff_validation_summary(
    heldout_row: dict[str, float],
    bootstrap_row: dict[str, float],
    null_row: dict[str, float],
    coverage_row: dict[str, float | str],
    constraint_count: int,
) -> HandoffValidationSummary:
    """Build a handoff validation summary from component diagnostics."""

    return HandoffValidationSummary(
        heldout_agreement_fraction=float(heldout_row["agreement_fraction"]),
        heldout_inversion_fraction=float(heldout_row["inversion_fraction"]),
        heldout_evaluable_fraction=float(heldout_row["evaluable_fraction"]),
        bootstrap_mean_agreement_fraction=float(
            bootstrap_row["mean_agreement_fraction"]
        ),
        bootstrap_stable_constraint_fraction=float(
            bootstrap_row["stable_constraint_fraction"]
        ),
        null_z_score=float(null_row["null_z_score"]),
        best_null_agreement_fraction=float(null_row["best_null_agreement_fraction"]),
        target_coverage_fraction=float(coverage_row["touched_target_fraction"]),
        pair_node_coverage_fraction=float(coverage_row["touched_pair_node_fraction"]),
        constraint_count=int(constraint_count),
    )


def decide_handoff_eligibility(
    summary: HandoffValidationSummary,
    gate: ConstraintValidationGate,
    *,
    min_target_coverage_fraction: float = 0.8,
    min_pair_node_coverage_fraction: float = 0.5,
) -> HandoffDecision:
    """Decide whether a response-comparison pool is eligible for handoff."""

    failed: list[str] = []
    if summary.constraint_count < gate.min_constraint_count:
        failed.append("too_few_constraints")
    if summary.heldout_evaluable_fraction < gate.min_evaluable_fraction:
        failed.append("low_heldout_evaluable_fraction")
    if summary.heldout_agreement_fraction < gate.min_agreement_fraction:
        failed.append("low_heldout_agreement")
    if summary.heldout_inversion_fraction > gate.max_inversion_fraction:
        failed.append("high_heldout_inversion")
    if (
        summary.bootstrap_mean_agreement_fraction < gate.min_bootstrap_confidence
        or summary.bootstrap_stable_constraint_fraction < gate.min_bootstrap_confidence
    ):
        failed.append("low_bootstrap_confidence")
    if summary.null_z_score < gate.min_null_z_score:
        failed.append("low_null_separation")
    if summary.target_coverage_fraction < min_target_coverage_fraction:
        failed.append("low_target_coverage")
    if summary.pair_node_coverage_fraction < min_pair_node_coverage_fraction:
        failed.append("low_pair_node_coverage")

    warnings: list[str] = []
    if summary.best_null_agreement_fraction >= summary.heldout_agreement_fraction:
        warnings.append("best_null_matches_or_exceeds_heldout")
    return HandoffDecision(
        eligible=not failed,
        failed_reasons=failed,
        warning_reasons=warnings,
    )


def split_constraint_indices(
    constraint_count: int,
    train_fraction: float,
    seed: int | None = None,
) -> tuple[NDArray[np.int_], NDArray[np.int_]]:
    """Split constraint indices for future representation experiments only."""

    if constraint_count < 0:
        raise ValueError("constraint_count must be nonnegative")
    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must be between 0 and 1")
    indices = np.arange(constraint_count, dtype=int)
    if constraint_count == 0:
        return indices, indices.copy()
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(indices)
    train_count = int(round(constraint_count * train_fraction))
    if constraint_count > 1:
        train_count = min(max(train_count, 1), constraint_count - 1)
    train = np.sort(permutation[:train_count])
    heldout = np.sort(permutation[train_count:])
    return train.astype(int), heldout.astype(int)


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def manifest_to_jsonable(
    manifest: ResponseConstraintHandoffManifest,
) -> dict[str, object]:
    """Convert a handoff manifest to canonical JSON-compatible data."""

    return _jsonable(asdict(manifest))


def manifest_digest(manifest_jsonable: dict[str, object]) -> str:
    """Return a SHA256 digest over canonical manifest JSON.

    The ``manifest_id`` field is excluded from the digest so the digest can be
    stored as the manifest identifier without a circular dependency.
    """

    payload = dict(manifest_jsonable)
    payload.pop("manifest_id", None)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_handoff_manifest(
    manifest: ResponseConstraintHandoffManifest,
    output_path: Path,
) -> Path:
    """Write a handoff manifest JSON file."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    jsonable = manifest_to_jsonable(manifest)
    output_path.write_text(
        json.dumps(jsonable, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return output_path


def read_handoff_manifest(path: Path) -> dict[str, object]:
    """Read a handoff manifest JSON file."""

    return json.loads(path.read_text(encoding="utf-8"))
