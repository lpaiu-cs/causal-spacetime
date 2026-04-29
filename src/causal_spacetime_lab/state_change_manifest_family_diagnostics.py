"""Family-level diagnostics for frozen-manifest representation fits."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_family import (
    ManifestFamilyAssignment,
)
from causal_spacetime_lab.state_change_manifest_representation import (
    ManifestRepresentationFit,
)


@dataclass(frozen=True)
class ManifestFitDiagnosticRow:
    """One manifest fit row annotated with family assignment."""

    manifest_id: str
    family_name: str
    family_kind: str
    eligible: bool
    fitted: bool
    reason_not_fit: str
    embedding_dim: int
    train_violation_rate: float
    heldout_violation_rate: float
    generalization_gap: float
    train_hinge_loss: float
    heldout_hinge_loss: float
    target_count: int
    train_constraint_count: int
    heldout_constraint_count: int


def manifest_fit_diagnostic_row(
    fit: ManifestRepresentationFit,
    assignment: ManifestFamilyAssignment,
) -> ManifestFitDiagnosticRow:
    """Attach family metadata to a latent representation fit."""

    gap = float(fit.heldout_violation_rate) - float(fit.train_violation_rate)
    return ManifestFitDiagnosticRow(
        manifest_id=fit.manifest_id,
        family_name=assignment.family_name,
        family_kind=assignment.family_kind,
        eligible=assignment.eligible,
        fitted=fit.fitted,
        reason_not_fit=fit.reason_not_fit,
        embedding_dim=int(fit.embedding_dim),
        train_violation_rate=float(fit.train_violation_rate),
        heldout_violation_rate=float(fit.heldout_violation_rate),
        generalization_gap=gap,
        train_hinge_loss=float(fit.train_hinge_loss),
        heldout_hinge_loss=float(fit.heldout_hinge_loss),
        target_count=int(fit.target_count),
        train_constraint_count=int(fit.train_constraint_count),
        heldout_constraint_count=int(fit.heldout_constraint_count),
    )


def fit_diagnostic_row_to_dict(
    row: ManifestFitDiagnosticRow,
) -> dict[str, float | str]:
    """Convert a family fit diagnostic row to CSV-safe values."""

    result = asdict(row)
    result["eligible"] = float(row.eligible)
    result["fitted"] = float(row.fitted)
    return result


def _finite(values: list[float]) -> np.ndarray:
    return np.asarray([value for value in values if np.isfinite(value)], dtype=float)


def summarize_family_fit_diagnostics(
    rows: list[ManifestFitDiagnosticRow],
) -> list[dict[str, float | str]]:
    """Summarize fit diagnostics by family and latent dimension."""

    grouped: dict[tuple[str, str, int], list[ManifestFitDiagnosticRow]] = defaultdict(
        list
    )
    for row in rows:
        grouped[(row.family_name, row.family_kind, row.embedding_dim)].append(row)

    summaries: list[dict[str, float | str]] = []
    for (family_name, family_kind, dim), group in sorted(grouped.items()):
        train = _finite([row.train_violation_rate for row in group if row.fitted])
        heldout = _finite([row.heldout_violation_rate for row in group if row.fitted])
        gap = _finite([row.generalization_gap for row in group if row.fitted])
        summaries.append(
            {
                "family_name": family_name,
                "family_kind": family_kind,
                "embedding_dim": float(dim),
                "manifest_count": float(len({row.manifest_id for row in group})),
                "fitted_count": float(sum(row.fitted for row in group)),
                "no_fit_count": float(sum(not row.fitted for row in group)),
                "mean_train_violation": (
                    float(np.mean(train)) if train.size else float("nan")
                ),
                "mean_heldout_violation": (
                    float(np.mean(heldout)) if heldout.size else float("nan")
                ),
                "mean_generalization_gap": (
                    float(np.mean(gap)) if gap.size else float("nan")
                ),
                "median_heldout_violation": (
                    float(np.median(heldout)) if heldout.size else float("nan")
                ),
                "best_heldout_violation": (
                    float(np.min(heldout)) if heldout.size else float("nan")
                ),
                "worst_heldout_violation": (
                    float(np.max(heldout)) if heldout.size else float("nan")
                ),
            }
        )
    return summaries


def failed_manifest_accounting_summary(
    assignments: list[ManifestFamilyAssignment],
) -> list[dict[str, float | str]]:
    """Report family counts and failure-reason counts."""

    rows: list[dict[str, float | str]] = []
    family_counts = Counter(assignment.family_name for assignment in assignments)
    for family_name, count in sorted(family_counts.items()):
        rows.append(
            {
                "row_type": "family_count",
                "family_name": family_name,
                "family_kind": next(
                    assignment.family_kind
                    for assignment in assignments
                    if assignment.family_name == family_name
                ),
                "reason": "",
                "count": float(count),
            }
        )
    reason_counts: Counter[str] = Counter()
    for assignment in assignments:
        if assignment.eligible and not assignment.failed_reasons:
            continue
        for reason in assignment.failed_reasons or ["ineligible_without_reason"]:
            reason_counts[reason] += 1
    for reason, count in sorted(reason_counts.items()):
        rows.append(
            {
                "row_type": "failure_reason",
                "family_name": "",
                "family_kind": "",
                "reason": reason,
                "count": float(count),
            }
        )
    return rows
