"""Shared helpers for cross-family robustness experiments."""

from __future__ import annotations

import json
from pathlib import Path

from manifest_family_experiment_helpers import (
    assignments_for_datasets,
    load_or_create_family_datasets,
)

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    CarryForwardRegistry,
    build_carry_forward_registry,
    registry_to_jsonable,
)
from causal_spacetime_lab.state_change_manifest_family import (
    ManifestFamilyAssignment,
)
from causal_spacetime_lab.state_change_manifest_family_outputs import (
    load_family_output_bundle,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
    FamilyRobustnessDecision,
    decide_family_robustness,
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_family_robustness_metrics import (
    aggregate_family_robustness_metrics,
)


def decision_to_row(
    decision: FamilyRobustnessDecision,
    *,
    missing_inputs: list[str] | None = None,
) -> dict[str, float | str]:
    """Convert a robustness decision to a flat CSV row."""

    row: dict[str, float | str] = {
        "family_name": decision.family_name,
        "family_kind": decision.family_kind,
        "decision": decision.decision,
        "passed": float(decision.passed),
        "failed_reasons": ";".join(decision.failed_reasons),
        "warning_reasons": ";".join(decision.warning_reasons),
        "missing_inputs": ";".join(missing_inputs or []),
    }
    row.update(decision.key_metrics)
    return row


def decision_from_row(row: dict[str, str]) -> FamilyRobustnessDecision:
    """Reconstruct a decision object from a CSV row."""

    metric_keys = [
        "manifest_count",
        "fitted_fraction",
        "no_fit_fraction",
        "mean_heldout_violation",
        "mean_generalization_gap",
        "stricter_threshold_pass_fraction",
        "destructive_null_gap",
        "symmetry_control_gap",
        "target_coverage_fraction",
        "pair_node_coverage_fraction",
        "restart_std",
        "latent_order_disagreement",
        "no_retuning_audit_pass",
        "failed_accounting_present",
    ]
    return FamilyRobustnessDecision(
        family_name=row.get("family_name", ""),
        family_kind=row.get("family_kind", ""),
        decision=row.get("decision", "blocked"),
        passed=float(row.get("passed", 0.0)) == 1.0,
        failed_reasons=[
            reason for reason in row.get("failed_reasons", "").split(";") if reason
        ],
        warning_reasons=[
            reason for reason in row.get("warning_reasons", "").split(";") if reason
        ],
        key_metrics={key: float(row.get(key, float("nan"))) for key in metric_keys},
    )


def load_metrics_and_decisions(
    output_dir: Path,
    criteria: CrossFamilyRobustnessCriteria | None = None,
) -> tuple[
    list[dict[str, float | str]],
    list[FamilyRobustnessDecision],
    list[str],
]:
    """Load family outputs and compute robustness metrics and decisions."""

    active_criteria = criteria or default_cross_family_robustness_criteria()
    bundle = load_family_output_bundle(output_dir)
    missing_inputs = [
        key for key, rows in bundle.items() if not rows
    ]
    metrics = aggregate_family_robustness_metrics(bundle)
    decisions = [
        decide_family_robustness(row, active_criteria)
        for row in metrics
    ]
    return metrics, decisions, missing_inputs


def load_default_assignments(output_dir: Path) -> list[ManifestFamilyAssignment]:
    """Load manifests and assign them to default families."""

    datasets = load_or_create_family_datasets(output_dir / "manifests", output_dir)
    return assignments_for_datasets(datasets)


def build_registry_from_output(
    output_dir: Path,
    criteria_name: str = "default",
) -> CarryForwardRegistry:
    """Build a carry-forward registry from current output CSVs."""

    _, decisions, _ = load_metrics_and_decisions(output_dir)
    assignments = load_default_assignments(output_dir)
    return build_carry_forward_registry(decisions, assignments, criteria_name)


def read_registry_json(path: Path) -> dict[str, object]:
    """Read a carry-forward registry JSON file."""

    return json.loads(path.read_text(encoding="utf-8"))


def registry_summary_rows(
    registry: CarryForwardRegistry,
) -> list[dict[str, float | str]]:
    """Return flat summary rows for a carry-forward registry."""

    payload = registry_to_jsonable(registry)
    rows: list[dict[str, float | str]] = []
    for record in payload["records"]:  # type: ignore[index]
        item = dict(record)  # type: ignore[arg-type]
        rows.append(
            {
                "registry_id": str(payload["registry_id"]),
                "family_name": str(item["family_name"]),
                "family_kind": str(item["family_kind"]),
                "decision": str(item["decision"]),
                "eligible_manifest_count": float(item["eligible_manifest_count"]),
                "ineligible_manifest_count": float(item["ineligible_manifest_count"]),
                "allowed_future_use": str(item["allowed_future_use"]),
                "failed_reasons": ";".join(item["failed_reasons"]),
                "warning_reasons": ";".join(item["warning_reasons"]),
            }
        )
    return rows
