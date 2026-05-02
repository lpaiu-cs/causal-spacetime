"""Shared helpers for Milestone 45 v4 protocol carry-forward experiments."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    CarryForwardRegistry,
    registry_from_jsonable,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
)

KEY_METRIC_COLUMNS = [
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


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows, returning empty rows when absent."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def data_path(output_dir: Path, filename: str) -> Path:
    """Return an M45 output data path."""

    return output_dir / "data" / filename


def decisions_from_csv(path: Path) -> list[FamilyRobustnessDecision]:
    """Reconstruct v4 protocol family decisions from CSV rows."""

    decisions: list[FamilyRobustnessDecision] = []
    for row in read_csv_rows(path):
        key_metrics = {
            key: float(row[key])
            for key in KEY_METRIC_COLUMNS
            if key in row and row[key] != ""
        }
        decisions.append(
            FamilyRobustnessDecision(
                family_name=row.get("family_name", ""),
                family_kind=row.get("family_kind", ""),
                decision=row.get("decision", ""),
                passed=float(row.get("passed", 0.0)) == 1.0,
                failed_reasons=[
                    item for item in row.get("failed_reasons", "").split(";") if item
                ],
                warning_reasons=[
                    item for item in row.get("warning_reasons", "").split(";") if item
                ],
                key_metrics=key_metrics,
            )
        )
    return decisions


def read_v4_protocol_registry(output_dir: Path) -> CarryForwardRegistry | None:
    """Read the v4 protocol carry-forward registry if present."""

    path = (
        output_dir
        / "carry_forward_v4_protocol"
        / "carry_forward_registry_v4_protocol.json"
    )
    if not path.exists():
        return None
    return registry_from_jsonable(json.loads(path.read_text(encoding="utf-8")))
