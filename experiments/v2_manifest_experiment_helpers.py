"""Shared helpers for Milestone 37 v2 manifest-generation experiments."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    build_future_manifest_run_spec_from_plan,
    future_run_spec_to_jsonable,
)
from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
    default_new_manifest_family_specs_v2,
    write_remediation_plan,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import (
    V2ManifestFamilySpec,
    load_v2_family_specs_from_design_csv,
    load_v2_family_specs_from_remediation_plan,
    v2_family_spec_table,
)

DEFAULT_OUTPUT_DIR = Path("outputs")
REQUIRED_V2_METRICS = [
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


def parse_int_values(values: list[str]) -> tuple[int, ...]:
    """Parse CLI int values, accepting repeated or comma-separated forms."""

    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def data_path(output_dir: Path, filename: str) -> Path:
    """Return a path under the output data directory."""

    return output_dir / "data" / filename


def figure_path(output_dir: Path, filename: str) -> Path:
    """Return a path under the output figure directory."""

    return output_dir / "figures" / filename


def remediation_path(output_dir: Path, filename: str) -> Path:
    """Return a path under the output remediation directory."""

    return output_dir / "remediation" / filename


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows, returning an empty list when the file is absent."""

    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def read_json_object(path: Path) -> dict[str, object]:
    """Read a JSON object, returning an empty object when absent."""

    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload if isinstance(payload, dict) else {}


def write_json_object(path: Path, payload: dict[str, object]) -> Path:
    """Write a JSON object."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def missing_input_rows(paths: list[Path]) -> list[dict[str, float | str]]:
    """Return explicit missing-input rows for absent paths."""

    return [
        {"input_path": str(path), "present": float(path.exists())}
        for path in paths
        if not path.exists()
    ]


def load_v2_specs_with_fallback(output_dir: Path) -> list[V2ManifestFamilySpec]:
    """Load v2 specs from M36 outputs or fall back to default planned specs."""

    design_path = data_path(output_dir, "new_manifest_family_design_v2.csv")
    plan_path = remediation_path(output_dir, "remediation_plan_m36.json")
    if design_path.exists():
        return load_v2_family_specs_from_design_csv(design_path)
    if plan_path.exists():
        return load_v2_family_specs_from_remediation_plan(plan_path)
    return [
        V2ManifestFamilySpec(
            family_name=str(row["family_name"]),
            family_kind=str(row["family_kind"]),
            profile_column_policy=str(row["profile_column_policy"]),
            target_inclusion_policy=str(row["target_inclusion_policy"]),
            pair_node_coverage_policy=str(row["pair_node_coverage_policy"]),
            null_taxonomy_policy=str(row["null_taxonomy_policy"]),
            restart_stability_required=bool(row["restart_stability_required"]),
            latent_order_stability_required=bool(
                row["latent_order_stability_required"]
            ),
            execution_status=str(row["execution_status"]),
        )
        for row in default_new_manifest_family_specs_v2()
    ]


def ensure_m36_plan_and_spec_for_tests(output_dir: Path) -> tuple[Path, Path]:
    """Create minimal M36 plan/spec artifacts when a test output dir is empty."""

    plan_path = remediation_path(output_dir, "remediation_plan_m36.json")
    spec_path = remediation_path(output_dir, "future_manifest_run_spec_m36.json")
    if plan_path.exists() and spec_path.exists():
        return plan_path, spec_path
    plan = build_preregistered_remediation_plan([], [])
    write_remediation_plan(plan, plan_path)
    spec = build_future_manifest_run_spec_from_plan(plan)
    write_json_object(spec_path, future_run_spec_to_jsonable(spec))
    write_csv(
        default_new_manifest_family_specs_v2(),
        data_path(output_dir, "new_manifest_family_design_v2.csv"),
        [
            "family_name",
            "family_kind",
            "profile_column_policy",
            "target_inclusion_policy",
            "pair_node_coverage_policy",
            "null_taxonomy_policy",
            "restart_stability_required",
            "latent_order_stability_required",
            "execution_status",
        ],
    )
    return plan_path, spec_path


def save_metric_bar(
    rows: list[dict[str, float | str]],
    *,
    label_key: str,
    value_key: str,
    path: Path,
    ylabel: str,
) -> Path:
    """Save a simple metric bar figure."""

    import matplotlib.pyplot as plt

    labels: list[str] = []
    values: list[float] = []
    for row in rows:
        try:
            value = float(row[value_key])
        except (KeyError, TypeError, ValueError):
            continue
        if not np.isfinite(value):
            continue
        labels.append(str(row.get(label_key, "")))
        values.append(value)

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


def v2_spec_rows(output_dir: Path) -> list[dict[str, float | str]]:
    """Return current v2 family specs as table rows."""

    return v2_family_spec_table(load_v2_specs_with_fallback(output_dir))
