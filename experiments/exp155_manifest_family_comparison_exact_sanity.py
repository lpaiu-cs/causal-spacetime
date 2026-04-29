"""Exact sanity checks for manifest-family comparison utilities."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_dataset
from causal_spacetime_lab.state_change_manifest_family import (
    ManifestFamilySpec,
    assign_manifest_to_family,
    default_manifest_family_specs,
)
from causal_spacetime_lab.state_change_manifest_family_comparison import (
    compare_manifest_family_against_thresholds,
)
from causal_spacetime_lab.state_change_manifest_family_config import (
    default_family_comparison_config,
)
from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    ManifestFitDiagnosticRow,
    failed_manifest_accounting_summary,
    summarize_family_fit_diagnostics,
)
from causal_spacetime_lab.state_change_manifest_null_taxonomy import (
    default_null_taxonomy,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _ineligible_dataset(output_dir: Path) -> object:
    path = build_exact_manifest(output_dir, "manifest_family_exact_base.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["manifest_id"] = "manifest_family_exact_ineligible"
    payload["handoff_decision"]["eligible"] = False
    payload["handoff_decision"]["failed_reasons"] = ["exact_failure"]
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return load_manifest_dataset(path)


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run exact family comparison sanity checks."""

    eligible = load_manifest_dataset(
        build_exact_manifest(output_dir, "manifest_family_exact_gap.json")
    )
    ineligible = _ineligible_dataset(output_dir)
    specs = default_manifest_family_specs()
    eligible_assignment = assign_manifest_to_family(eligible, specs)
    ineligible_assignment = assign_manifest_to_family(ineligible, specs)
    has_permuted_targets = any(
        entry.null_type == "permuted_targets" for entry in default_null_taxonomy()
    )
    fit_rows = [
        ManifestFitDiagnosticRow(
            manifest_id="m1",
            family_name="eligible_rank_gap",
            family_kind="structured",
            eligible=True,
            fitted=True,
            reason_not_fit="",
            embedding_dim=1,
            train_violation_rate=0.05,
            heldout_violation_rate=0.10,
            generalization_gap=0.05,
            train_hinge_loss=0.01,
            heldout_hinge_loss=0.02,
            target_count=6,
            train_constraint_count=10,
            heldout_constraint_count=5,
        )
    ]
    summaries = summarize_family_fit_diagnostics(fit_rows)
    threshold_rows = compare_manifest_family_against_thresholds(
        fit_rows,
        default_family_comparison_config(),
    )
    failure_rows = failed_manifest_accounting_summary(
        [eligible_assignment, ineligible_assignment]
    )
    checks = [
        (
            "family_spec_validation",
            ManifestFamilySpec("x", "structured", "ok").family_kind == "structured",
            "ok",
        ),
        (
            "eligible_assignment",
            eligible_assignment.family_name == "eligible_rank_gap",
            eligible_assignment.family_name,
        ),
        (
            "ineligible_assignment",
            ineligible_assignment.family_name == "ineligible_reported",
            ineligible_assignment.family_name,
        ),
        ("family_summary_nonempty", bool(summaries), len(summaries)),
        (
            "null_taxonomy_contains_permuted_targets",
            has_permuted_targets,
            "permuted_targets",
        ),
        (
            "default_config_thresholds",
            default_family_comparison_config().heldout_violation_threshold == 0.20,
            default_family_comparison_config(),
        ),
        (
            "threshold_pass",
            threshold_rows[0]["threshold_pass"] == 1.0,
            threshold_rows[0]["threshold_pass"],
        ),
        (
            "failed_accounting_counts_reasons",
            any(row["row_type"] == "failure_reason" for row in failure_rows),
            failure_rows,
        ),
    ]
    return [
        {"check": name, "passed": float(passed), "value": str(value)}
        for name, passed, value in checks
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact family comparison sanity CSV."""

    path = output_dir / "data" / "manifest_family_comparison_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote manifest family comparison exact sanity: {output_path}")


if __name__ == "__main__":
    main()
