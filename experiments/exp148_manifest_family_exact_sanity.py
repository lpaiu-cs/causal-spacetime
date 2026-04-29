"""Exact sanity checks for manifest-family assignment and null taxonomy."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from manifest_family_experiment_helpers import ensure_failed_control_manifest
from manifest_representation_experiment_helpers import build_exact_manifest

from causal_spacetime_lab.state_change_manifest_dataset import load_manifest_dataset
from causal_spacetime_lab.state_change_manifest_family import (
    assign_manifest_to_family,
    default_manifest_family_specs,
)
from causal_spacetime_lab.state_change_manifest_null_taxonomy import (
    classify_null_type,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run exact family-assignment checks."""

    eligible_path = build_exact_manifest(output_dir, "manifest_family_gap.json")
    eligible = load_manifest_dataset(eligible_path)
    ineligible_payload = json.loads(eligible_path.read_text(encoding="utf-8"))
    ineligible_payload["manifest_id"] = "manifest_family_ineligible"
    ineligible_payload["handoff_decision"]["eligible"] = False
    ineligible_payload["handoff_decision"]["failed_reasons"] = ["low_test_signal"]
    ineligible_path = output_dir / "manifests" / "manifest_family_ineligible.json"
    ineligible_path.write_text(
        json.dumps(ineligible_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    ineligible = load_manifest_dataset(ineligible_path)
    failed = ensure_failed_control_manifest(output_dir)

    specs = default_manifest_family_specs()
    eligible_assignment = assign_manifest_to_family(eligible, specs)
    ineligible_assignment = assign_manifest_to_family(ineligible, specs)
    failed_assignment = assign_manifest_to_family(failed, specs)
    checks = [
        (
            "eligible_gap_family",
            eligible_assignment.family_name == "eligible_rank_gap",
            eligible_assignment.family_name,
        ),
        (
            "ineligible_family",
            ineligible_assignment.family_name == "ineligible_reported",
            ineligible_assignment.family_name,
        ),
        (
            "failed_control_family",
            failed_assignment.family_name == "failed_synthetic_controls",
            failed_assignment.family_name,
        ),
        (
            "permuted_targets_symmetry_control",
            classify_null_type("permuted_targets") == "symmetry_control",
            classify_null_type("permuted_targets"),
        ),
        (
            "shuffled_sides_destructive_null",
            classify_null_type("shuffled_sides") == "destructive_null",
            classify_null_type("shuffled_sides"),
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
    """Write exact family sanity CSV."""

    path = output_dir / "data" / "manifest_family_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote manifest family exact sanity: {output_path}")


if __name__ == "__main__":
    main()
