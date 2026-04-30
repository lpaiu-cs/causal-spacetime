"""Final exact sanity checks for Milestone 37 v2 manifest generation."""

from __future__ import annotations

import argparse
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    REQUIRED_V2_METRICS,
    data_path,
    load_v2_specs_with_fallback,
    read_csv_rows,
    read_json_object,
    remediation_path,
)


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run final deterministic M37 sanity checks."""

    specs = load_v2_specs_with_fallback(output_dir)
    manifest_dir = output_dir / "manifests_v2"
    metric_rows = read_csv_rows(
        data_path(output_dir, "v2_cross_family_robustness_metrics.csv")
    )
    plan = read_json_object(remediation_path(output_dir, "remediation_plan_m36.json"))
    planned_specs = plan.get("new_manifest_family_specs", [])
    plan_still_planned = True
    if isinstance(planned_specs, list):
        plan_still_planned = all(
            isinstance(row, dict) and row.get("execution_status") == "planned_only"
            for row in planned_specs
        )
    checks = [
        ("v2_specs_exist", bool(specs)),
        (
            "v2_manifests_written_to_manifests_v2",
            manifest_dir.exists() and bool(list(manifest_dir.glob("*.json"))),
        ),
        (
            "v2_required_metrics_include_all_14",
            bool(metric_rows)
            and all(metric in metric_rows[0] for metric in REQUIRED_V2_METRICS),
        ),
        (
            "no_v2_carry_forward_decision_file",
            not data_path(
                output_dir,
                "v2_cross_family_robustness_decisions.csv",
            ).exists(),
        ),
        (
            "no_v2_stress_test_outputs",
            not any((output_dir / "data").glob("v2_*stress*.csv")),
        ),
        ("m36_plan_specs_remain_planned_only", plan_still_planned),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write final M37 sanity rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_manifest_generation_final_sanity.csv"),
        ["check", "passed"],
    )


def parse_args() -> Path:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 final generation sanity.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args().output_dir


def main() -> None:
    output_dir = parse_args()
    path = write_outputs(run_experiment(output_dir), output_dir)
    print(f"Wrote v2 manifest generation final sanity: {path}")


if __name__ == "__main__":
    main()
