"""Final exact sanity checks for Milestone 38 v2 carry-forward evaluation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

DEFAULT_OUTPUT_DIR = Path("outputs")


def parse_args() -> Path:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 carry-forward final sanity.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args().output_dir


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run final sanity checks."""

    registry_path = output_dir / "carry_forward_v2" / "carry_forward_registry_v2.json"
    registry = {}
    if registry_path.exists():
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    plan_rows = read_csv_rows(data_path(output_dir, "v2_stress_test_handoff_plan.csv"))
    blocked_rows_blocked = all(
        float(row.get("allowed", 0.0)) == 0.0
        for row in plan_rows
        if row.get("decision") in {"blocked", "report_only", "failed_control"}
    )
    decision_rows = read_csv_rows(
        data_path(output_dir, "v2_cross_family_robustness_decisions.csv")
    )
    checks = [
        (
            "v2_decisions_file_exists",
            data_path(output_dir, "v2_cross_family_robustness_decisions.csv").exists(),
        ),
        ("v2_registry_json_serializes", bool(registry.get("registry_id"))),
        ("blocked_families_not_allowed_in_plan", blocked_rows_blocked),
        (
            "no_stress_test_result_output",
            not any((output_dir / "data").glob("v2_*stress_test_result*.csv")),
        ),
        (
            "diagnostic_complete_flag_exists",
            bool(decision_rows) and "diagnostic_complete" in decision_rows[0],
        ),
        (
            "registry_forbidden_interpretations_present",
            bool(registry.get("forbidden_interpretations")),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write final sanity rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_carry_forward_final_sanity.csv"),
        ["check", "passed"],
    )


def main() -> None:
    output_dir = parse_args()
    path = write_outputs(run_experiment(output_dir), output_dir)
    print(f"Wrote v2 carry-forward final sanity: {path}")


if __name__ == "__main__":
    main()
