"""Audit that failure decomposition did not retune thresholds or run stress tests."""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from carry_forward_failure_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_remediation_design import (
    default_remediation_proposals,
)
from causal_spacetime_lab.state_change_manifest_stop_condition import (
    evaluate_stress_test_stop_condition,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _criteria_table_matches_default(output_dir: Path) -> bool:
    table = read_csv_rows(
        data_path(output_dir, "cross_family_robustness_criteria_table.csv")
    )
    default = asdict(default_cross_family_robustness_criteria())
    exported = {row.get("criterion", ""): row.get("value", "") for row in table}
    if not exported:
        return False
    for key, value in default.items():
        if key == "name":
            continue
        exported_value = exported.get(key, "")
        if isinstance(value, bool):
            if str(exported_value).lower() != str(value).lower():
                return False
            continue
        try:
            if float(exported_value) != float(value):
                return False
        except (TypeError, ValueError):
            return False
    return True


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run no-retuning compliance checks."""

    decisions = read_csv_rows(
        data_path(output_dir, "cross_family_robustness_decisions.csv")
    )
    stop_report = evaluate_stress_test_stop_condition(decisions)
    future_stress_outputs = list(
        (output_dir / "data").glob("representability_stress_test_*.csv")
    ) + list((output_dir / "data").glob("future_stress_test_*.csv"))
    proposals = default_remediation_proposals()
    execution_change_allowed = any(
        proposal.requires_new_preregistration
        and proposal.allowed_in_current_milestone
        for proposal in proposals
    )
    checks = [
        (
            "criteria_match_exported_table",
            _criteria_table_matches_default(output_dir),
        ),
        ("no_future_stress_test_outputs", not future_stress_outputs),
        ("execution_changing_remediations_not_allowed", not execution_change_allowed),
        (
            "stop_condition_respected",
            (not stop_report.stress_tests_allowed)
            if stop_report.carry_forward_count == 0
            and stop_report.provisional_count == 0
            else True,
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write no-retuning audit CSV."""

    path = output_dir / "data" / "failure_decomposition_no_retuning_audit.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote failure-decomposition no-retuning audit: {path}")


if __name__ == "__main__":
    main()
