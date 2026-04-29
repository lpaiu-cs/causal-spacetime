"""Final exact sanity checks for carry-forward failure decomposition."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_diagnostic_completeness import (
    diagnostic_completeness_for_family,
)
from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
)
from causal_spacetime_lab.state_change_manifest_remediation_design import (
    default_remediation_proposals,
)
from causal_spacetime_lab.state_change_manifest_stop_condition import (
    evaluate_stress_test_stop_condition,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run final deterministic checks."""

    validation_raised = False
    try:
        CriterionFailureRecord(
            family_name="x",
            family_kind="structured",
            criterion_name="bad",
            criterion_type="invalid",
            observed_value=0.0,
            threshold_value=1.0,
            status="pass",
            root_cause_category="unknown",
            explanation="invalid",
        )
    except ValueError:
        validation_raised = True
    completeness = diagnostic_completeness_for_family(
        {"family_name": "missing", "manifest_count": 1.0}
    )
    stop_report = evaluate_stress_test_stop_condition(
        [
            {"family_name": "blocked", "decision": "blocked"},
            {"family_name": "reported", "decision": "report_only"},
        ]
    )
    proposals = default_remediation_proposals()
    data_changing_allowed = any(
        proposal.requires_new_preregistration
        and proposal.allowed_in_current_milestone
        for proposal in proposals
    )
    report_path = output_dir / "data" / "carry_forward_failure_report_card.csv"
    report_text = (
        report_path.read_text(encoding="utf-8") if report_path.exists() else ""
    )
    checks = [
        ("criterion_record_validation", validation_raised),
        (
            "diagnostic_completeness_detects_missing",
            completeness.missing_metric_count > 0,
        ),
        (
            "stop_without_carry_forward_or_provisional",
            not stop_report.stress_tests_allowed,
        ),
        ("new_data_remediations_not_allowed_now", not data_changing_allowed),
        (
            "report_card_mentions_noncarry_families",
            "blocked" in report_text
            or "report_only" in report_text
            or "failed_control" in report_text,
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write final sanity CSV."""

    path = output_dir / "data" / "carry_forward_failure_final_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote carry-forward failure final sanity: {path}")


if __name__ == "__main__":
    main()
