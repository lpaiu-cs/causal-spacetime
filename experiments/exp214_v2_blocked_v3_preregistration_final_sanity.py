"""Final exact sanity checks for M39 v2 blocked audit and v3 preregistration."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    V2BlockingCriterionRecord,
)
from causal_spacetime_lab.state_change_manifest_v2_criterion_margins import (
    criterion_margins_from_blocking_records,
)
from causal_spacetime_lab.state_change_manifest_v3_audit import (
    audit_v3_preregistration_only,
    check_no_v3_execution_outputs,
)
from causal_spacetime_lab.state_change_manifest_v3_design import (
    default_v3_manifest_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v3_preregistration import (
    build_v3_preregistration_spec,
)


def run_experiment(output_dir: Path = Path("outputs")) -> list[dict[str, float | str]]:
    pass_record = V2BlockingCriterionRecord(
        family_name="f",
        family_kind="structured",
        criterion_name="manifest_count",
        observed_value=5.0,
        threshold_value=3.0,
        criterion_direction="min_required",
        blocking_type="not_blocking",
        status="pass",
        explanation="pass",
    )
    fail_record = V2BlockingCriterionRecord(
        family_name="f",
        family_kind="structured",
        criterion_name="mean_heldout_violation",
        observed_value=0.5,
        threshold_value=0.2,
        criterion_direction="max_allowed",
        blocking_type="measured_blocking",
        status="fail",
        explanation="fail",
    )
    margins = criterion_margins_from_blocking_records([pass_record, fail_record])
    designs = default_v3_manifest_family_designs()
    spec = build_v3_preregistration_spec(designs)
    audit = audit_v3_preregistration_only(spec)
    execution = check_no_v3_execution_outputs(output_dir)
    report_rows = read_csv_rows(
        data_path(output_dir, "v2_blocked_decision_report_card.csv")
    )
    blocked_reported = any(row.get("decision") == "blocked" for row in report_rows)
    return [
        {"check": "blocking_record_allowed_type", "passed": 1.0},
        {
            "check": "criterion_margin_sign_convention",
            "passed": float(margins[0].margin > 0 and margins[1].margin < 0),
        },
        {
            "check": "structured_v3_planned_count_at_least_3",
            "passed": float(
                all(
                    design.planned_manifest_count >= 3
                    for design in designs
                    if design.family_kind == "structured"
                )
            ),
        },
        {
            "check": "v3_preregistration_execution_disallowed",
            "passed": float(not spec.execution_allowed_in_current_milestone),
        },
        {
            "check": "v3_no_execution_outputs_absent",
            "passed": float(all(float(value) == 1.0 for value in execution.values())),
        },
        {
            "check": "v3_preregistration_audit_passes",
            "passed": float(all(float(value) == 1.0 for value in audit.values())),
        },
        {
            "check": "report_card_includes_blocked_v2_families",
            "passed": float(blocked_reported),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v2_blocked_v3_preregistration_final_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v2 blocked/v3 preregistration final sanity: {path}")


if __name__ == "__main__":
    main()
