"""Stop-condition utilities for carry-forward robustness decisions."""

from __future__ import annotations

from dataclasses import asdict, dataclass

ALLOWED_MODES = {
    "stress_tests_allowed",
    "provisional_only_with_caveats",
    "report_only_controls",
    "stop",
}


@dataclass(frozen=True)
class StopConditionReport:
    """Summary of whether future stress tests are allowed."""

    carry_forward_count: int
    provisional_count: int
    blocked_count: int
    report_only_count: int
    failed_control_count: int
    stress_tests_allowed: bool
    reason: str
    allowed_mode: str

    def __post_init__(self) -> None:
        if self.allowed_mode not in ALLOWED_MODES:
            allowed = ", ".join(sorted(ALLOWED_MODES))
            raise ValueError(f"allowed_mode must be one of: {allowed}")


def evaluate_stress_test_stop_condition(
    decision_rows: list[dict[str, str | float]],
) -> StopConditionReport:
    """Evaluate the stress-test stop condition from decision rows."""

    counts = {
        "carry_forward": 0,
        "provisional": 0,
        "blocked": 0,
        "report_only": 0,
        "failed_control": 0,
    }
    for row in decision_rows:
        decision = str(row.get("decision", ""))
        if decision in counts:
            counts[decision] += 1
    if counts["carry_forward"] > 0:
        allowed = True
        mode = "stress_tests_allowed"
        reason = "at_least_one_carry_forward_family"
    elif counts["provisional"] > 0:
        allowed = True
        mode = "provisional_only_with_caveats"
        reason = "only_provisional_family_available"
    elif decision_rows and (
        counts["report_only"] + counts["blocked"] + counts["failed_control"]
        == len(decision_rows)
    ):
        allowed = False
        mode = "report_only_controls"
        reason = "no_carry_forward_or_provisional_family"
    else:
        allowed = False
        mode = "stop"
        reason = "no_decision_rows"
    return StopConditionReport(
        carry_forward_count=counts["carry_forward"],
        provisional_count=counts["provisional"],
        blocked_count=counts["blocked"],
        report_only_count=counts["report_only"],
        failed_control_count=counts["failed_control"],
        stress_tests_allowed=allowed,
        reason=reason,
        allowed_mode=mode,
    )


def stop_condition_to_row(report: StopConditionReport) -> dict[str, float | str]:
    """Convert stop condition report to CSV row."""

    row = asdict(report)
    row["stress_tests_allowed"] = float(report.stress_tests_allowed)
    return row
