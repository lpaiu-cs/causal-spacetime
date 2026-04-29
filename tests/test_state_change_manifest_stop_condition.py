from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_stop_condition import (
    evaluate_stress_test_stop_condition,
    stop_condition_to_row,
)


def test_evaluate_stress_test_stop_condition_stop_case() -> None:
    report = evaluate_stress_test_stop_condition(
        [
            {"family_name": "blocked", "decision": "blocked"},
            {"family_name": "reported", "decision": "report_only"},
        ]
    )

    assert not report.stress_tests_allowed
    assert report.allowed_mode == "report_only_controls"


def test_evaluate_stress_test_stop_condition_carry_forward_case() -> None:
    report = evaluate_stress_test_stop_condition(
        [{"family_name": "ok", "decision": "carry_forward"}]
    )

    assert report.stress_tests_allowed
    assert report.allowed_mode == "stress_tests_allowed"
    assert stop_condition_to_row(report)["stress_tests_allowed"] == 1.0
