from __future__ import annotations

from causal_spacetime_lab.state_change_response_preregistration import (
    default_preregistration_rules,
    preregistration_rule_table,
)


def test_default_preregistration_rules_contains_no_metric_overclaim_rule() -> None:
    rules = default_preregistration_rules()

    assert any(
        rule.rule_name == "no_metric_interpretation_without_calibration"
        for rule in rules
    )


def test_preregistration_rule_table() -> None:
    rows = preregistration_rule_table()

    assert rows
    assert {"rule_name", "description", "violation_example"} <= set(rows[0])
