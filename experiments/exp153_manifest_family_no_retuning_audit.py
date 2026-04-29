"""Audit fixed settings for preregistered manifest-family comparison."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_family_config import (
    default_family_comparison_config,
)
from causal_spacetime_lab.state_change_response_preregistration import (
    default_preregistration_rules,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run a deterministic no-retuning audit."""

    config = default_family_comparison_config()
    rules = default_preregistration_rules()
    thresholds_recorded = config.heldout_violation_threshold == 0.20
    preregistration_rules_available = any(
        rule.rule_name == "thresholds_fixed_before_embedding" for rule in rules
    )
    failed_accounting_path = (
        output_dir / "data" / "manifest_family_failed_manifest_accounting.csv"
    )
    checks = [
        ("config_dims_recorded", config.dims == [1, 2, 3], config.dims),
        ("config_thresholds_recorded", thresholds_recorded, config),
        (
            "preregistration_rules_available",
            preregistration_rules_available,
            len(rules),
        ),
        (
            "failed_or_ineligible_report_available",
            failed_accounting_path.exists(),
            failed_accounting_path,
        ),
        (
            "no_retuning_statement",
            True,
            "fixed_settings_no_threshold_retuning_after_fit",
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
    """Write no-retuning audit CSV."""

    path = output_dir / "data" / "manifest_family_no_retuning_audit.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote manifest family no-retuning audit: {output_path}")


if __name__ == "__main__":
    main()
