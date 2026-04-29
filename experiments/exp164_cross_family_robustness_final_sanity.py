"""Final exact sanity checks for cross-family robustness handoff."""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    build_carry_forward_registry,
    forbidden_carry_forward_interpretations,
    registry_digest,
    registry_to_jsonable,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    FamilyRobustnessDecision,
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_stress_handoff import (
    build_stress_test_plan,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, float | str]]:
    """Run final deterministic checks."""

    criteria_before = asdict(default_cross_family_robustness_criteria())
    blocked = FamilyRobustnessDecision(
        family_name="blocked_family",
        family_kind="structured",
        decision="blocked",
        passed=False,
        failed_reasons=["high_heldout_violation"],
        warning_reasons=[],
        key_metrics={"mean_heldout_violation": 0.4},
    )
    provisional = FamilyRobustnessDecision(
        family_name="provisional_family",
        family_kind="structured",
        decision="provisional",
        passed=False,
        failed_reasons=["missing_target_coverage"],
        warning_reasons=["missing_target_coverage"],
        key_metrics={"mean_heldout_violation": 0.1},
    )
    registry = build_carry_forward_registry([blocked, provisional])
    payload = registry_to_jsonable(registry)
    plan = build_stress_test_plan(registry)
    criteria_after = asdict(default_cross_family_robustness_criteria())
    report_path = output_dir / "data" / "cross_family_robustness_report_card.csv"
    report_text = (
        report_path.read_text(encoding="utf-8") if report_path.exists() else ""
    )
    has_forbidden_interpretations = bool(forbidden_carry_forward_interpretations())
    checks = [
        ("default_criteria_exist", bool(criteria_before)),
        ("forbidden_interpretations_nonempty", has_forbidden_interpretations),
        ("registry_json_serializes", bool(registry_digest(payload))),
        ("stress_plan_blocks_blocked", not plan[0].allowed),
        ("threshold_sensitivity_preserves_defaults", criteria_before == criteria_after),
        (
            "report_card_includes_noncarry_families",
            "blocked" in report_text or "provisional" in report_text,
        ),
    ]
    return [
        {"check": name, "passed": float(passed)}
        for name, passed in checks
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write final exact sanity CSV."""

    path = output_dir / "data" / "cross_family_robustness_final_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote cross-family robustness final sanity: {path}")


if __name__ == "__main__":
    main()
