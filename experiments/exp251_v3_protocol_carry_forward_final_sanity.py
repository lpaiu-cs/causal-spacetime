"""Final sanity checks for M42 v3 protocol carry-forward evaluation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_carry_forward_experiment_helpers import data_path, read_csv_rows


def parse_args() -> Path:
    parser = argparse.ArgumentParser(description="V3 protocol final sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args().output_dir


def run_experiment(output_dir: Path = Path("outputs")) -> list[dict[str, float | str]]:
    registry_path = (
        output_dir
        / "carry_forward_v3_protocol"
        / "carry_forward_registry_v3_protocol.json"
    )
    registry = {}
    if registry_path.exists():
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    plan_rows = read_csv_rows(
        data_path(output_dir, "v3_protocol_stress_test_handoff_plan.csv")
    )
    blocked_rows_blocked = all(
        float(row.get("allowed", 0.0)) == 0.0
        for row in plan_rows
        if row.get("decision") in {"blocked", "report_only", "failed_control"}
    )
    decisions = read_csv_rows(
        data_path(output_dir, "v3_protocol_cross_family_robustness_decisions.csv")
    )
    docs_text = "\n".join(
        [
            Path("README.md").read_text(encoding="utf-8"),
            Path("docs/theory/rejected_claims.md").read_text(encoding="utf-8"),
            Path("docs/theory/handoff_manifest_provenance.md").read_text(
                encoding="utf-8"
            ),
        ]
    )
    checks = [
        (
            "v3_protocol_decisions_file_exists",
            data_path(
                output_dir, "v3_protocol_cross_family_robustness_decisions.csv"
            ).exists(),
        ),
        ("v3_protocol_registry_json_serializes", bool(registry.get("registry_id"))),
        ("blocked_report_only_failed_control_blocked_in_plan", blocked_rows_blocked),
        (
            "no_v3_protocol_stress_test_result_output",
            not any((output_dir / "data").glob("v3_protocol_*stress_test_result*.csv")),
        ),
        (
            "diagnostic_complete_flag_exists",
            bool(decisions) and "diagnostic_complete" in decisions[0],
        ),
        (
            "preconditions_passed_flag_exists",
            bool(decisions) and "preconditions_passed" in decisions[0],
        ),
        (
            "top_down_hybrid_not_interpreted_as_evidence",
            (
                "Top-down/hybrid provenance is allowed only as "
                "preregistered manifest provenance"
            )
            in docs_text,
        ),
        (
            "registry_forbidden_interpretations_present",
            bool(registry.get("forbidden_interpretations")),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def main() -> None:
    output_dir = parse_args()
    path = write_csv(
        run_experiment(output_dir),
        data_path(output_dir, "v3_protocol_carry_forward_final_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v3 protocol carry-forward final sanity: {path}")


if __name__ == "__main__":
    main()
