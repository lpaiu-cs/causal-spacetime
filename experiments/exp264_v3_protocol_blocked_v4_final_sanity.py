"""Final sanity checks for M43 v3 blocked-decision audit and v4 preregistration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_blocked_v4_experiment_helpers import data_path, read_csv_rows


def parse_args() -> Path:
    parser = argparse.ArgumentParser(description="M43 final sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args().output_dir


def _json_execution_allowed_false(path: Path) -> bool:
    if not path.exists():
        return False
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload.get("execution_allowed_in_current_milestone") is False


def _blocked_v3_families_remain_blocked(output_dir: Path) -> bool:
    rows = read_csv_rows(
        data_path(output_dir, "v3_protocol_cross_family_robustness_decisions.csv")
    )
    structured_rows = [
        row for row in rows if row.get("family_kind") == "structured"
    ]
    return bool(structured_rows) and all(
        row.get("decision") == "blocked" for row in structured_rows
    )


def _counterfactuals_report_only(output_dir: Path) -> bool:
    filenames = [
        "v3_protocol_heldout_counterfactual.csv",
        "v3_protocol_latent_order_counterfactual.csv",
        "v3_protocol_single_fix_counterfactual.csv",
    ]
    for filename in filenames:
        rows = read_csv_rows(data_path(output_dir, filename))
        if not rows:
            return False
        if any(
            row.get("output_note") != "report_only_not_decision_changing"
            for row in rows
        ):
            return False
    return True


def run_experiment(output_dir: Path = Path("outputs")) -> list[dict[str, float | str]]:
    prereg_path = (
        output_dir / "remediation" / "v4_protocol_preregistration_spec_m43.json"
    )
    data_dir = output_dir / "data"
    checks = [
        (
            "v3_blocking_audit_exists",
            data_path(output_dir, "v3_protocol_blocked_root_cause_audit.csv").exists(),
        ),
        (
            "criterion_margin_rows_exist",
            bool(
                read_csv_rows(
                    data_path(output_dir, "v3_protocol_criterion_margin_report.csv")
                )
            ),
        ),
        ("counterfactuals_report_only", _counterfactuals_report_only(output_dir)),
        (
            "v4_preregistration_execution_not_allowed",
            _json_execution_allowed_false(prereg_path),
        ),
        (
            "no_v4_production_manifests",
            not (output_dir / "manifests_v4").exists()
            or not list((output_dir / "manifests_v4").glob("*.json")),
        ),
        (
            "blocked_v3_families_remain_blocked",
            _blocked_v3_families_remain_blocked(output_dir),
        ),
        (
            "no_stress_test_result_output",
            not any(data_dir.glob("*stress_test_result*.csv")),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def main() -> None:
    output_dir = parse_args()
    path = write_csv(
        run_experiment(output_dir),
        data_path(output_dir, "v3_protocol_blocked_v4_final_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote M43 final sanity: {path}")


if __name__ == "__main__":
    main()
