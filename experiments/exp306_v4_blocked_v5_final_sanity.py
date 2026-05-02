"""Final sanity checks for M46 v4 blocked/v5 preregistration."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import read_csv_rows
from v4_blocked_v5_experiment_helpers import data_path


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="M46 final sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def _json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _counterfactuals_report_only(output_dir: Path) -> bool:
    paths = [
        data_path(output_dir, "v4_heldout_counterfactual.csv"),
        data_path(output_dir, "v4_latent_order_counterfactual.csv"),
        data_path(output_dir, "v4_destructive_null_counterfactual.csv"),
        data_path(output_dir, "v4_single_fix_counterfactual.csv"),
    ]
    return all(
        rows
        and all(
            row.get("note", row.get("output_note", ""))
            == "report_only_not_decision_changing"
            for row in rows
        )
        for rows in [read_csv_rows(path) for path in paths]
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    prereg = _json(
        config.output_dir / "remediation" / "v5_protocol_preregistration_spec_m46.json"
    )
    decisions = read_csv_rows(
        data_path(
            config.output_dir,
            "v4_protocol_cross_family_robustness_decisions.csv",
        )
    )
    blocked_v4_remain_blocked = all(
        row.get("decision") == "blocked"
        for row in decisions
        if row.get("family_kind") == "structured"
    )
    data_dir = config.output_dir / "data"
    checks = [
        (
            "v4_blocking_audit_exists",
            data_path(config.output_dir, "v4_blocked_root_cause_audit.csv").exists(),
        ),
        (
            "criterion_margin_rows_exist",
            bool(
                read_csv_rows(
                    data_path(config.output_dir, "v4_criterion_margin_report.csv")
                )
            ),
        ),
        (
            "v3_to_v4_delta_rows_exist",
            bool(
                read_csv_rows(
                    data_path(config.output_dir, "v3_to_v4_metric_delta_audit.csv")
                )
            ),
        ),
        (
            "counterfactuals_report_only",
            _counterfactuals_report_only(config.output_dir),
        ),
        (
            "v5_execution_allowed_false",
            prereg.get("execution_allowed_in_current_milestone") is False,
        ),
        (
            "v5_remediation_iteration_risk_audit_exists",
            data_path(
                config.output_dir,
                "v5_remediation_iteration_risk_audit.csv",
            ).exists(),
        ),
        (
            "no_v5_production_manifests",
            not (config.output_dir / "manifests_v5").exists()
            or not list((config.output_dir / "manifests_v5").glob("*.json")),
        ),
        ("blocked_v4_families_remain_blocked", blocked_v4_remain_blocked),
        (
            "no_stress_test_output",
            not any(data_dir.glob("*stress_test_result*.csv")),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v4_blocked_v5_final_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote M46 final sanity: {path}")


if __name__ == "__main__":
    main()
