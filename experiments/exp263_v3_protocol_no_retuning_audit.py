"""Audit that M43 did not retune, refit, change decisions, or run stress tests."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_blocked_v4_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="M43 no-retuning audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def _criteria_export_matches_default(output_dir: Path) -> bool:
    rows = read_csv_rows(
        data_path(output_dir, "cross_family_robustness_criteria_table.csv")
    )
    if not rows:
        return True
    exported = {row.get("criterion", ""): row.get("value", "") for row in rows}
    for key, value in asdict(default_cross_family_robustness_criteria()).items():
        if key not in exported:
            return False
        if str(exported[key]) == str(value):
            continue
        try:
            if float(exported[key]) != float(value):
                return False
        except (TypeError, ValueError):
            return False
    return True


def _all_counterfactual_rows_report_only(output_dir: Path) -> bool:
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    data_dir = config.output_dir / "data"
    checks = [
        (
            "m42_decision_files_remain_inputs",
            data_path(
                config.output_dir,
                "v3_protocol_cross_family_robustness_decisions.csv",
            ).exists()
            and data_path(
                config.output_dir,
                "v3_protocol_cross_family_robustness_decision_metrics.csv",
            ).exists(),
        ),
        (
            "no_v4_production_manifests",
            not (config.output_dir / "manifests_v4").exists()
            or not list((config.output_dir / "manifests_v4").glob("*.json")),
        ),
        ("no_v4_fit_outputs", not any(data_dir.glob("v4_*fit*.csv"))),
        (
            "no_v4_carry_forward_registry",
            not (config.output_dir / "carry_forward_v4").exists()
            and not (config.output_dir / "carry_forward_v4_protocol").exists(),
        ),
        (
            "no_stress_test_result_outputs",
            not any(data_dir.glob("*stress_test_result*.csv")),
        ),
        (
            "counterfactuals_report_only",
            _all_counterfactual_rows_report_only(config.output_dir),
        ),
        (
            "fixed_m34_criteria_unchanged",
            _criteria_export_matches_default(config.output_dir),
        ),
    ]
    return [
        {
            "check": name,
            "passed": float(passed),
            "detail": "M43 read-only blocked-decision audit",
        }
        for name, passed in checks
    ]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v3_protocol_blocked_no_retuning_audit.csv"),
        ["check", "passed", "detail"],
    )
    print(f"Wrote M43 no-retuning audit: {path}")


if __name__ == "__main__":
    main()
