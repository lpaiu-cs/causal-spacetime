"""Verify M46 did not retune thresholds, rerun fits, or run stress tests."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import read_csv_rows
from v4_blocked_v5_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="M46 no-retuning audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def _criteria_match_export_if_present(output_dir: Path) -> bool:
    rows = read_csv_rows(
        data_path(output_dir, "cross_family_robustness_criteria_table.csv")
    )
    if not rows:
        return True
    exported = {row["criterion"]: row["value"] for row in rows}
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


def _counterfactuals_report_only(output_dir: Path) -> bool:
    paths = [
        data_path(output_dir, "v4_heldout_counterfactual.csv"),
        data_path(output_dir, "v4_latent_order_counterfactual.csv"),
        data_path(output_dir, "v4_destructive_null_counterfactual.csv"),
        data_path(output_dir, "v4_single_fix_counterfactual.csv"),
    ]
    for path in paths:
        rows = read_csv_rows(path)
        if not rows:
            return False
        if not all(
            row.get("note", row.get("output_note", ""))
            == "report_only_not_decision_changing"
            for row in rows
        ):
            return False
    return True


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    data_dir = config.output_dir / "data"
    checks = [
        (
            "m45_decision_files_remain_inputs",
            data_path(
                config.output_dir,
                "v4_protocol_cross_family_robustness_decisions.csv",
            ).exists(),
        ),
        (
            "no_v5_production_manifests",
            not (config.output_dir / "manifests_v5").exists()
            or not list((config.output_dir / "manifests_v5").glob("*.json")),
        ),
        ("no_v5_fit_outputs", not any(data_dir.glob("v5_*fit*.csv"))),
        (
            "no_v5_carry_forward_registry",
            not (config.output_dir / "carry_forward_v5").exists()
            and not (config.output_dir / "carry_forward_v5_protocol").exists(),
        ),
        (
            "no_stress_test_outputs",
            not any(data_dir.glob("*stress_test_result*.csv")),
        ),
        (
            "counterfactual_outputs_report_only",
            _counterfactuals_report_only(config.output_dir),
        ),
        (
            "fixed_m34_criteria_unchanged",
            _criteria_match_export_if_present(config.output_dir),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v4_blocked_no_retuning_audit.csv"),
        ["check", "passed"],
    )
    print(f"Wrote M46 no-retuning audit: {path}")


if __name__ == "__main__":
    main()
