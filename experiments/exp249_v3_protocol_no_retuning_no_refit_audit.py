"""Audit M42 for no retuning, no refit, no regeneration, and no stress tests."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol no-retuning audit.")
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    data_dir = config.output_dir / "data"
    checks = [
        (
            "m41_v3_protocol_bundle_exists",
            data_path(
                config.output_dir,
                "v3_protocol_cross_family_robustness_metrics.csv",
            ).exists(),
        ),
        (
            "m42_did_not_write_new_manifest_generation_file",
            not data_path(
                config.output_dir, "v3_protocol_manifest_generation_m42.csv"
            ).exists(),
        ),
        (
            "m42_did_not_write_refit_output",
            not data_path(
                config.output_dir,
                "v3_protocol_manifest_family_fit_summary_m42.csv",
            ).exists(),
        ),
        (
            "no_v3_protocol_stress_test_result_output",
            not any(data_dir.glob("v3_protocol_*stress_test_result*.csv")),
        ),
        (
            "threshold_sensitivity_labeled_sensitivity_only",
            all(
                row.get("analysis_label") == "sensitivity_only_not_retuning"
                for row in read_csv_rows(
                    data_path(
                        config.output_dir, "v3_protocol_threshold_sensitivity.csv"
                    )
                )
            ),
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
        data_path(config.output_dir, "v3_protocol_carry_forward_no_retuning_audit.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v3 protocol no-retuning/no-refit audit: {path}")


if __name__ == "__main__":
    main()
