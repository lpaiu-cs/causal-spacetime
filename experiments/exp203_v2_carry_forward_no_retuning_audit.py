"""Audit M38 for no threshold retuning, no refit, and no stress-test run."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 no-retuning audit."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 no-retuning/no-stress audit.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _criteria_table_matches_default(output_dir: Path) -> bool:
    path = data_path(output_dir, "cross_family_robustness_criteria_table.csv")
    rows = read_csv_rows(path)
    if not rows:
        return True
    exported = {row["criterion"]: row["value"] for row in rows}
    default = asdict(default_cross_family_robustness_criteria())
    for key, value in default.items():
        if key not in exported:
            return False
        if str(exported[key]) != str(value):
            try:
                if float(exported[key]) != float(value):
                    return False
            except (TypeError, ValueError):
                return False
    return True


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run v2 no-retuning/no-stress audit checks."""

    data_dir = config.output_dir / "data"
    checks = [
        (
            "fixed_criteria_match_export_if_present",
            _criteria_table_matches_default(config.output_dir),
        ),
        (
            "v2_decisions_exist",
            data_path(
                config.output_dir,
                "v2_cross_family_robustness_decisions.csv",
            ).exists(),
        ),
        (
            "m37_v2_manifests_still_present",
            (config.output_dir / "manifests_v2").exists()
            and bool(list((config.output_dir / "manifests_v2").glob("*.json"))),
        ),
        (
            "no_m38_manifest_generation_output",
            not data_path(config.output_dir, "v2_manifest_generation_m38.csv").exists(),
        ),
        (
            "no_v2_stress_test_result_output",
            not any(data_dir.glob("v2_*stress_test_result*.csv")),
        ),
        (
            "threshold_sensitivity_labeled_sensitivity",
            all(
                row.get("analysis_label") == "sensitivity_only_not_threshold_retuning"
                for row in read_csv_rows(
                    data_path(
                        config.output_dir,
                        "v2_carry_forward_threshold_sensitivity.csv",
                    )
                )
            ),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write v2 no-retuning audit rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_carry_forward_no_retuning_audit.csv"),
        ["check", "passed"],
    )


def main() -> None:
    config = parse_args()
    path = write_outputs(run_experiment(config), config.output_dir)
    print(f"Wrote v2 carry-forward no-retuning audit: {path}")


if __name__ == "__main__":
    main()
