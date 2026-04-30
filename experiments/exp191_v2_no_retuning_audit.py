"""Audit that M37 used fixed M36 inputs and avoided carry-forward decisions."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    data_path,
    remediation_path,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 no-retuning audit."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 no-retuning audit.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run no-retuning and no-stress audit checks."""

    output_dir = config.output_dir
    checks = [
        (
            "remediation_plan_exists",
            remediation_path(output_dir, "remediation_plan_m36.json").exists(),
        ),
        (
            "future_run_spec_exists",
            remediation_path(output_dir, "future_manifest_run_spec_m36.json").exists(),
        ),
        ("v2_manifest_dir_exists", (output_dir / "manifests_v2").exists()),
        (
            "v2_manifests_do_not_replace_v1_dir",
            (output_dir / "manifests_v2").resolve()
            != (output_dir / "manifests").resolve(),
        ),
        (
            "no_v2_carry_forward_decision_file",
            not data_path(
                output_dir,
                "v2_cross_family_robustness_decisions.csv",
            ).exists(),
        ),
        (
            "no_v2_stress_test_output",
            not any((output_dir / "data").glob("v2_*stress*.csv")),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 no-retuning audit rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_no_retuning_audit.csv"),
        ["check", "passed"],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    print(f"Wrote v2 no-retuning audit: {output_path}")


if __name__ == "__main__":
    main()
