"""Verify M40 protocol patch did not execute v3 generation."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_audit import (
    check_no_v3_execution_outputs,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="M40 no-execution audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> list[dict[str, float | str]]:
    checks = check_no_v3_execution_outputs(config.output_dir)
    prereg_path = (
        config.output_dir
        / "remediation"
        / "v3_protocol_patched_preregistration_m40.json"
    )
    prereg_disallows_execution = False
    if prereg_path.exists():
        payload = json.loads(prereg_path.read_text(encoding="utf-8"))
        prereg_disallows_execution = not bool(
            payload.get("execution_allowed_in_current_milestone", True)
        )
    rows = [
        {"check": key, "passed": float(value)}
        for key, value in checks.items()
    ]
    rows.append(
        {
            "check": "patched_preregistration_disallows_execution",
            "passed": float(prereg_disallows_execution),
        }
    )
    return rows


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "protocol_patch_no_execution_audit.csv"),
        ["check", "passed"],
    )
    print(f"Wrote protocol patch no-execution audit: {path}")


if __name__ == "__main__":
    main()
