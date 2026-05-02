"""Verify M46 did not execute v5."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_blocked_v5_experiment_helpers import data_path, read_json

from causal_spacetime_lab.state_change_manifest_v5_audit import (
    audit_v5_preregistration_only,
    check_no_v5_execution_outputs,
)
from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v5_preregistration import (
    build_v5_protocol_preregistration_spec,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V5 no-execution audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    prereg_path = (
        config.output_dir / "remediation" / "v5_protocol_preregistration_spec_m46.json"
    )
    if not read_json(prereg_path):
        return [
            {
                "check": "v5_preregistration_exists",
                "passed": 0.0,
                "detail": "missing_v5_preregistration_spec",
            }
        ]
    spec = build_v5_protocol_preregistration_spec(
        default_v5_protocol_family_designs()
    )
    audit = audit_v5_preregistration_only(spec)
    no_outputs = check_no_v5_execution_outputs(config.output_dir)
    rows = [
        {"check": key, "passed": float(value), "detail": ""}
        for key, value in {**audit, **no_outputs}.items()
    ]
    rows.insert(0, {"check": "v5_preregistration_exists", "passed": 1.0, "detail": ""})
    return rows


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v5_no_execution_audit.csv"),
        ["check", "passed", "detail"],
    )
    print(f"Wrote v5 no-execution audit: {path}")


if __name__ == "__main__":
    main()
