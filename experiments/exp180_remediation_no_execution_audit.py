"""Audit that Milestone 36 did not execute remediation."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from remediation_plan_experiment_helpers import (
    read_or_build_remediation_plan,
    remediation_plan_from_jsonable,
    write_table,
)

from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    FutureManifestRunSpec,
)
from causal_spacetime_lab.state_change_manifest_remediation_audit import (
    future_run_spec_audit,
    remediation_plan_execution_audit,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for no-execution audit."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Remediation no-execution audit.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _read_plan(output_dir: Path):
    path = output_dir / "remediation" / "remediation_plan_m36.json"
    if path.exists():
        return remediation_plan_from_jsonable(json.loads(path.read_text()))
    return read_or_build_remediation_plan(output_dir)


def _read_spec(output_dir: Path) -> FutureManifestRunSpec | None:
    path = output_dir / "remediation" / "future_manifest_run_spec_m36.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    return FutureManifestRunSpec(**payload)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run no-execution audit checks."""

    plan = _read_plan(config.output_dir)
    spec = _read_spec(config.output_dir)
    plan_audit = remediation_plan_execution_audit(plan)
    spec_audit = future_run_spec_audit(spec) if spec else {"passed": 0.0}
    remediation_dir = config.output_dir / "remediation"
    production_manifest_exports = list(remediation_dir.glob("response_handoff_*.json"))
    new_fit_outputs = list(remediation_dir.glob("*fit*.json"))
    stress_outputs = list((config.output_dir / "data").glob("remediation_stress_*.csv"))
    checks = [
        ("no_new_production_handoff_manifests", not production_manifest_exports),
        ("no_new_representation_fit_outputs", not new_fit_outputs),
        ("no_stress_test_outputs", not stress_outputs),
        ("remediation_plan_execution_not_allowed", plan_audit["passed"] == 1.0),
        ("future_run_spec_execution_not_allowed", spec_audit["passed"] == 1.0),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write no-execution audit CSV."""

    return write_table(
        rows,
        output_dir,
        "remediation_no_execution_audit.csv",
        ["check", "passed"],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote remediation no-execution audit: {path}")


if __name__ == "__main__":
    main()
