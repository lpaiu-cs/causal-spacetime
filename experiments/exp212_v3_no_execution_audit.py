"""Audit that Milestone 39 did not execute v3."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_audit import (
    audit_v3_preregistration_only,
    check_no_v3_execution_outputs,
)
from causal_spacetime_lab.state_change_manifest_v3_design import (
    V3ManifestFamilyDesign,
    default_v3_manifest_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v3_preregistration import (
    V3PreregistrationSpec,
    build_v3_preregistration_spec,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 no-execution audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _read_spec(output_dir: Path) -> V3PreregistrationSpec:
    path = output_dir / "remediation" / "v3_preregistration_spec_m39.json"
    if not path.exists():
        return build_v3_preregistration_spec(default_v3_manifest_family_designs())
    payload = json.loads(path.read_text(encoding="utf-8"))
    families = [
        V3ManifestFamilyDesign(**row)
        for row in payload.get("planned_families", [])
        if isinstance(row, dict)
    ]
    return V3PreregistrationSpec(
        spec_id=str(payload.get("spec_id", "")),
        created_by_milestone=str(payload.get("created_by_milestone", "")),
        based_on_outputs=list(payload.get("based_on_outputs", [])),
        fixed_criteria_source=str(payload.get("fixed_criteria_source", "")),
        planned_families=families,
        required_metrics=list(payload.get("required_metrics", [])),
        required_output_files=list(payload.get("required_output_files", [])),
        forbidden_actions=list(payload.get("forbidden_actions", [])),
        forbidden_interpretations=list(payload.get("forbidden_interpretations", [])),
        execution_allowed_in_current_milestone=bool(
            payload.get("execution_allowed_in_current_milestone", True)
        ),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    spec = _read_spec(config.output_dir)
    rows = audit_v3_preregistration_only(spec)
    rows.update(check_no_v3_execution_outputs(config.output_dir))
    return [{"check": key, "passed": value} for key, value in rows.items()]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v3_no_execution_audit.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v3 no-execution audit: {path}")


if __name__ == "__main__":
    main()

