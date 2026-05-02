"""Audit M44 for no retuning, no carry-forward, and no stress tests."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 protocol no-retuning audit.")
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


def _later_m45_registry_exists(output_dir: Path) -> bool:
    registry_path = (
        output_dir
        / "carry_forward_v4_protocol"
        / "carry_forward_registry_v4_protocol.json"
    )
    if not registry_path.exists():
        return False
    try:
        payload = json.loads(registry_path.read_text())
    except json.JSONDecodeError:
        return False
    return payload.get("created_by_milestone") == "Milestone 45"


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    data_dir = config.output_dir / "data"
    later_m45_registry_exists = _later_m45_registry_exists(config.output_dir)
    carry_forward_decision_exists = data_path(
        config.output_dir,
        "v4_protocol_cross_family_robustness_decisions.csv",
    ).exists()
    carry_forward_registry_exists = (
        (config.output_dir / "carry_forward_v4").exists()
        or (config.output_dir / "carry_forward_v4_protocol").exists()
    )
    checks = [
        (
            "m43_v4_preregistration_exists",
            (
                config.output_dir
                / "remediation"
                / "v4_protocol_preregistration_spec_m43.json"
            ).exists(),
        ),
        (
            "v4_manifests_written_under_manifests_v4",
            bool(list((config.output_dir / "manifests_v4").glob("*.json"))),
        ),
        (
            "no_v4_carry_forward_decision_file",
            not carry_forward_decision_exists or later_m45_registry_exists,
        ),
        (
            "no_v4_carry_forward_registry",
            not carry_forward_registry_exists or later_m45_registry_exists,
        ),
        (
            "no_v4_stress_test_result_output",
            not any(data_dir.glob("v4_*stress_test_result*.csv")),
        ),
        (
            "fixed_m34_criteria_unchanged",
            _criteria_match_export_if_present(config.output_dir),
        ),
        (
            "m43_no_execution_audit_not_retroactively_invalidated",
            data_path(config.output_dir, "v4_no_execution_audit.csv").exists(),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v4_protocol_no_retuning_audit.csv"),
        ["check", "passed"],
    )
    print(f"Wrote v4 protocol no-retuning audit: {path}")


if __name__ == "__main__":
    main()
