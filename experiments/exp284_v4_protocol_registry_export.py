"""Export the v4 protocol carry-forward registry."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_protocol_carry_forward_experiment_helpers import (
    data_path,
    decisions_from_csv,
)

from causal_spacetime_lab.state_change_manifest_v4_registry import (
    build_v4_protocol_carry_forward_registry,
    v4_protocol_registry_summary_rows,
    write_v4_protocol_carry_forward_registry,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    manifest_dir: Path = Path("outputs/manifests_v4")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Export v4 protocol registry.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--manifest-dir", type=Path, default=Path("outputs/manifests_v4")
    )
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir, args.manifest_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[Path | None, list[dict[str, float | str]]]:
    decisions = decisions_from_csv(
        data_path(
            config.output_dir, "v4_protocol_cross_family_robustness_decisions.csv"
        )
    )
    if not decisions:
        return None, [
            {
                "registry_id": "",
                "criteria_name": "fixed_m34_v4_protocol",
                "family_name": "",
                "family_kind": "",
                "decision": "missing_decisions",
                "manifest_count": 0.0,
                "eligible_manifest_count": 0.0,
                "ineligible_manifest_count": 0.0,
                "failed_reasons": "missing_v4_protocol_decisions",
                "warning_reasons": "",
                "allowed_future_use": "none",
            }
        ]
    registry = build_v4_protocol_carry_forward_registry(
        decisions,
        config.manifest_dir,
    )
    return (
        write_v4_protocol_carry_forward_registry(registry, config.output_dir),
        v4_protocol_registry_summary_rows(registry),
    )


def main() -> None:
    config = parse_args()
    registry_path, rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v4_protocol_carry_forward_registry_export.csv"),
        ["registry_id", "criteria_name", "family_name", "decision"],
    )
    print(f"Wrote v4 protocol registry summary: {path}")
    if registry_path is not None:
        print(f"Wrote v4 protocol registry: {registry_path}")


if __name__ == "__main__":
    main()
