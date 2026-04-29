"""Export the frozen carry-forward registry."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from cross_family_robustness_experiment_helpers import (
    build_registry_from_output,
    registry_summary_rows,
)
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    write_carry_forward_registry,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for registry export."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Carry-forward registry export.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[Path, list[dict[str, float | str]]]:
    """Build and write carry-forward registry."""

    registry = build_registry_from_output(config.output_dir)
    registry_path = write_carry_forward_registry(
        registry,
        config.output_dir / "carry_forward" / "carry_forward_registry.json",
    )
    return registry_path, registry_summary_rows(registry)


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write registry export summary CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "carry_forward_registry_export.csv",
        [
            "registry_id",
            "family_name",
            "family_kind",
            "decision",
            "eligible_manifest_count",
            "ineligible_manifest_count",
            "allowed_future_use",
            "failed_reasons",
            "warning_reasons",
        ],
    )


def main() -> None:
    config = parse_args()
    registry_path, rows = run_experiment(config)
    csv_path = write_outputs(rows, config.output_dir)
    print(f"Wrote carry-forward registry export: {csv_path}")
    print(f"Wrote carry-forward registry JSON: {registry_path}")


if __name__ == "__main__":
    main()
