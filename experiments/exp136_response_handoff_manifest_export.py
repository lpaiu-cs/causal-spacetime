"""Export eligible response-comparison handoff manifests."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from pairwise_response_experiment_helpers import (
    default_pairwise_protocols,
    profile_from_synthetic_config,
)

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import write_handoff_manifest
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
    select_eligible_manifests,
    summarize_handoff_manifests,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for handoff manifest export."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32, 40)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    max_manifests: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Export handoff manifests.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument(
        "--emission-positions",
        nargs="+",
        default=["8", "16", "24", "32", "40"],
    )
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument("--max-manifests", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        max_manifests=args.max_manifests,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[Path]]:
    """Build and export eligible handoff manifests."""

    profile = profile_from_synthetic_config(
        config.reference_length,
        config.emission_positions,
        config.layer_delay_ranks,
        config.targets_per_layer,
        repetition=0,
        seed=config.seed,
    )
    manifests = []
    for protocol in default_pairwise_protocols():
        for min_margin in (0.0, 0.05, 0.10):
            manifests.append(
                build_candidate_handoff_manifest(
                    profile,
                    protocol,
                    ConstraintValidationGate("export"),
                    max_constraints=2000,
                    min_margin=min_margin,
                    constraint_seed=config.seed,
                    bootstrap_count=10,
                    null_repetitions=3,
                    source_label="manifest_export",
                )
            )
    eligible = select_eligible_manifests(manifests)[: config.max_manifests]
    manifest_dir = config.output_dir / "manifests"
    paths = [
        write_handoff_manifest(
            manifest,
            manifest_dir / f"response_handoff_{manifest.manifest_id[:12]}.json",
        )
        for manifest in eligible
    ]
    rows = summarize_handoff_manifests(eligible)
    for row, path in zip(rows, paths, strict=True):
        row["manifest_path"] = str(path)
    return rows, paths


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write manifest export CSV."""

    path = output_dir / "data" / "response_handoff_manifest_export.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else [
        "manifest_id",
        "protocol_name",
        "method",
        "min_margin",
        "eligible",
        "failed_reasons",
        "constraint_count",
        "heldout_agreement",
        "bootstrap_agreement",
        "null_z_score",
        "target_coverage",
        "pair_node_coverage",
        "manifest_path",
    ]
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    config = parse_args()
    rows, paths = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    print(f"Wrote response handoff manifest export: {data_path}")
    for path in paths:
        print(f"Wrote manifest: {path}")


if __name__ == "__main__":
    main()
