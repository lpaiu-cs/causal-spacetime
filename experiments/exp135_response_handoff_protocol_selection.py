"""Compare protocols for eligible response handoff manifests."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import (
    default_pairwise_protocols,
    profile_from_synthetic_config,
)

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
    summarize_handoff_manifests,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for protocol handoff selection."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32, 40)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    min_margins: tuple[float, ...] = (0.0, 0.05, 0.10)
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Handoff protocol selection.")
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
    parser.add_argument("--min-margins", nargs="+", default=["0.0", "0.05", "0.10"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        min_margins=_parse_float_values(args.min_margins),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build handoff manifests across protocols and margins."""

    gate = ConstraintValidationGate("protocol_selection")
    manifests = []
    metadata: list[tuple[int, float]] = []
    for repetition in range(config.repetitions):
        profile = profile_from_synthetic_config(
            config.reference_length,
            config.emission_positions,
            config.layer_delay_ranks,
            config.targets_per_layer,
            repetition,
            config.seed,
        )
        for protocol in default_pairwise_protocols():
            for min_margin in config.min_margins:
                manifests.append(
                    build_candidate_handoff_manifest(
                        profile,
                        protocol,
                        gate,
                        max_constraints=2000,
                        min_margin=min_margin,
                        constraint_seed=config.seed + 1000 * repetition,
                        bootstrap_count=10,
                        null_repetitions=3,
                        source_label=f"protocol_selection_rep{repetition}",
                    )
                )
                metadata.append((repetition, min_margin))
    rows = summarize_handoff_manifests(manifests)
    for row, (repetition, min_margin) in zip(rows, metadata, strict=True):
        row["repetition"] = float(repetition)
        row["min_margin"] = float(min_margin)
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write protocol selection CSV."""

    path = output_dir / "data" / "response_handoff_protocol_selection.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save protocol selection figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    labels = sorted({str(row["protocol_name"]) for row in rows})
    pass_rates = [
        float(
            np.mean(
                [
                    float(row["eligible"])
                    for row in rows
                    if row["protocol_name"] == label
                ]
            )
        )
        for label in labels
    ]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, pass_rates)
    ax.set_ylabel("eligible manifest fraction")
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_handoff_protocol_eligibility.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    values = [
        float(
            np.nanmean(
                [
                    float(row["heldout_agreement"])
                    for row in rows
                    if row["protocol_name"] == label
                ]
            )
        )
        for label in labels
    ]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, values)
    ax.set_ylabel("held-out agreement")
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_handoff_metric_summary.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote response handoff protocol selection: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
