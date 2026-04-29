"""Held-out protocol validation for response-comparison constraint pools."""

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
    heldout_protocol_constraint_validation,
    validation_gate_pass_fail,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for held-out protocol validation."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32, 40)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    train_fractions: tuple[float, ...] = (0.5, 0.7)
    max_constraints: int = 5000
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

    parser = argparse.ArgumentParser(description="Held-out constraint validation.")
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
    parser.add_argument("--train-fractions", nargs="+", default=["0.5", "0.7"])
    parser.add_argument("--max-constraints", type=int, default=5000)
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
        train_fractions=_parse_float_values(args.train_fractions),
        max_constraints=args.max_constraints,
        min_margins=_parse_float_values(args.min_margins),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run held-out protocol-column validation."""

    rows: list[dict[str, float | str]] = []
    gate = ConstraintValidationGate("heldout_default")
    protocols = default_pairwise_protocols()
    for repetition in range(config.repetitions):
        profile = profile_from_synthetic_config(
            config.reference_length,
            config.emission_positions,
            config.layer_delay_ranks,
            config.targets_per_layer,
            repetition,
            config.seed,
        )
        for protocol in protocols:
            for train_fraction in config.train_fractions:
                for min_margin in config.min_margins:
                    row_seed = config.seed + 1000 * repetition + int(min_margin * 1000)
                    validation = heldout_protocol_constraint_validation(
                        profile,
                        protocol,
                        train_fraction=train_fraction,
                        max_constraints=config.max_constraints,
                        min_margin=min_margin,
                        seed=row_seed,
                    )
                    gate_row = validation_gate_pass_fail(validation, gate)
                    rows.append(
                        {
                            "repetition": float(repetition),
                            "protocol_name": protocol.name,
                            "method": protocol.method,
                            "train_fraction": float(train_fraction),
                            "min_margin": float(min_margin),
                            **validation,
                            **gate_row,
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write held-out validation CSV."""

    path = output_dir / "data" / "response_constraint_heldout_protocol_validation.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_by_protocol(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    labels = sorted({str(row["protocol_name"]) for row in rows})
    values = [
        float(
            np.nanmean(
                [float(row[key]) for row in rows if row["protocol_name"] == label]
            )
        )
        for label in labels
    ]
    return labels, values


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save held-out validation figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "agreement_fraction",
            "response_constraint_heldout_agreement.png",
            "held-out agreement fraction",
        ),
        (
            "passed",
            "response_constraint_heldout_gate_pass_rate.png",
            "gate pass rate",
        ),
    ]:
        labels, values = _mean_by_protocol(rows, key)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(labels, values)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", labelrotation=30)
        ax.grid(True, axis="y", alpha=0.3)
        path = figure_dir / filename
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
    print(f"Wrote response constraint held-out validation: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
