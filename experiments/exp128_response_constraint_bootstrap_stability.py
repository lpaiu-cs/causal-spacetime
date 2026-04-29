"""Bootstrap stability for response-comparison constraint pools."""

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

from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_constraint_validation import (
    bootstrap_constraint_stability,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    pairwise_response_dissimilarity,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for bootstrap stability."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32, 40)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    max_constraints: int = 5000
    min_margin: float = 0.05
    bootstrap_count: int = 50
    sample_fraction: float = 0.8
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Constraint bootstrap stability.")
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
    parser.add_argument("--max-constraints", type=int, default=5000)
    parser.add_argument("--min-margin", type=float, default=0.05)
    parser.add_argument("--bootstrap-count", type=int, default=50)
    parser.add_argument("--sample-fraction", type=float, default=0.8)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        max_constraints=args.max_constraints,
        min_margin=args.min_margin,
        bootstrap_count=args.bootstrap_count,
        sample_fraction=args.sample_fraction,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run bootstrap stability diagnostics."""

    rows: list[dict[str, float | str]] = []
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
            dissimilarity = pairwise_response_dissimilarity(profile, protocol)
            pool = build_constraint_pool_from_dissimilarity(
                dissimilarity,
                max_constraints=config.max_constraints,
                min_margin=config.min_margin,
                seed=config.seed + repetition,
            )
            rows.append(
                {
                    "repetition": float(repetition),
                    "protocol_name": protocol.name,
                    "method": protocol.method,
                    "constraint_count": float(pool.constraints.shape[0]),
                    "min_margin": float(config.min_margin),
                    **bootstrap_constraint_stability(
                        profile,
                        protocol,
                        pool,
                        bootstrap_count=config.bootstrap_count,
                        sample_fraction=config.sample_fraction,
                        seed=config.seed + 1000 * repetition,
                    ),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write bootstrap stability CSV."""

    path = output_dir / "data" / "response_constraint_bootstrap_stability.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _means(
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
    """Save bootstrap stability figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "mean_agreement_fraction",
            "response_constraint_bootstrap_agreement.png",
            "mean bootstrap agreement",
        ),
        (
            "stable_constraint_fraction",
            "response_constraint_stable_fraction.png",
            "stable constraint fraction",
        ),
    ]:
        labels, values = _means(rows, key)
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
    print(f"Wrote response constraint bootstrap stability: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
