"""Protocol-choice dependence for pairwise response-profile comparison."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import (
    default_pairwise_protocols,
    synthetic_random_profile,
)

from causal_spacetime_lab.state_change_response_pairwise import (
    pairwise_response_dissimilarity,
    pairwise_response_order_inversion_rate,
    response_pair_comparison_constraints,
)
from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    compare_pairwise_protocols,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for protocol-choice dependence diagnostics."""

    target_count: int = 40
    protocol_count: int = 4
    reachable_probability: float = 0.8
    unique_rank_count: int = 6
    repetitions: int = 10
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Pairwise protocol-choice dependence.")
    parser.add_argument("--target-count", type=int, default=40)
    parser.add_argument("--protocol-count", type=int, default=4)
    parser.add_argument("--reachable-probability", type=float, default=0.8)
    parser.add_argument("--unique-rank-count", type=int, default=6)
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        target_count=args.target_count,
        protocol_count=args.protocol_count,
        reachable_probability=args.reachable_probability,
        unique_rank_count=args.unique_rank_count,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _constraint_disagreement(
    constraints_a: np.ndarray,
    constraints_b: np.ndarray,
) -> float:
    set_a = {tuple(int(value) for value in row) for row in constraints_a}
    set_b = {tuple(int(value) for value in row) for row in constraints_b}
    if not set_a and not set_b:
        return 0.0
    union = set_a | set_b
    intersection = set_a & set_b
    return float(1.0 - len(intersection) / len(union))


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run protocol-choice dependence diagnostics."""

    rows: list[dict[str, float | str]] = []
    protocols = default_pairwise_protocols()
    for repetition in range(config.repetitions):
        profile = synthetic_random_profile(
            config.target_count,
            config.protocol_count,
            config.reachable_probability,
            config.unique_rank_count,
            config.seed + repetition,
        )
        comparisons = compare_pairwise_protocols(profile, protocols)
        dissimilarities = {
            protocol.name: pairwise_response_dissimilarity(profile, protocol)
            for protocol in protocols
        }
        constraints = {
            name: response_pair_comparison_constraints(
                dissimilarity,
                num_constraints=200,
                seed=config.seed + repetition,
            )
            for name, dissimilarity in dissimilarities.items()
        }
        for comparison in comparisons:
            protocol_a = str(comparison["protocol_a"])
            protocol_b = str(comparison["protocol_b"])
            rows.append(
                {
                    "repetition": float(repetition),
                    **comparison,
                    "direct_inversion_rate": pairwise_response_order_inversion_rate(
                        dissimilarities[protocol_a],
                        dissimilarities[protocol_b],
                    ),
                    "constraint_disagreement": _constraint_disagreement(
                        constraints[protocol_a],
                        constraints[protocol_b],
                    ),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "pairwise_response_protocol_choice_dependence.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save protocol-choice figure."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = [f"{row['protocol_a']} vs {row['protocol_b']}" for row in rows]
    values = [float(row["order_inversion_rate"]) for row in rows]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(range(len(values)), values)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("order inversion rate")
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "pairwise_response_protocol_choice_inversion.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return [path]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote pairwise response protocol-choice dependence: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
