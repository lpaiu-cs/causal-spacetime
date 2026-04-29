"""Compare admissible pairwise response-profile protocols."""

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

from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    compare_pairwise_protocols,
    pairwise_protocol_admissibility_report,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for pairwise response protocol comparison."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
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

    parser = argparse.ArgumentParser(description="Pairwise response protocols.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument("--emission-positions", nargs="+", default=["8", "16", "24"])
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run pairwise protocol comparison diagnostics."""

    rows: list[dict[str, float | str]] = []
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
            rows.append(
                {
                    "row_type": "admissibility",
                    "repetition": float(repetition),
                    **pairwise_protocol_admissibility_report(profile, protocol),
                    "protocol_a": "",
                    "protocol_b": "",
                    "common_valid_pair_count": float("nan"),
                    "order_inversion_rate": float("nan"),
                    "tie_change_fraction": float("nan"),
                }
            )
        for comparison in compare_pairwise_protocols(profile, protocols):
            rows.append(
                {
                    "row_type": "protocol_comparison",
                    "repetition": float(repetition),
                    "protocol_name": "",
                    "method": "",
                    "missing_policy": "",
                    "target_count": float(profile.target_event_ids.size),
                    "valid_pair_count": float("nan"),
                    "valid_pair_fraction": float("nan"),
                    "finite_value_fraction": float("nan"),
                    "diagonal_zero": float("nan"),
                    "symmetric": float("nan"),
                    "mean_dissimilarity": float("nan"),
                    "tie_fraction": float("nan"),
                    "strict_pair_order_fraction": float("nan"),
                    **comparison,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "pairwise_response_protocol_comparison.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save protocol comparison figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    admissibility = [row for row in rows if row["row_type"] == "admissibility"]
    labels = sorted({str(row["protocol_name"]) for row in admissibility})
    valid_means = [
        float(
            np.nanmean(
                [
                    float(row["valid_pair_fraction"])
                    for row in admissibility
                    if row["protocol_name"] == label
                ]
            )
        )
        for label in labels
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels, valid_means)
    ax.set_ylabel("valid pair fraction")
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "pairwise_response_valid_pair_fraction.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    comparisons = [row for row in rows if row["row_type"] == "protocol_comparison"]
    pair_labels = [
        f"{row['protocol_a']} vs {row['protocol_b']}" for row in comparisons
    ]
    values = [float(row["order_inversion_rate"]) for row in comparisons]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(range(len(values)), values)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(pair_labels, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("order inversion rate")
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "pairwise_response_protocol_disagreement.png"
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
    print(f"Wrote pairwise response protocol comparison: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
