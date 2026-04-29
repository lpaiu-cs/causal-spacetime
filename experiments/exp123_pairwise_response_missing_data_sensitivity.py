"""Missing-data sensitivity for pairwise response-profile comparison."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import synthetic_random_profile

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)
from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    pairwise_protocol_admissibility_report,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for missing-data sensitivity diagnostics."""

    target_count: int = 60
    protocol_count: int = 5
    reachable_probabilities: tuple[float, ...] = (0.3, 0.5, 0.7, 1.0)
    unique_rank_count: int = 5
    repetitions: int = 10
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Pairwise response missing data.")
    parser.add_argument("--target-count", type=int, default=60)
    parser.add_argument("--protocol-count", type=int, default=5)
    parser.add_argument(
        "--reachable-probabilities",
        nargs="+",
        default=["0.3", "0.5", "0.7", "1.0"],
    )
    parser.add_argument("--unique-rank-count", type=int, default=5)
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        target_count=args.target_count,
        protocol_count=args.protocol_count,
        reachable_probabilities=_parse_float_values(args.reachable_probabilities),
        unique_rank_count=args.unique_rank_count,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _protocols() -> list[PairwiseResponseComparisonProtocol]:
    return [
        PairwiseResponseComparisonProtocol(
            "common_reachable",
            "rank_gap_mean",
            missing_policy="common_reachable",
        ),
        PairwiseResponseComparisonProtocol(
            "require_all_reachable",
            "rank_gap_mean",
            missing_policy="require_all_reachable",
        ),
        PairwiseResponseComparisonProtocol(
            "penalize_mismatch",
            "combined_gap_and_mismatch",
            missing_policy="penalize_mismatch",
        ),
    ]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run missing-data sensitivity diagnostics."""

    rows: list[dict[str, float | str]] = []
    for probability in config.reachable_probabilities:
        for repetition in range(config.repetitions):
            profile = synthetic_random_profile(
                config.target_count,
                config.protocol_count,
                probability,
                config.unique_rank_count,
                config.seed + repetition + int(probability * 1000),
            )
            for protocol in _protocols():
                rows.append(
                    {
                        "reachable_probability": float(probability),
                        "repetition": float(repetition),
                        **pairwise_protocol_admissibility_report(profile, protocol),
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "pairwise_response_missing_data_sensitivity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_rows(
    rows: list[dict[str, float | str]],
    protocol: str,
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["reachable_probability"]) for row in rows})
    values = [
        float(
            np.nanmean(
                [
                    float(row[key])
                    for row in rows
                    if row["protocol_name"] == protocol
                    and float(row["reachable_probability"]) == probability
                ]
            )
        )
        for probability in probabilities
    ]
    return probabilities, values


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save missing-data figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    protocols = sorted({str(row["protocol_name"]) for row in rows})
    for key, filename, ylabel in [
        (
            "valid_pair_fraction",
            "pairwise_response_missing_valid_pairs.png",
            "valid pair fraction",
        ),
        ("tie_fraction", "pairwise_response_missing_ties.png", "tie fraction"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for protocol in protocols:
            xs, ys = _mean_rows(rows, protocol, key)
            ax.plot(xs, ys, marker="o", label=protocol)
        ax.set_xlabel("reachable probability")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend()
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
    print(f"Wrote pairwise response missing-data sensitivity: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
