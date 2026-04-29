"""Null baselines for pairwise response-profile comparison."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import (
    constraint_count,
    profile_from_synthetic_config,
)

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
    pairwise_response_order_inversion_rate,
)
from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    pairwise_protocol_admissibility_report,
)
from causal_spacetime_lab.state_change_response_pairwise_nulls import (
    permute_target_profiles,
    random_profile_with_same_marginals,
    shuffle_profile_delays_within_protocols,
    shuffle_profile_reachability_masks,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for pairwise response null baselines."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    repetitions: int = 5
    null_repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Pairwise response null baselines.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument("--emission-positions", nargs="+", default=["8", "16", "24"])
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--null-repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        repetitions=args.repetitions,
        null_repetitions=args.null_repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run null baseline diagnostics."""

    protocol = PairwiseResponseComparisonProtocol("gap_mean", "rank_gap_mean")
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
        structured = pairwise_response_dissimilarity(profile, protocol)
        profiles = [("structured", profile)]
        for null_rep in range(config.null_repetitions):
            seed = config.seed + 1000 * repetition + null_rep
            profiles.extend(
                [
                    (
                        "shuffle_delays",
                        shuffle_profile_delays_within_protocols(profile, seed),
                    ),
                    (
                        "shuffle_reachability",
                        shuffle_profile_reachability_masks(profile, seed + 11),
                    ),
                    (
                        "permute_profiles",
                        permute_target_profiles(profile, seed=seed + 22),
                    ),
                    (
                        "random_same_marginals",
                        random_profile_with_same_marginals(profile, seed + 33),
                    ),
                ]
            )
        for baseline_type, candidate in profiles:
            report = pairwise_protocol_admissibility_report(candidate, protocol)
            dissimilarity = pairwise_response_dissimilarity(candidate, protocol)
            rows.append(
                {
                    "baseline_type": baseline_type,
                    "repetition": float(repetition),
                    "order_inversion_vs_structured": (
                        pairwise_response_order_inversion_rate(
                            structured,
                            dissimilarity,
                        )
                    ),
                    "constraint_count": constraint_count(
                        candidate,
                        protocol,
                        config.seed + repetition,
                    ),
                    **report,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "pairwise_response_null_baselines.csv"
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
    labels = sorted({str(row["baseline_type"]) for row in rows})
    values = [
        float(
            np.nanmean(
                [
                    float(row[key])
                    for row in rows
                    if row["baseline_type"] == label
                ]
            )
        )
        for label in labels
    ]
    return labels, values


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save null baseline figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "strict_pair_order_fraction",
            "pairwise_response_null_strict_fraction.png",
            "strict pair order fraction",
        ),
        (
            "constraint_count",
            "pairwise_response_null_constraint_count.png",
            "constraint count",
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
    print(f"Wrote pairwise response null baselines: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
