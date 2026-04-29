"""Coverage diagnostics for response-comparison constraint pools."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import synthetic_random_profile

from causal_spacetime_lab.state_change_response_constraint_coverage import (
    constraint_pool_summary,
)
from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for constraint-pool coverage diagnostics."""

    target_counts: tuple[int, ...] = (20, 40, 80)
    protocol_counts: tuple[int, ...] = (3, 5, 8)
    unique_rank_count: int = 6
    reachable_probability: float = 0.8
    max_constraints_values: tuple[int, ...] = (500, 2000, 8000)
    min_margins: tuple[float, ...] = (0.0, 0.05)
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

    parser = argparse.ArgumentParser(description="Constraint-pool coverage.")
    parser.add_argument("--target-counts", nargs="+", default=["20", "40", "80"])
    parser.add_argument("--protocol-counts", nargs="+", default=["3", "5", "8"])
    parser.add_argument("--unique-rank-count", type=int, default=6)
    parser.add_argument("--reachable-probability", type=float, default=0.8)
    parser.add_argument(
        "--max-constraints-values",
        nargs="+",
        default=["500", "2000", "8000"],
    )
    parser.add_argument("--min-margins", nargs="+", default=["0.0", "0.05"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        target_counts=_parse_int_values(args.target_counts),
        protocol_counts=_parse_int_values(args.protocol_counts),
        unique_rank_count=args.unique_rank_count,
        reachable_probability=args.reachable_probability,
        max_constraints_values=_parse_int_values(args.max_constraints_values),
        min_margins=_parse_float_values(args.min_margins),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run constraint-pool coverage diagnostics."""

    protocol = PairwiseResponseComparisonProtocol("gap_mean", "rank_gap_mean")
    rows: list[dict[str, float | str]] = []
    for target_count in config.target_counts:
        for protocol_count in config.protocol_counts:
            for max_constraints in config.max_constraints_values:
                for min_margin in config.min_margins:
                    for repetition in range(config.repetitions):
                        profile = synthetic_random_profile(
                            target_count,
                            protocol_count,
                            config.reachable_probability,
                            config.unique_rank_count,
                            config.seed
                            + repetition
                            + 100 * target_count
                            + 10 * protocol_count,
                        )
                        pool = build_constraint_pool_from_dissimilarity(
                            pairwise_response_dissimilarity(profile, protocol),
                            max_constraints=max_constraints,
                            min_margin=min_margin,
                            seed=config.seed + repetition,
                        )
                        rows.append(
                            {
                                "target_count_setting": float(target_count),
                                "protocol_count": float(protocol_count),
                                "max_constraints": float(max_constraints),
                                "margin_threshold": float(min_margin),
                                "repetition": float(repetition),
                                **constraint_pool_summary(pool, target_count),
                            }
                        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write coverage CSV."""

    path = output_dir / "data" / "response_constraint_pool_coverage.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_by_max_constraints(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    labels = sorted({str(int(float(row["max_constraints"]))) for row in rows})
    values = [
        float(
            np.nanmean(
                [
                    float(row[key])
                    for row in rows
                    if str(int(float(row["max_constraints"]))) == label
                ]
            )
        )
        for label in labels
    ]
    return labels, values


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save coverage figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "touched_target_fraction",
            "response_constraint_target_coverage.png",
            "touched target fraction",
        ),
        (
            "touched_pair_node_fraction",
            "response_constraint_pair_node_coverage.png",
            "touched pair-node fraction",
        ),
    ]:
        labels, values = _mean_by_max_constraints(rows, key)
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(labels, values)
        ax.set_xlabel("max constraints")
        ax.set_ylabel(ylabel)
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
    print(f"Wrote response constraint pool coverage: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
