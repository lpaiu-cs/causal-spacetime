"""Single-reference response-order underdetermination sweep."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_response_underdetermination import (
    single_reference_underdetermination_report,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for single-reference underdetermination diagnostics."""

    target_counts: tuple[int, ...] = (8, 16, 32)
    unique_rank_counts: tuple[int, ...] = (3, 5, 8)
    layout_count: int = 50
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

    parser = argparse.ArgumentParser(
        description="Single-reference response-order underdetermination."
    )
    parser.add_argument("--target-counts", nargs="+", default=["8", "16", "32"])
    parser.add_argument("--unique-rank-counts", nargs="+", default=["3", "5", "8"])
    parser.add_argument("--layout-count", type=int, default=50)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        target_counts=_parse_int_values(args.target_counts),
        unique_rank_counts=_parse_int_values(args.unique_rank_counts),
        layout_count=args.layout_count,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _synthetic_delay_ranks(
    target_count: int,
    unique_rank_count: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ranks = np.resize(np.arange(1, unique_rank_count + 1), target_count)
    rng.shuffle(ranks)
    return ranks.astype(int)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run the underdetermination sweep."""

    rows: list[dict[str, float]] = []
    for target_count in config.target_counts:
        for unique_rank_count in config.unique_rank_counts:
            if unique_rank_count > target_count:
                continue
            for repetition in range(config.repetitions):
                seed = config.seed + 1000 * target_count + 100 * unique_rank_count
                seed += repetition
                delays = _synthetic_delay_ranks(target_count, unique_rank_count, seed)
                reachable = np.ones(target_count, dtype=bool)
                report = single_reference_underdetermination_report(
                    delays,
                    reachable,
                    layout_count=config.layout_count,
                    seed=seed + 10,
                )
                rows.append(
                    {
                        "target_count": float(target_count),
                        "unique_rank_count": float(unique_rank_count),
                        "repetition": float(repetition),
                        **report,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = (
        output_dir
        / "data"
        / "single_reference_response_order_underdetermination.csv"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_rows(
    rows: list[dict[str, float]],
    group_key: str,
    value_key: str,
) -> tuple[list[float], list[float]]:
    groups = sorted({float(row[group_key]) for row in rows})
    means = [
        float(
            np.mean(
                [
                    row[value_key]
                    for row in rows
                    if float(row[group_key]) == group
                    and np.isfinite(row[value_key])
                ]
            )
        )
        for group in groups
    ]
    return groups, means


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for group_key, filename, xlabel in [
        (
            "target_count",
            "single_reference_pair_distance_disagreement.png",
            "target count",
        ),
        (
            "unique_rank_count",
            "single_reference_unique_rank_effect.png",
            "unique response ranks",
        ),
    ]:
        xs, ys = _mean_rows(rows, group_key, "mean_pair_distance_order_disagreement")
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("mean pair-distance-order disagreement")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
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
    print(f"Wrote single-reference underdetermination: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
