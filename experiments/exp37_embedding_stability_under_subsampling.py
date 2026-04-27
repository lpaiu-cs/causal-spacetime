"""Embedding stability under independent order-constraint subsampling."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.embedding_stability import (
    fit_embeddings_on_constraint_subsets,
    pairwise_order_stability,
    pairwise_procrustes_stability,
)
from causal_spacetime_lab.ordinal_embedding import (
    sample_quadruplet_constraints_from_coords,
)
from causal_spacetime_lab.representation_validation import shuffle_constraint_sides

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for subsampling stability diagnostics."""

    true_dim: int = 2
    n_points: int = 40
    total_constraints: int = 10000
    subset_sizes: tuple[int, ...] = (500, 1000, 3000, 6000)
    num_subsets: int = 5
    repetitions: int = 5
    seed: int = 0
    steps: int = 600
    restarts: int = 2
    learning_rate: float = 0.05
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str], name: str) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Embedding stability diagnostics.")
    parser.add_argument("--true-dim", type=int, default=2)
    parser.add_argument("--n-points", type=int, default=40)
    parser.add_argument("--total-constraints", type=int, default=10000)
    parser.add_argument(
        "--subset-sizes",
        nargs="+",
        default=["500", "1000", "3000", "6000"],
    )
    parser.add_argument("--num-subsets", type=int, default=5)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        true_dim=args.true_dim,
        n_points=args.n_points,
        total_constraints=args.total_constraints,
        subset_sizes=_parse_int_values(args.subset_sizes, "subset size"),
        num_subsets=args.num_subsets,
        repetitions=args.repetitions,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run stability diagnostics for structured and shuffled constraints."""

    rng = np.random.default_rng(config.seed)
    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        points = rng.normal(size=(config.n_points, config.true_dim))
        structured = sample_quadruplet_constraints_from_coords(
            points,
            config.total_constraints,
            seed=config.seed + repetition,
        )
        pools = {
            "structured": structured,
            "shuffled": shuffle_constraint_sides(
                structured,
                seed=config.seed + 10_000 + repetition,
            ),
        }
        for model, pool in pools.items():
            for subset_size in config.subset_sizes:
                embeddings = fit_embeddings_on_constraint_subsets(
                    config.n_points,
                    config.true_dim,
                    pool,
                    num_subsets=config.num_subsets,
                    subset_size=subset_size,
                    seed=config.seed + 100_000 * subset_size + repetition,
                    steps=config.steps,
                    restarts=config.restarts,
                    learning_rate=config.learning_rate,
                )
                procrustes = pairwise_procrustes_stability(embeddings)
                order = pairwise_order_stability(
                    embeddings,
                    seed=config.seed + subset_size + repetition,
                )
                rows.append(
                    {
                        "model": model,
                        "true_dim": float(config.true_dim),
                        "n_points": float(config.n_points),
                        "total_constraints": float(config.total_constraints),
                        "subset_size": float(subset_size),
                        "num_subsets": float(config.num_subsets),
                        "repetition": float(repetition),
                        **procrustes,
                        **order,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write subsampling stability rows."""

    output_path = output_dir / "data" / "embedding_stability_under_subsampling.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(
    rows: list[dict[str, float | str]],
    key: str,
    ylabel: str,
    output_path: Path,
) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for model in sorted({str(row["model"]) for row in rows}):
        subset = [row for row in rows if row["model"] == model]
        x_values = sorted({float(row["subset_size"]) for row in subset})
        y_values = []
        for x in x_values:
            values = [
                float(row[key])
                for row in subset
                if float(row["subset_size"]) == x
            ]
            y_values.append(float(np.nanmean(values)))
        ax.plot(x_values, y_values, marker="o", label=model)
    ax.set_xscale("log")
    ax.set_xlabel("Constraint subset size")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel + " vs subset size")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save subsampling stability figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "mean_procrustes_rmse",
            "Mean pairwise Procrustes RMSE",
            figure_dir / "embedding_procrustes_stability_vs_constraints.png",
        ),
        _plot(
            rows,
            "mean_order_disagreement",
            "Mean pairwise order disagreement",
            figure_dir / "embedding_order_stability_vs_constraints.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote embedding stability data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
