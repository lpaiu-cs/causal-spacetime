"""Ordinal embedding recovery from pure distance-order constraints."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    procrustes_align,
    quadruplet_violation_rate,
    sample_quadruplet_constraints_from_coords,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for ordinal embedding recovery."""

    true_dims: tuple[int, ...] = (1, 2)
    n_points_values: tuple[int, ...] = (20, 40)
    constraint_counts: tuple[int, ...] = (1000, 3000, 8000)
    repetitions: int = 5
    seed: int = 0
    steps: int = 400
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

    parser = argparse.ArgumentParser(description="Recover embeddings from order data.")
    parser.add_argument("--true-dims", nargs="+", default=["1", "2"])
    parser.add_argument("--n-points-values", nargs="+", default=["20", "40"])
    parser.add_argument(
        "--constraint-counts",
        nargs="+",
        default=["1000", "3000", "8000"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        true_dims=_parse_int_values(args.true_dims, "dimension"),
        n_points_values=_parse_int_values(args.n_points_values, "n_points"),
        constraint_counts=_parse_int_values(args.constraint_counts, "constraints"),
        repetitions=args.repetitions,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run controlled ordinal embedding recovery."""

    rows: list[dict[str, float]] = []
    rng = np.random.default_rng(config.seed)
    for true_dim in config.true_dims:
        for n_points in config.n_points_values:
            for constraint_count in config.constraint_counts:
                for repetition in range(config.repetitions):
                    points = rng.normal(size=(n_points, true_dim))
                    constraints = sample_quadruplet_constraints_from_coords(
                        points,
                        constraint_count,
                        seed=config.seed
                        + 1_000_000 * true_dim
                        + 10_000 * n_points
                        + 100 * constraint_count
                        + repetition,
                    )
                    embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
                        n_points,
                        true_dim,
                        constraints,
                        steps=config.steps,
                        learning_rate=config.learning_rate,
                        restarts=config.restarts,
                        seed=config.seed + repetition + 17 * true_dim,
                        batch_size=min(2048, max(1, constraints.shape[0])),
                    )
                    aligned, alignment = procrustes_align(embedding, points)
                    rows.append(
                        {
                            "true_dim": float(true_dim),
                            "embedding_dim": float(true_dim),
                            "n_points": float(n_points),
                            "constraint_count": float(constraint_count),
                            "repetition": float(repetition),
                            "final_hinge_loss": diagnostics["final_loss"],
                            "constraint_violation_rate": quadruplet_violation_rate(
                                embedding,
                                constraints,
                            ),
                            "distance_order_error": embedding_distance_order_error(
                                aligned,
                                points,
                                seed=config.seed + repetition,
                            ),
                            "procrustes_rmse": alignment["rmse"],
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write ordinal embedding recovery rows."""

    output_path = output_dir / "data" / "ordinal_embedding_recovery.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by(
    rows: list[dict[str, float]],
    y_key: str,
) -> dict[tuple[float, float], float]:
    grouped: dict[tuple[float, float], list[float]] = {}
    for row in rows:
        key = (row["true_dim"], row["constraint_count"])
        grouped.setdefault(key, []).append(row[y_key])
    return {key: float(np.nanmean(values)) for key, values in grouped.items()}


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save violation and RMSE figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = (
        figure_dir / "ordinal_embedding_violation_vs_constraints.png",
        figure_dir / "ordinal_embedding_rmse_vs_constraints.png",
    )
    for path, key, label in (
        (paths[0], "constraint_violation_rate", "Constraint violation rate"),
        (paths[1], "procrustes_rmse", "Procrustes RMSE"),
    ):
        means = _mean_by(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for true_dim in sorted({row["true_dim"] for row in rows}):
            x_values = sorted(
                {constraint for dim, constraint in means if dim == true_dim}
            )
            y_values = [means[(true_dim, x)] for x in x_values]
            ax.plot(x_values, y_values, marker="o", label=f"true dim={int(true_dim)}")
        ax.set_xscale("log")
        ax.set_xlabel("Distance-order constraint count")
        ax.set_ylabel(label)
        ax.set_title(label + " vs constraints")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote ordinal embedding recovery data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
