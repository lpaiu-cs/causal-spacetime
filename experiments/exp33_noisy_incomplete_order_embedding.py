"""Robustness of ordinal embeddings under noisy incomplete constraints."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal_embedding import (
    fit_ordinal_embedding_gradient_descent,
    procrustes_align,
    quadruplet_violation_rate,
    sample_quadruplet_constraints_from_coords,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for noisy incomplete ordinal embedding."""

    true_dim: int = 2
    n_points: int = 40
    constraint_counts: tuple[int, ...] = (1000, 3000, 8000)
    flip_probabilities: tuple[float, ...] = (0.0, 0.02, 0.05, 0.10)
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


def _parse_float_values(values: list[str], name: str) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Noisy ordinal embedding diagnostic.")
    parser.add_argument("--true-dim", type=int, default=2)
    parser.add_argument("--n-points", type=int, default=40)
    parser.add_argument(
        "--constraint-counts",
        nargs="+",
        default=["1000", "3000", "8000"],
    )
    parser.add_argument(
        "--flip-probabilities",
        nargs="+",
        default=["0.0", "0.02", "0.05", "0.10"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        true_dim=args.true_dim,
        n_points=args.n_points,
        constraint_counts=_parse_int_values(args.constraint_counts, "constraint"),
        flip_probabilities=_parse_float_values(args.flip_probabilities, "flip"),
        repetitions=args.repetitions,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def _flip_constraints(
    constraints: np.ndarray,
    probability: float,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    noisy = constraints.copy()
    mask = rng.uniform(0.0, 1.0, constraints.shape[0]) < probability
    left = noisy[mask, 0:2].copy()
    noisy[mask, 0:2] = noisy[mask, 2:4]
    noisy[mask, 2:4] = left
    return noisy


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run noisy incomplete embedding diagnostics."""

    rows: list[dict[str, float]] = []
    rng = np.random.default_rng(config.seed)
    for repetition in range(config.repetitions):
        points = rng.normal(size=(config.n_points, config.true_dim))
        for constraint_count in config.constraint_counts:
            clean_constraints = sample_quadruplet_constraints_from_coords(
                points,
                constraint_count,
                seed=config.seed + 10_000 * constraint_count + repetition,
            )
            validation_constraints = sample_quadruplet_constraints_from_coords(
                points,
                min(4000, max(1000, constraint_count)),
                seed=config.seed + 1_000_000 + 10_000 * constraint_count + repetition,
            )
            for flip_probability in config.flip_probabilities:
                noisy_constraints = _flip_constraints(
                    clean_constraints,
                    flip_probability,
                    config.seed + 7_000_000 + repetition,
                )
                embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
                    config.n_points,
                    config.true_dim,
                    noisy_constraints,
                    steps=config.steps,
                    learning_rate=config.learning_rate,
                    restarts=config.restarts,
                    seed=config.seed + repetition,
                    batch_size=min(2048, max(1, noisy_constraints.shape[0])),
                )
                _, alignment = procrustes_align(embedding, points)
                rows.append(
                    {
                        "true_dim": float(config.true_dim),
                        "n_points": float(config.n_points),
                        "constraint_count": float(constraint_count),
                        "flip_probability": float(flip_probability),
                        "repetition": float(repetition),
                        "final_hinge_loss": diagnostics["final_loss"],
                        "noisy_violation_rate": quadruplet_violation_rate(
                            embedding,
                            noisy_constraints,
                        ),
                        "clean_validation_violation_rate": quadruplet_violation_rate(
                            embedding,
                            validation_constraints,
                        ),
                        "procrustes_rmse": alignment["rmse"],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write noisy incomplete embedding rows."""

    output_path = output_dir / "data" / "noisy_incomplete_order_embedding.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(rows: list[dict[str, float]], output_dir: Path, key: str, name: str) -> Path:
    output_path = output_dir / "figures" / name
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for flip_probability in sorted({row["flip_probability"] for row in rows}):
        subset = [row for row in rows if row["flip_probability"] == flip_probability]
        x_values = sorted({row["constraint_count"] for row in subset})
        y_values = [
            float(
                np.nanmean(
                    [row[key] for row in subset if row["constraint_count"] == x]
                )
            )
            for x in x_values
        ]
        ax.plot(x_values, y_values, marker="o", label=f"flip={flip_probability:g}")
    ax.set_xscale("log")
    ax.set_xlabel("Constraint count")
    ax.set_ylabel(key.replace("_", " "))
    ax.set_title(key.replace("_", " ").title())
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save noisy embedding figures."""

    return (
        _plot(
            rows,
            output_dir,
            "clean_validation_violation_rate",
            "noisy_order_embedding_violation.png",
        ),
        _plot(rows, output_dir, "procrustes_rmse", "noisy_order_embedding_rmse.png"),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote noisy incomplete order embedding data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
