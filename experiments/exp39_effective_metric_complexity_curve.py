"""Effective metric complexity curves for structured and null order data."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal_embedding import (
    fit_ordinal_embedding_gradient_descent,
    sample_quadruplet_constraints_from_coords,
)
from causal_spacetime_lab.representation_validation import (
    representation_generalization_report,
    shuffle_constraint_sides,
    split_constraints,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for effective metric complexity curves."""

    true_dim: int = 2
    candidate_dims: tuple[int, ...] = (1, 2, 3, 4, 5)
    n_points: int = 50
    num_constraints: int = 8000
    repetitions: int = 5
    seed: int = 0
    steps: int = 600
    restarts: int = 2
    learning_rate: float = 0.05
    complexity_lambda: float = 0.005
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

    parser = argparse.ArgumentParser(description="Effective metric complexity curves.")
    parser.add_argument("--true-dim", type=int, default=2)
    parser.add_argument(
        "--candidate-dims",
        nargs="+",
        default=["1", "2", "3", "4", "5"],
    )
    parser.add_argument("--n-points", type=int, default=50)
    parser.add_argument("--num-constraints", type=int, default=8000)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--complexity-lambda", type=float, default=0.005)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        true_dim=args.true_dim,
        candidate_dims=_parse_int_values(args.candidate_dims, "candidate dimension"),
        n_points=args.n_points,
        num_constraints=args.num_constraints,
        repetitions=args.repetitions,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        complexity_lambda=args.complexity_lambda,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run complexity-curve diagnostics."""

    rng = np.random.default_rng(config.seed)
    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        points = rng.normal(size=(config.n_points, config.true_dim))
        structured = sample_quadruplet_constraints_from_coords(
            points,
            config.num_constraints,
            seed=config.seed + repetition,
        )
        pools = {
            "structured": structured,
            "shuffled": shuffle_constraint_sides(
                structured,
                seed=config.seed + 10_000 + repetition,
            ),
        }
        for model, constraints in pools.items():
            train, test = split_constraints(
                constraints,
                train_fraction=0.7,
                seed=config.seed + repetition,
            )
            for candidate_dim in config.candidate_dims:
                embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
                    config.n_points,
                    candidate_dim,
                    train,
                    steps=config.steps,
                    restarts=config.restarts,
                    learning_rate=config.learning_rate,
                    seed=config.seed + 100_000 * candidate_dim + repetition,
                    batch_size=min(2048, max(1, train.shape[0])),
                )
                report = representation_generalization_report(embedding, train, test)
                rows.append(
                    {
                        "model": model,
                        "true_dim": float(config.true_dim),
                        "candidate_dim": float(candidate_dim),
                        "n_points": float(config.n_points),
                        "num_constraints": float(config.num_constraints),
                        "parameter_count": float(config.n_points * candidate_dim),
                        "repetition": float(repetition),
                        "final_train_loss": diagnostics["final_loss"],
                        **report,
                        "complexity_penalized_score": (
                            report["test_violation_rate"]
                            + config.complexity_lambda * candidate_dim
                        ),
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write complexity-curve rows."""

    output_path = output_dir / "data" / "effective_metric_complexity_curve.csv"
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
        dims = sorted({float(row["candidate_dim"]) for row in subset})
        values = []
        for dim in dims:
            dim_values = [
                float(row[key])
                for row in subset
                if float(row["candidate_dim"]) == dim
            ]
            values.append(float(np.nanmean(dim_values)))
        ax.plot(dims, values, marker="o", label=model)
    ax.set_xlabel("Candidate embedding dimension")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel + " vs dimension")
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
    """Save complexity-curve figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "test_violation_rate",
            "Held-out violation rate",
            figure_dir / "effective_metric_complexity_curve.png",
        ),
        _plot(
            rows,
            "complexity_penalized_score",
            "Complexity-penalized score",
            figure_dir / "effective_metric_penalized_score.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote effective metric complexity curve data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
