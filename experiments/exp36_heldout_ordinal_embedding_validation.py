"""Held-out ordinal embedding validation against null baselines."""

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
    flip_constraint_fraction,
    random_quadruplet_constraints,
    representation_generalization_report,
    shuffle_constraint_sides,
    split_constraints,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for held-out ordinal embedding validation."""

    true_dim: int = 2
    n_points_values: tuple[int, ...] = (30, 50)
    constraint_counts: tuple[int, ...] = (1000, 3000, 8000)
    repetitions: int = 5
    train_fraction: float = 0.7
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

    parser = argparse.ArgumentParser(
        description="Held-out validation for ordinal embeddings."
    )
    parser.add_argument("--true-dim", type=int, default=2)
    parser.add_argument("--n-points-values", nargs="+", default=["30", "50"])
    parser.add_argument(
        "--constraint-counts",
        nargs="+",
        default=["1000", "3000", "8000"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--train-fraction", type=float, default=0.7)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        true_dim=args.true_dim,
        n_points_values=_parse_int_values(args.n_points_values, "n_points"),
        constraint_counts=_parse_int_values(args.constraint_counts, "constraints"),
        repetitions=args.repetitions,
        train_fraction=args.train_fraction,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def _model_constraints(
    structured: np.ndarray,
    n_points: int,
    model: str,
    seed: int,
) -> np.ndarray:
    if model == "structured":
        return structured
    if model == "shuffled":
        return shuffle_constraint_sides(structured, seed=seed)
    if model == "random":
        return random_quadruplet_constraints(n_points, structured.shape[0], seed=seed)
    if model == "noisy_0.05":
        return flip_constraint_fraction(structured, 0.05, seed=seed)
    if model == "noisy_0.10":
        return flip_constraint_fraction(structured, 0.10, seed=seed)
    raise ValueError(f"unknown model: {model}")


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run held-out validation against structured and null constraints."""

    rng = np.random.default_rng(config.seed)
    rows: list[dict[str, float | str]] = []
    model_names = ("structured", "shuffled", "random", "noisy_0.05", "noisy_0.10")
    for n_points in config.n_points_values:
        for constraint_count in config.constraint_counts:
            for repetition in range(config.repetitions):
                points = rng.normal(size=(n_points, config.true_dim))
                structured = sample_quadruplet_constraints_from_coords(
                    points,
                    constraint_count,
                    seed=config.seed + 100_000 * n_points + 100 * repetition,
                )
                for model_index, model in enumerate(model_names):
                    constraints = _model_constraints(
                        structured,
                        n_points,
                        model,
                        seed=config.seed + 1_000_000 * model_index + repetition,
                    )
                    train, test = split_constraints(
                        constraints,
                        config.train_fraction,
                        seed=config.seed + 10_000 * model_index + repetition,
                    )
                    embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
                        n_points,
                        config.true_dim,
                        train,
                        steps=config.steps,
                        restarts=config.restarts,
                        learning_rate=config.learning_rate,
                        seed=config.seed + 500_000 * model_index + repetition,
                        batch_size=min(2048, max(1, train.shape[0])),
                    )
                    report = representation_generalization_report(
                        embedding,
                        train,
                        test,
                    )
                    rows.append(
                        {
                            "model": model,
                            "true_dim": float(config.true_dim),
                            "n_points": float(n_points),
                            "constraint_count": float(constraint_count),
                            "repetition": float(repetition),
                            "final_train_loss": diagnostics["final_loss"],
                            **report,
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write held-out validation rows."""

    output_path = output_dir / "data" / "heldout_ordinal_embedding_validation.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_model(rows: list[dict[str, float | str]], key: str) -> dict[str, float]:
    result: dict[str, float] = {}
    for model in sorted({str(row["model"]) for row in rows}):
        values = [float(row[key]) for row in rows if row["model"] == model]
        result[model] = float(np.nanmean(values))
    return result


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save held-out validation figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = (
        figure_dir / "heldout_violation_by_model.png",
        figure_dir / "heldout_generalization_gap.png",
    )
    for path, key, ylabel in (
        (paths[0], "test_violation_rate", "Held-out violation rate"),
        (paths[1], "generalization_gap", "Generalization gap"),
    ):
        means = _mean_by_model(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        labels = list(means)
        ax.bar(labels, [means[label] for label in labels])
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel + " by constraint model")
        ax.tick_params(axis="x", rotation=25)
        ax.grid(True, axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote held-out ordinal embedding validation data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
