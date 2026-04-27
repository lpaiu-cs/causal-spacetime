"""Ordinal stress curves for effective embedding dimension diagnostics."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal_embedding import (
    ordinal_embedding_dimension_curve,
    sample_quadruplet_constraints_from_coords,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for ordinal dimension selection diagnostics."""

    true_dims: tuple[int, ...] = (1, 2, 3)
    candidate_dims: tuple[int, ...] = (1, 2, 3, 4)
    n_points: int = 40
    num_constraints: int = 6000
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

    parser = argparse.ArgumentParser(description="Fit ordinal embeddings by dimension.")
    parser.add_argument("--true-dims", nargs="+", default=["1", "2", "3"])
    parser.add_argument("--candidate-dims", nargs="+", default=["1", "2", "3", "4"])
    parser.add_argument("--n-points", type=int, default=40)
    parser.add_argument("--num-constraints", type=int, default=6000)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        true_dims=_parse_int_values(args.true_dims, "true_dim"),
        candidate_dims=_parse_int_values(args.candidate_dims, "candidate_dim"),
        n_points=args.n_points,
        num_constraints=args.num_constraints,
        repetitions=args.repetitions,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run ordinal dimension curve diagnostics."""

    rows: list[dict[str, float]] = []
    rng = np.random.default_rng(config.seed)
    for true_dim in config.true_dims:
        for repetition in range(config.repetitions):
            points = rng.normal(size=(config.n_points, true_dim))
            constraints = sample_quadruplet_constraints_from_coords(
                points,
                config.num_constraints,
                seed=config.seed + 10_000 * true_dim + repetition,
            )
            curve = ordinal_embedding_dimension_curve(
                config.n_points,
                constraints,
                list(config.candidate_dims),
                steps=config.steps,
                restarts=config.restarts,
                learning_rate=config.learning_rate,
                seed=config.seed + 100_000 * true_dim + repetition,
                batch_size=min(2048, max(1, constraints.shape[0])),
            )
            for row in curve:
                rows.append(
                    {
                        "true_dim": float(true_dim),
                        "candidate_dim": row["embedding_dim"],
                        "n_points": float(config.n_points),
                        "num_constraints": float(config.num_constraints),
                        "repetition": float(repetition),
                        "final_hinge_loss": row["final_loss"],
                        "violation_rate": row["violation_rate"],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write dimension selection rows."""

    output_path = output_dir / "data" / "ordinal_dimension_selection.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot_metric(
    rows: list[dict[str, float]],
    output_dir: Path,
    key: str,
    name: str,
) -> Path:
    output_path = output_dir / "figures" / name
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for true_dim in sorted({row["true_dim"] for row in rows}):
        subset = [row for row in rows if row["true_dim"] == true_dim]
        candidate_dims = sorted({row["candidate_dim"] for row in subset})
        means = [
            float(
                np.nanmean(
                    [row[key] for row in subset if row["candidate_dim"] == dim]
                )
            )
            for dim in candidate_dims
        ]
        ax.plot(candidate_dims, means, marker="o", label=f"true dim={int(true_dim)}")
        ax.axvline(true_dim, color="gray", linewidth=0.7, alpha=0.25)
    ax.set_xlabel("Candidate embedding dimension")
    ax.set_ylabel(key.replace("_", " "))
    ax.set_title(key.replace("_", " ").title())
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save dimension stress and violation curves."""

    return (
        _plot_metric(
            rows,
            output_dir,
            "final_hinge_loss",
            "ordinal_dimension_stress_curve.png",
        ),
        _plot_metric(
            rows,
            output_dir,
            "violation_rate",
            "ordinal_dimension_violation_curve.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote ordinal dimension selection data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
