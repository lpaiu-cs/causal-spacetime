"""Myrheim-Meyer dimension reconstruction in flat Alexandrov intervals."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_minkowski
from causal_spacetime_lab.dimension import (
    estimate_myrheim_meyer_dimension,
    myrheim_meyer_relation_fraction,
    relation_fraction,
)
from causal_spacetime_lab.sprinkling import sprinkle_minkowski_causal_diamond

DEFAULT_DIMS = (2, 3, 4)
DEFAULT_N_VALUES = (300, 600, 1200, 2400)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for Myrheim-Meyer dimension reconstruction."""

    spacetime_dims: tuple[int, ...] = DEFAULT_DIMS
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    repetitions: int = 5
    T: float = 2.0
    seed: int = 0
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
        description=(
            "Estimate spacetime dimension from Myrheim-Meyer causal-order "
            "relation fractions in flat Alexandrov intervals."
        )
    )
    parser.add_argument(
        "--dims",
        nargs="+",
        default=[str(dim) for dim in DEFAULT_DIMS],
        help="Spacetime dimensions as space-separated values or comma lists.",
    )
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(n) for n in DEFAULT_N_VALUES],
        help="Event counts as space-separated values or comma lists.",
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    return ExperimentConfig(
        spacetime_dims=_parse_int_values(args.dims, "dimension"),
        n_values=_parse_int_values(args.n_values, "N"),
        repetitions=args.repetitions,
        T=args.T,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if any(dim < 2 for dim in config.spacetime_dims):
        raise ValueError("all spacetime dimensions must be at least 2")
    if any(n <= 1 for n in config.n_values):
        raise ValueError("all N values must be greater than 1")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.T <= 0:
        raise ValueError("T must be positive")


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    """Run dimension reconstruction and return detailed and summary rows."""

    _validate_config(config)
    result_rows: list[dict[str, float]] = []

    for dim_index, dimension in enumerate(config.spacetime_dims):
        expected_fraction = myrheim_meyer_relation_fraction(float(dimension))
        for n_index, n_events in enumerate(config.n_values):
            for repetition in range(config.repetitions):
                seed = config.seed + 100_000 * dim_index + 1_000 * n_index + repetition
                events = sprinkle_minkowski_causal_diamond(
                    n_events,
                    spacetime_dim=dimension,
                    T=config.T,
                    seed=seed,
                )
                causal_matrix = causal_matrix_minkowski(events)
                observed_fraction = relation_fraction(causal_matrix)
                estimated_dimension = estimate_myrheim_meyer_dimension(
                    observed_fraction
                )
                result_rows.append(
                    {
                        "spacetime_dim": float(dimension),
                        "N": float(n_events),
                        "repetition": float(repetition),
                        "relation_fraction": observed_fraction,
                        "expected_relation_fraction": expected_fraction,
                        "estimated_dim": estimated_dimension,
                        "error": estimated_dimension - dimension,
                    }
                )

    summary_rows = summarize_results(result_rows)
    return result_rows, summary_rows


def summarize_results(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    """Aggregate dimension reconstruction rows by dimension and event count."""

    summary_rows: list[dict[str, float]] = []
    keys = sorted({(int(row["spacetime_dim"]), int(row["N"])) for row in rows})
    for dimension, n_events in keys:
        group = [
            row
            for row in rows
            if int(row["spacetime_dim"]) == dimension and int(row["N"]) == n_events
        ]
        estimates = np.asarray([row["estimated_dim"] for row in group])
        errors = estimates - dimension
        relation_fractions = np.asarray([row["relation_fraction"] for row in group])
        summary_rows.append(
            {
                "spacetime_dim": float(dimension),
                "N": float(n_events),
                "repetitions": float(len(group)),
                "mean_estimated_dim": float(np.mean(estimates)),
                "std_estimated_dim": float(np.std(estimates, ddof=1))
                if len(group) > 1
                else 0.0,
                "bias": float(np.mean(errors)),
                "rmse": float(np.sqrt(np.mean(errors * errors))),
                "mean_relation_fraction": float(np.mean(relation_fractions)),
                "expected_relation_fraction": myrheim_meyer_relation_fraction(
                    float(dimension)
                ),
            }
        )
    return summary_rows


def _write_csv(rows: list[dict[str, float]], output_path: Path) -> Path:
    if not rows:
        raise RuntimeError("no rows to write")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def write_outputs(
    result_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write dimension reconstruction CSV outputs."""

    data_dir = output_dir / "data"
    result_path = data_dir / "dimension_reconstruction_results.csv"
    summary_path = data_dir / "dimension_reconstruction_summary.csv"
    return _write_csv(result_rows, result_path), _write_csv(summary_rows, summary_path)


def save_dimension_estimate_plot(
    summary_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save estimated dimension versus N."""

    figure_path = output_dir / "figures" / "dimension_estimate_vs_N.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    dimensions = sorted({int(row["spacetime_dim"]) for row in summary_rows})
    for dimension in dimensions:
        rows = [row for row in summary_rows if int(row["spacetime_dim"]) == dimension]
        n_values = np.asarray([row["N"] for row in rows])
        estimates = np.asarray([row["mean_estimated_dim"] for row in rows])
        std = np.asarray([row["std_estimated_dim"] for row in rows])
        ax.errorbar(n_values, estimates, yerr=std, marker="o", label=f"D={dimension}")
        ax.axhline(dimension, color="0.65", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Event count N")
    ax.set_ylabel("Estimated dimension")
    ax.set_title("Myrheim-Meyer dimension estimate versus N")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_relation_fraction_plot(
    result_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> Path:
    """Save relation fraction curve with observed points."""

    figure_path = config.output_dir / "figures" / "relation_fraction_vs_dimension.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    max_dim = max(max(config.spacetime_dims) + 1.0, 5.0)
    curve_dims = np.linspace(1.05, max_dim, 300)
    curve_values = np.asarray(
        [myrheim_meyer_relation_fraction(float(dim)) for dim in curve_dims]
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(curve_dims, curve_values, color="black", label="Myrheim-Meyer curve")
    for dimension in sorted({int(row["spacetime_dim"]) for row in result_rows}):
        rows = [row for row in result_rows if int(row["spacetime_dim"]) == dimension]
        x = np.full(len(rows), dimension, dtype=float)
        y = np.asarray([row["relation_fraction"] for row in rows])
        ax.scatter(x, y, s=16, alpha=0.45, label=f"Observed D={dimension}")
    ax.set_xlabel("Spacetime dimension D")
    ax.set_ylabel("Ordered relation fraction")
    ax.set_title("Order fraction as a dimension observable")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_dimension_error_plot(
    result_rows: list[dict[str, float]],
    output_dir: Path,
) -> Path:
    """Save estimated dimension error versus N."""

    figure_path = output_dir / "figures" / "dimension_error_vs_N.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for dimension in sorted({int(row["spacetime_dim"]) for row in result_rows}):
        rows = [row for row in result_rows if int(row["spacetime_dim"]) == dimension]
        n_values = np.asarray([row["N"] for row in rows])
        errors = np.asarray([row["error"] for row in rows])
        ax.scatter(n_values, errors, s=18, alpha=0.55, label=f"D={dimension}")
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("Event count N")
    ax.set_ylabel("estimated_dim - true_dim")
    ax.set_title("Dimension reconstruction error")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    return figure_path


def save_figures(
    result_rows: list[dict[str, float]],
    summary_rows: list[dict[str, float]],
    config: ExperimentConfig,
) -> tuple[Path, Path, Path]:
    """Save dimension reconstruction figures."""

    return (
        save_dimension_estimate_plot(summary_rows, config.output_dir),
        save_relation_fraction_plot(result_rows, config),
        save_dimension_error_plot(result_rows, config.output_dir),
    )


def main() -> None:
    config = parse_args()
    result_rows, summary_rows = run_experiment(config)
    result_path, summary_path = write_outputs(
        result_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_figures(result_rows, summary_rows, config)
    print(f"Wrote results: {result_path}")
    print(f"Wrote summary: {summary_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")
    print(
        "Dimension is estimated from causal-order statistics in known flat "
        "Alexandrov intervals; this is a controlled validation, not a proof of "
        "spacetime emergence."
    )


if __name__ == "__main__":
    main()

