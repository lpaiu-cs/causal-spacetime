"""Slice-local embedding validation without cross-slice transport."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.observer_slice_pipeline import (
    reconstruct_stationary_oriented_slices_1p1,
)
from causal_spacetime_lab.ordinal import pair_distance_values_1d
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    procrustes_align,
)
from causal_spacetime_lab.slice_local_embedding import (
    fit_slice_local_ordinal_embeddings,
)
from causal_spacetime_lab.slice_local_validation import (
    summarize_slice_local_validation,
    validate_slice_local_embeddings_against_true_positions,
)
from causal_spacetime_lab.sliced_distance_order import (
    pair_slice_labels,
    quadruplet_constraints_from_sliced_pair_distances,
)
from causal_spacetime_lab.spatial_slices import same_slice_unordered_pairs
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for slice-local embedding validation."""

    T: float = 2.0
    n_values: tuple[int, ...] = (500, 1000)
    tick_values: tuple[int, ...] = (64, 128)
    bin_width_values: tuple[int, ...] = (2, 4)
    constraint_count: int = 5000
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
    steps: int = 600
    restarts: int = 2
    learning_rate: float = 0.1
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

    parser = argparse.ArgumentParser(description="Slice-local embedding validation.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N-values", nargs="+", default=["500", "1000"])
    parser.add_argument("--tick-values", nargs="+", default=["64", "128"])
    parser.add_argument("--bin-width-values", nargs="+", default=["2", "4"])
    parser.add_argument("--constraint-count", type=int, default=5000)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.N_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        bin_width_values=_parse_int_values(args.bin_width_values, "bin width"),
        constraint_count=args.constraint_count,
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def _constraints_for_result(
    result: dict[str, np.ndarray],
    constraint_count: int,
    seed: int,
) -> np.ndarray:
    pairs = same_slice_unordered_pairs(
        result["slice_labels"],
        max_pairs_per_slice=500,
        seed=seed,
    )
    if pairs.shape[0] < 2:
        return np.empty((0, 4), dtype=int)
    distances = pair_distance_values_1d(result["reconstructed_X"], pairs)
    pair_labels = pair_slice_labels(pairs, result["slice_labels"])
    return quadruplet_constraints_from_sliced_pair_distances(
        pairs,
        distances,
        pair_labels,
        constraint_count,
        seed=seed + 1,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run slice-local embedding validation."""

    rows: list[dict[str, float]] = []
    for n_events in config.n_values:
        for tick_count in config.tick_values:
            for bin_width in config.bin_width_values:
                for repetition in range(config.repetitions):
                    support = sprinkle_1p1_causal_diamond(
                        n_events,
                        T=config.T,
                        seed=config.seed + 100_000 * n_events + repetition,
                    )
                    result = reconstruct_stationary_oriented_slices_1p1(
                        support,
                        config.T,
                        tick_count,
                        config.beacon_separation,
                        bin_width=bin_width,
                    )
                    constraints = _constraints_for_result(
                        result,
                        config.constraint_count,
                        seed=config.seed + 10_000 * tick_count + repetition,
                    )
                    local_embeddings = fit_slice_local_ordinal_embeddings(
                        n_events,
                        result["slice_labels"],
                        constraints,
                        steps=config.steps,
                        restarts=config.restarts,
                        learning_rate=config.learning_rate,
                        seed=config.seed + repetition,
                    )
                    local_rows = validate_slice_local_embeddings_against_true_positions(
                        local_embeddings,
                        support[:, 1],
                    )
                    summary = summarize_slice_local_validation(local_rows)
                    if constraints.shape[0] >= 2:
                        global_embedding, _ = fit_ordinal_embedding_gradient_descent(
                            n_events,
                            1,
                            constraints,
                            steps=config.steps,
                            restarts=config.restarts,
                            learning_rate=config.learning_rate,
                            seed=config.seed + 50_000 + repetition,
                            batch_size=min(2048, constraints.shape[0]),
                        )
                        aligned, diagnostics = procrustes_align(
                            global_embedding,
                            support[:, 1, None],
                        )
                        global_rmse = diagnostics["rmse"]
                        global_order_error = embedding_distance_order_error(
                            aligned,
                            support[:, 1, None],
                            seed=config.seed + repetition,
                        )
                    else:
                        global_rmse = float("nan")
                        global_order_error = float("nan")
                    rows.append(
                        {
                            "N": float(n_events),
                            "tick_count": float(tick_count),
                            "bin_width": float(bin_width),
                            "repetition": float(repetition),
                            "constraint_count": float(constraints.shape[0]),
                            "global_no_transport_rmse": global_rmse,
                            "global_no_transport_order_error": global_order_error,
                            **summary,
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write slice-local validation rows."""

    output_path = output_dir / "data" / "slice_local_embedding_validation.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(rows: list[dict[str, float]], key: str, ylabel: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for tick_count in sorted({row["tick_count"] for row in rows}):
        subset = [row for row in rows if row["tick_count"] == tick_count]
        widths = sorted({row["bin_width"] for row in subset})
        values = [
            float(np.nanmean([row[key] for row in subset if row["bin_width"] == width]))
            for width in widths
        ]
        ax.plot(widths, values, marker="o", label=f"ticks={int(tick_count)}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Radar-time bin width")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save slice-local validation figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "mean_local_order_error",
            "Mean local distance-order error",
            figure_dir / "slice_local_order_error_vs_ticks.png",
        ),
        _plot(
            rows,
            "mean_local_rmse",
            "Mean local RMSE",
            figure_dir / "slice_local_rmse_vs_ticks.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote slice-local embedding validation data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
