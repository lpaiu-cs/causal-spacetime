"""Ordinal embedding from observer-derived distance-order constraints."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    append_oriented_protocol_chains,
    common_safe_tau_range_for_oriented_protocol_1p1,
    reconstruct_oriented_chart_from_order,
)
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    procrustes_align,
    quadruplet_violation_rate,
    sample_quadruplet_constraints_from_distance_values,
    squared_distance_matrix,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for observer-derived distance-order embedding."""

    T: float = 2.0
    n_events: int = 400
    tick_values: tuple[int, ...] = (32, 64, 128)
    constraint_counts: tuple[int, ...] = (1000, 3000)
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
    steps: int = 1000
    restarts: int = 3
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

    parser = argparse.ArgumentParser(
        description="Fit 1D embeddings from observer-derived distance order."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=400)
    parser.add_argument("--tick-values", nargs="+", default=["32", "64", "128"])
    parser.add_argument("--constraint-counts", nargs="+", default=["1000", "3000"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--restarts", type=int, default=3)
    parser.add_argument("--learning-rate", type=float, default=0.1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_values=_parse_int_values(args.tick_values, "tick"),
        constraint_counts=_parse_int_values(args.constraint_counts, "constraint"),
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def _reconstruct_positions(
    support_events: np.ndarray,
    tick_count: int,
    config: ExperimentConfig,
) -> tuple[np.ndarray, np.ndarray]:
    spec = ObserverProtocolSpec(
        name="stationary_oriented",
        beta=0.0,
        origin_lab_time=0.0,
        origin_lab_position=0.0,
        beacon_separation=config.beacon_separation,
    )
    tau_range = common_safe_tau_range_for_oriented_protocol_1p1(config.T, spec)
    chain_events: list[np.ndarray] = []
    indices, _ = append_oriented_protocol_chains(
        chain_events,
        spec,
        tick_count,
        tau_range,
        support_events.shape[0],
    )
    combined = np.vstack((support_events, *chain_events))
    chart = reconstruct_oriented_chart_from_order(
        causal_matrix_1p1(combined),
        np.arange(support_events.shape[0], dtype=int),
        indices.primary,
        indices.beacon,
        indices.clock_times,
        config.beacon_separation,
        spec.name,
    )
    accessible = chart.accessible
    return support_events[accessible, 1], chart.reconstructed_coords[accessible, 1]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run observer-derived distance-order embedding diagnostics."""

    rows: list[dict[str, float]] = []
    for repetition in range(config.repetitions):
        support_events = sprinkle_1p1_causal_diamond(
            config.n_events,
            T=config.T,
            seed=config.seed + repetition,
        )
        for tick_count in config.tick_values:
            true_x, reconstructed_x = _reconstruct_positions(
                support_events,
                tick_count,
                config,
            )
            accessible_count = true_x.size
            if accessible_count < 4:
                continue
            distance_matrix = squared_distance_matrix(reconstructed_x[:, None])
            for constraint_count in config.constraint_counts:
                constraints = sample_quadruplet_constraints_from_distance_values(
                    distance_matrix,
                    constraint_count,
                    seed=config.seed + 10_000 * tick_count + repetition,
                )
                embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
                    accessible_count,
                    1,
                    constraints,
                    steps=config.steps,
                    learning_rate=config.learning_rate,
                    restarts=config.restarts,
                    seed=config.seed + 100_000 * tick_count + repetition,
                    batch_size=min(2048, max(1, constraints.shape[0])),
                )
                aligned, alignment = procrustes_align(embedding, true_x[:, None])
                rows.append(
                    {
                        "tick_count": float(tick_count),
                        "constraint_count": float(constraint_count),
                        "repetition": float(repetition),
                        "accessible_count": float(accessible_count),
                        "accessible_fraction": float(
                            accessible_count / config.n_events
                        ),
                        "final_hinge_loss": diagnostics["final_loss"],
                        "violation_rate": quadruplet_violation_rate(
                            embedding,
                            constraints,
                        ),
                        "distance_order_error_against_true_x": (
                            embedding_distance_order_error(
                                aligned,
                                true_x[:, None],
                                seed=config.seed + repetition,
                            )
                        ),
                        "alignment_rmse": alignment["rmse"],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write observer-derived embedding rows."""

    output_path = output_dir / "data" / "observer_distance_order_embedding.csv"
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
    for constraint_count in sorted({row["constraint_count"] for row in rows}):
        subset = [row for row in rows if row["constraint_count"] == constraint_count]
        ticks = sorted({row["tick_count"] for row in subset})
        values = [
            float(np.nanmean([row[key] for row in subset if row["tick_count"] == tick]))
            for tick in ticks
        ]
        ax.plot(ticks, values, marker="o", label=f"constraints={int(constraint_count)}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel(key.replace("_", " "))
    ax.set_title(key.replace("_", " ").title())
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save observer-derived embedding figures."""

    return (
        _plot(
            rows,
            output_dir,
            "violation_rate",
            "observer_distance_order_embedding_violation_vs_ticks.png",
        ),
        _plot(
            rows,
            output_dir,
            "alignment_rmse",
            "observer_distance_order_embedding_rmse_vs_ticks.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote observer distance-order embedding data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
