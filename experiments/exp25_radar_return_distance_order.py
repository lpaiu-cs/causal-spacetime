"""Radar-return distance order from causal order and observer tick order."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer import (
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
)
from causal_spacetime_lab.ordinal import order_agreement_rate, order_inversion_rate
from causal_spacetime_lab.radar_order import radar_tick_brackets_from_order

DEFAULT_TICK_VALUES = (16, 32, 64, 128)
DEFAULT_TARGET_COUNTS = (50, 100, 200)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for radar-return order validation."""

    T: float = 2.0
    tick_values: tuple[int, ...] = DEFAULT_TICK_VALUES
    target_counts: tuple[int, ...] = DEFAULT_TARGET_COUNTS
    repetitions: int = 5
    seed: int = 0
    emission_tick_fraction: float = 0.25
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
        description="Reconstruct radar-return distance order from tick order."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--tick-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_TICK_VALUES],
    )
    parser.add_argument(
        "--target-counts",
        nargs="+",
        default=[str(value) for value in DEFAULT_TARGET_COUNTS],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--emission-tick-fraction", type=float, default=0.25)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        tick_values=_parse_int_values(args.tick_values, "tick"),
        target_counts=_parse_int_values(args.target_counts, "target"),
        repetitions=args.repetitions,
        seed=args.seed,
        emission_tick_fraction=args.emission_tick_fraction,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if any(value < 2 for value in config.tick_values):
        raise ValueError("all tick_values must be at least 2")
    if any(value <= 0 for value in config.target_counts):
        raise ValueError("all target_counts must be positive")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if not 0.0 <= config.emission_tick_fraction < 1.0:
        raise ValueError("emission_tick_fraction must satisfy 0 <= f < 1")


def _target_events(
    clock_times: np.ndarray,
    tick_count: int,
    target_count: int,
    emission_tick_fraction: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    emission_index = min(
        tick_count - 2,
        max(0, int(np.floor(emission_tick_fraction * (tick_count - 1)))),
    )
    tau_emit = float(clock_times[emission_index])
    max_radius = 0.5 * (float(clock_times[-1]) - tau_emit)
    if max_radius <= 0.0:
        raise ValueError("emission tick leaves no return range")
    rng = np.random.default_rng(seed)
    unique_count = (target_count + 1) // 2
    radii = np.sort(rng.uniform(0.03 * max_radius, 0.95 * max_radius, unique_count))
    signs = np.tile(np.asarray([1.0, -1.0]), unique_count)[:target_count]
    repeated_radii = np.repeat(radii, 2)[:target_count]
    events = np.column_stack(
        (
            tau_emit + repeated_radii,
            signs * repeated_radii,
        )
    )
    return events.astype(np.float64, copy=False), repeated_radii, emission_index


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run radar-return order validation and return target-level rows."""

    _validate_config(config)
    rows: list[dict[str, float]] = []
    for tick_index, tick_count in enumerate(config.tick_values):
        observer_events, clock_times = make_stationary_observer_chain_1p1(
            config.T,
            tick_count,
        )
        observer_indices = observer_chain_indices(0, tick_count)
        for target_index, target_count in enumerate(config.target_counts):
            for repetition in range(config.repetitions):
                targets, true_distances, emission_index = _target_events(
                    clock_times,
                    tick_count,
                    target_count,
                    config.emission_tick_fraction,
                    config.seed
                    + 1_000_000 * tick_index
                    + 10_000 * target_index
                    + repetition,
                )
                combined = np.vstack((observer_events, targets))
                target_indices = np.arange(
                    tick_count,
                    tick_count + target_count,
                    dtype=int,
                )
                causal_matrix = causal_matrix_1p1(combined)
                _, successor, accessible = radar_tick_brackets_from_order(
                    causal_matrix,
                    observer_indices,
                    target_indices,
                )
                finite = accessible & (successor >= 0)
                inversion = order_inversion_rate(
                    true_distances[finite],
                    successor[finite].astype(float),
                    ignore_ties=True,
                )
                agreement = order_agreement_rate(
                    true_distances[finite],
                    successor[finite].astype(float),
                    ignore_ties=True,
                )
                for row, target in enumerate(target_indices):
                    rows.append(
                        {
                            "tick_count": float(tick_count),
                            "target_count": float(target_count),
                            "repetition": float(repetition),
                            "target_index": float(target),
                            "emission_tick": float(emission_index),
                            "true_distance": float(true_distances[row]),
                            "successor_tick_position": float(successor[row]),
                            "accessible": float(accessible[row]),
                            "order_inversion_rate": float(inversion),
                            "order_agreement_rate": float(agreement),
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write radar-return order rows."""

    output_path = output_dir / "data" / "radar_return_distance_order.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_inversion_plot(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save inversion rate versus tick count."""

    output_path = output_dir / "figures" / "radar_return_order_inversion_vs_ticks.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for target_count in sorted({row["target_count"] for row in rows}):
        summary = []
        for tick_count in sorted({row["tick_count"] for row in rows}):
            subset = [
                row
                for row in rows
                if row["target_count"] == target_count
                and row["tick_count"] == tick_count
            ]
            values = {
                (row["repetition"], row["order_inversion_rate"]) for row in subset
            }
            summary.append(float(np.nanmean([value for _, value in values])))
        ax.plot(
            sorted({row["tick_count"] for row in rows}),
            summary,
            marker="o",
            label=f"targets={int(target_count)}",
        )
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Distance-order inversion/error rate")
    ax.set_title("Radar-return order preservation")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_scatter_plot(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save true distance versus successor tick position scatter."""

    output_path = output_dir / "figures" / "radar_return_order_scatter.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    max_tick = max(row["tick_count"] for row in rows)
    subset = [row for row in rows if row["tick_count"] == max_tick]
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    ax.scatter(
        [row["true_distance"] for row in subset],
        [row["successor_tick_position"] for row in subset],
        s=12,
        alpha=0.55,
    )
    ax.set_xlabel("True unsigned distance")
    ax.set_ylabel("Return tick position")
    ax.set_title("Radar-return order at highest tick count")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save radar-return order figures."""

    return save_inversion_plot(rows, output_dir), save_scatter_plot(rows, output_dir)


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote radar-return distance order data: {output_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
