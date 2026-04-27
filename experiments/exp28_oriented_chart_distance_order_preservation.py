"""Oriented chart distance-order preservation diagnostics."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.metric_representation import unordered_pair_indices
from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    append_oriented_protocol_chains,
    common_safe_tau_range_for_oriented_protocol_1p1,
    reconstruct_oriented_chart_from_order,
)
from causal_spacetime_lab.ordinal import pair_distance_order_inversion_rate
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (300, 600, 1200)
DEFAULT_TICK_VALUES = (32, 64, 128)
DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for oriented chart distance-order validation."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    tick_values: tuple[int, ...] = DEFAULT_TICK_VALUES
    repetitions: int = 5
    seed: int = 0
    beacon_separation: float = 0.15
    pair_count: int = 500
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
        description="Check oriented radar chart distance-order preservation."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_VALUES],
    )
    parser.add_argument(
        "--tick-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_TICK_VALUES],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--pair-count", type=int, default=500)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        repetitions=args.repetitions,
        seed=args.seed,
        beacon_separation=args.beacon_separation,
        pair_count=args.pair_count,
        output_dir=args.output_dir,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0.0:
        raise ValueError("T must be positive")
    if any(value <= 0 for value in config.n_values):
        raise ValueError("all n_values must be positive")
    if any(value < 2 for value in config.tick_values):
        raise ValueError("all tick_values must be at least 2")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.beacon_separation <= 0.0:
        raise ValueError("beacon_separation must be positive")
    if config.pair_count <= 0:
        raise ValueError("pair_count must be positive")


def _sample_pairs(
    accessible_count: int,
    pair_count: int,
    seed: int,
) -> np.ndarray:
    all_pairs = unordered_pair_indices(accessible_count)
    if all_pairs.shape[0] <= pair_count:
        return all_pairs
    rng = np.random.default_rng(seed)
    selected = rng.choice(all_pairs.shape[0], size=pair_count, replace=False)
    return all_pairs[selected]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run oriented chart distance-order diagnostics."""

    _validate_config(config)
    rows: list[dict[str, float]] = []
    spec = ObserverProtocolSpec(
        name="stationary_oriented",
        beta=0.0,
        origin_lab_time=0.0,
        origin_lab_position=0.0,
        beacon_separation=config.beacon_separation,
    )
    for n_index, n_events in enumerate(config.n_values):
        for repetition in range(config.repetitions):
            support_events = sprinkle_1p1_causal_diamond(
                n_events,
                T=config.T,
                seed=config.seed + 10_000 * n_index + repetition,
            )
            target_indices = np.arange(n_events, dtype=int)
            for tick_index, tick_count in enumerate(config.tick_values):
                tau_range = common_safe_tau_range_for_oriented_protocol_1p1(
                    config.T,
                    spec,
                )
                chain_events: list[np.ndarray] = []
                indices, _ = append_oriented_protocol_chains(
                    chain_events,
                    spec,
                    tick_count,
                    tau_range,
                    n_events,
                )
                combined = np.vstack((support_events, *chain_events))
                causal_matrix = causal_matrix_1p1(combined)
                chart = reconstruct_oriented_chart_from_order(
                    causal_matrix,
                    target_indices,
                    indices.primary,
                    indices.beacon,
                    indices.clock_times,
                    config.beacon_separation,
                    spec.name,
                )
                accessible = chart.accessible
                accessible_count = int(np.count_nonzero(accessible))
                if accessible_count >= 3:
                    true_x = support_events[accessible, 1]
                    estimated_x = chart.reconstructed_coords[accessible, 1]
                    pairs = _sample_pairs(
                        accessible_count,
                        config.pair_count,
                        config.seed
                        + 1_000_000 * n_index
                        + 10_000 * tick_index
                        + repetition,
                    )
                    inversion = pair_distance_order_inversion_rate(
                        true_x,
                        estimated_x,
                        pairs,
                    )
                    used_pairs = pairs.shape[0]
                else:
                    inversion = float("nan")
                    used_pairs = 0
                rows.append(
                    {
                        "N": float(n_events),
                        "tick_count": float(tick_count),
                        "repetition": float(repetition),
                        "accessible_count": float(accessible_count),
                        "accessible_fraction": float(accessible_count / n_events),
                        "pair_count": float(used_pairs),
                        "distance_order_inversion_rate": float(inversion),
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write oriented chart distance-order rows."""

    output_path = (
        output_dir / "data" / "oriented_chart_distance_order_preservation.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save distance-order inversion versus tick count."""

    output_path = (
        output_dir
        / "figures"
        / "oriented_chart_distance_order_inversion_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    tick_values = sorted({row["tick_count"] for row in rows})
    for n_events in sorted({row["N"] for row in rows}):
        means = []
        for tick_count in tick_values:
            subset = [
                row
                for row in rows
                if row["N"] == n_events and row["tick_count"] == tick_count
            ]
            means.append(
                float(
                    np.nanmean(
                        [row["distance_order_inversion_rate"] for row in subset]
                    )
                )
            )
        ax.plot(tick_values, means, marker="o", label=f"N={int(n_events)}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Pair-distance order inversion/error rate")
    ax.set_title("Oriented chart spatial distance-order preservation")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    figure_path = save_plot(rows, config.output_dir)
    print(f"Wrote oriented chart distance-order data: {output_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
