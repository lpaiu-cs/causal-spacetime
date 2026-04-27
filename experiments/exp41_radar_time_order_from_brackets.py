"""Radar-time order reconstruction from observer tick brackets."""

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
from causal_spacetime_lab.ordinal import order_inversion_rate
from causal_spacetime_lab.radar_order import radar_tick_brackets_from_order
from causal_spacetime_lab.spatial_slices import radar_time_rank_from_tick_brackets
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for radar-time order reconstruction."""

    T: float = 2.0
    n_values: tuple[int, ...] = (300, 600, 1200)
    tick_values: tuple[int, ...] = (16, 32, 64, 128)
    repetitions: int = 5
    seed: int = 0
    bin_width: int = 2
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

    parser = argparse.ArgumentParser(description="Radar-time order from brackets.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--n-values", nargs="+", default=["300", "600", "1200"])
    parser.add_argument("--tick-values", nargs="+", default=["16", "32", "64", "128"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--bin-width", type=int, default=2)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        repetitions=args.repetitions,
        seed=args.seed,
        bin_width=args.bin_width,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run radar-time order reconstruction diagnostics."""

    rows: list[dict[str, float]] = []
    for n_events in config.n_values:
        for tick_count in config.tick_values:
            for repetition in range(config.repetitions):
                support = sprinkle_1p1_causal_diamond(
                    n_events,
                    T=config.T,
                    seed=config.seed + 100_000 * n_events + repetition,
                )
                observer_events, _ = make_stationary_observer_chain_1p1(
                    config.T,
                    tick_count,
                )
                combined = np.vstack((support, observer_events))
                observer_indices = observer_chain_indices(n_events, tick_count)
                targets = np.arange(n_events, dtype=int)
                predecessor, successor, accessible = radar_tick_brackets_from_order(
                    causal_matrix_1p1(combined),
                    observer_indices,
                    targets,
                )
                ranks = radar_time_rank_from_tick_brackets(
                    predecessor,
                    successor,
                    accessible,
                )
                if np.count_nonzero(accessible) >= 2:
                    inversion = order_inversion_rate(
                        support[accessible, 0],
                        ranks[accessible],
                    )
                else:
                    inversion = float("nan")
                rows.append(
                    {
                        "N": float(n_events),
                        "tick_count": float(tick_count),
                        "repetition": float(repetition),
                        "accessible_count": float(np.count_nonzero(accessible)),
                        "accessible_fraction": float(np.mean(accessible)),
                        "radar_time_order_inversion_rate": inversion,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write radar-time order rows."""

    output_path = output_dir / "data" / "radar_time_order_from_brackets.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save radar-time inversion figure."""

    output_path = output_dir / "figures" / "radar_time_order_inversion_vs_ticks.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for n_events in sorted({row["N"] for row in rows}):
        subset = [row for row in rows if row["N"] == n_events]
        ticks = sorted({row["tick_count"] for row in subset})
        values = [
            float(
                np.nanmean(
                    [
                        row["radar_time_order_inversion_rate"]
                        for row in subset
                        if row["tick_count"] == tick
                    ]
                )
            )
            for tick in ticks
        ]
        ax.plot(ticks, values, marker="o", label=f"N={int(n_events)}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Radar-time order inversion rate")
    ax.set_title("Radar-Time Order From Tick Brackets")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_path = save_figure(rows, config.output_dir)
    print(f"Wrote radar-time order data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
