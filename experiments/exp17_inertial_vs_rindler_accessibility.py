"""Compare inertial and Rindler two-way radar accessibility."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.discrete_radar import discrete_radar_coordinates_from_order
from causal_spacetime_lab.horizon_validation import (
    finite_coverage_mask_from_analytic_ticks,
)
from causal_spacetime_lab.observer import (
    make_stationary_observer_chain_1p1,
    observer_chain_indices,
)
from causal_spacetime_lab.rindler import (
    analytic_rindler_radar_ticks_1p1,
    make_rindler_observer_chain_1p1,
    safe_tau_range_for_rindler_chain_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

OUTPUT_DATA = Path("outputs/data/inertial_vs_rindler_accessibility.csv")
OUTPUT_FIGURE = Path("outputs/figures/inertial_vs_rindler_accessibility.png")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for inertial versus Rindler accessibility comparison."""

    T: float = 4.0
    n_events: int = 800
    tick_count: int = 128
    acceleration: float = 2.0
    seed: int = 0
    direction: int = 1
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Compare inertial and Rindler two-way radar accessibility."
    )
    parser.add_argument("--T", type=float, default=4.0)
    parser.add_argument("--n-events", type=int, default=800)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument("--acceleration", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--direction", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.n_events,
        tick_count=args.tick_count,
        acceleration=args.acceleration,
        seed=args.seed,
        direction=args.direction,
        output_dir=args.output_dir,
    )


def _accessible_from_chain(
    support_events: np.ndarray,
    observer_events: np.ndarray,
    clock_times: np.ndarray,
) -> np.ndarray:
    combined = np.vstack((support_events, observer_events))
    causal_matrix = causal_matrix_1p1(combined)
    observer_indices = observer_chain_indices(support_events.shape[0], clock_times.size)
    target_indices = np.arange(support_events.shape[0], dtype=int)
    results = discrete_radar_coordinates_from_order(
        causal_matrix,
        observer_indices,
        target_indices,
        clock_times,
    )
    return np.asarray([result.accessible for result in results], dtype=bool)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run the inertial versus Rindler accessibility comparison."""

    support_events = sprinkle_1p1_causal_diamond(
        config.n_events,
        T=config.T,
        seed=config.seed,
    )
    inertial_events, inertial_clocks = make_stationary_observer_chain_1p1(
        config.T,
        config.tick_count,
    )
    inertial_accessible = _accessible_from_chain(
        support_events,
        inertial_events,
        inertial_clocks,
    )

    tau_min, tau_max = safe_tau_range_for_rindler_chain_1p1(
        config.T,
        config.acceleration,
        direction=config.direction,
    )
    rindler_events, rindler_clocks = make_rindler_observer_chain_1p1(
        config.acceleration,
        config.tick_count,
        tau_min,
        tau_max,
        direction=config.direction,
    )
    rindler_accessible = _accessible_from_chain(
        support_events,
        rindler_events,
        rindler_clocks,
    )
    tau_minus, tau_plus, wedge_accessible = analytic_rindler_radar_ticks_1p1(
        support_events,
        config.acceleration,
        direction=config.direction,
    )
    finite_coverage = finite_coverage_mask_from_analytic_ticks(
        tau_minus,
        tau_plus,
        tau_min,
        tau_max,
    )

    rows: list[dict[str, float]] = []
    for index, event in enumerate(support_events):
        rows.append(
            {
                "target_index": float(index),
                "t": float(event[0]),
                "x": float(event[1]),
                "inertial_accessible": float(inertial_accessible[index]),
                "rindler_accessible": float(rindler_accessible[index]),
                "rindler_wedge_accessible": float(wedge_accessible[index]),
                "rindler_finite_coverage_accessible": float(finite_coverage[index]),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float]],
    output_path: Path,
) -> Path:
    """Write comparison CSV rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(rows: list[dict[str, float]], output_path: Path) -> Path:
    """Save inertial versus Rindler accessibility map."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    t = np.asarray([row["t"] for row in rows], dtype=float)
    x = np.asarray([row["x"] for row in rows], dtype=float)
    inertial = np.asarray([row["inertial_accessible"] for row in rows], dtype=bool)
    rindler = np.asarray([row["rindler_accessible"] for row in rows], dtype=bool)
    colors = np.where(rindler, "tab:green", np.where(inertial, "tab:blue", "0.7"))

    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    ax.scatter(t, x, c=colors, s=12, alpha=0.7)
    line = np.linspace(float(np.min(t)), float(np.max(t)), 200)
    ax.plot(line, line, color="black", linewidth=0.9)
    ax.plot(line, -line, color="black", linewidth=0.9)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    ax.set_title("Inertial versus Rindler two-way accessibility")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(
        rows,
        config.output_dir / OUTPUT_DATA.relative_to("outputs"),
    )
    figure_path = save_plot(
        rows,
        config.output_dir / OUTPUT_FIGURE.relative_to("outputs"),
    )
    print(f"Wrote comparison data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
