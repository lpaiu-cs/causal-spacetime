"""Finite-speed lattice counterexample to Lorentz-invariance from speed alone."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.lattice import (
    cumulative_counts_by_time,
    edge_displacements,
    lattice_cumulative_counts,
    regular_lattice_causal_graph_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_forward_cone

T_STEPS = 30
SEED = 20260424
OUTPUT_DATA = Path("outputs/data/finite_speed_lattice_growth.csv")
OUTPUT_CONE_FIGURE = Path("outputs/figures/finite_speed_lattice_cones.png")
OUTPUT_GROWTH_FIGURE = Path("outputs/figures/finite_speed_lattice_count_growth.png")


def run_experiment() -> tuple[list[dict[str, float]], np.ndarray, np.ndarray]:
    """Build lattice and sprinkled comparison data."""

    graph = regular_lattice_causal_graph_1p1(T_STEPS)
    lattice_total = graph.events.shape[0]
    sprinkled_interior = sprinkle_1p1_forward_cone(
        lattice_total - 1,
        T_STEPS,
        seed=SEED,
    )
    sprinkled_events = np.vstack(([0.0, 0.0], sprinkled_interior))

    times = np.arange(T_STEPS + 1, dtype=float)
    lattice_counts = lattice_cumulative_counts(T_STEPS)
    sprinkled_counts = cumulative_counts_by_time(sprinkled_events, times)
    rho = (lattice_total - 1) / float(T_STEPS * T_STEPS)
    continuum_expected = 1.0 + rho * times * times

    rows: list[dict[str, float]] = []
    for time, lattice_count, sprinkled_count, expected_count in zip(
        times,
        lattice_counts,
        sprinkled_counts,
        continuum_expected,
        strict=True,
    ):
        rows.append(
            {
                "time": float(time),
                "lattice_cumulative_count": float(lattice_count),
                "sprinkled_cumulative_count": float(sprinkled_count),
                "continuum_expected_count": float(expected_count),
            }
        )

    return rows, graph.events.astype(float), sprinkled_events


def write_summary(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write event-count growth summary to CSV."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_cone_plot(
    lattice_events: np.ndarray,
    sprinkled_events: np.ndarray,
    output_path: Path = OUTPUT_CONE_FIGURE,
) -> Path:
    """Save side-by-side causal cone plots for lattice and sprinkle."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    graph = regular_lattice_causal_graph_1p1(T_STEPS)

    fig, (ax_lattice, ax_sprinkle) = plt.subplots(
        1,
        2,
        figsize=(11.0, 5.2),
        sharex=True,
        sharey=True,
    )

    for source, target in graph.edges:
        t0, x0 = graph.events[source]
        t1, x1 = graph.events[target]
        ax_lattice.plot([x0, x1], [t0, t1], color="0.70", linewidth=0.45)
    ax_lattice.scatter(
        lattice_events[:, 1],
        lattice_events[:, 0],
        s=10,
        color="#1f77b4",
        zorder=3,
    )
    ax_lattice.set_title("Regular finite-speed lattice")

    ax_sprinkle.scatter(
        sprinkled_events[:, 1],
        sprinkled_events[:, 0],
        s=10,
        color="#2ca02c",
        alpha=0.85,
    )
    ax_sprinkle.set_title("Poisson sprinkle in continuum cone")

    for ax in (ax_lattice, ax_sprinkle):
        ax.plot([-T_STEPS, 0], [T_STEPS, 0], color="black", linewidth=1.0)
        ax.plot([T_STEPS, 0], [T_STEPS, 0], color="black", linewidth=1.0)
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.2)
        ax.set_aspect("equal", adjustable="box")
    ax_lattice.set_ylabel("t")
    fig.suptitle("Finite signal speed does not remove lattice-preferred structure")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_growth_plot(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_GROWTH_FIGURE,
) -> Path:
    """Save cumulative event-count growth comparison."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    times = np.asarray([row["time"] for row in rows])
    lattice_counts = np.asarray([row["lattice_cumulative_count"] for row in rows])
    sprinkled_counts = np.asarray([row["sprinkled_cumulative_count"] for row in rows])
    expected_counts = np.asarray([row["continuum_expected_count"] for row in rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(times, lattice_counts, marker="o", markersize=3, label="Lattice cone")
    ax.plot(times, sprinkled_counts, marker="s", markersize=3, label="Poisson sample")
    ax.plot(times, expected_counts, linestyle="--", label="Continuum expectation")
    ax.set_xlabel("Coordinate time t")
    ax.set_ylabel("Cumulative event count")
    ax.set_title("Causal cone count growth")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def print_direction_summary() -> None:
    """Print the available edge directions in the lattice graph."""

    graph = regular_lattice_causal_graph_1p1(T_STEPS)
    displacements = edge_displacements(graph)
    unique, counts = np.unique(displacements, axis=0, return_counts=True)
    pairs = ", ".join(
        f"(dt={int(dt)}, dx={int(dx)}): {int(count)}"
        for (dt, dx), count in zip(unique, counts, strict=True)
    )
    print(f"Lattice edge directions: {pairs}")


def main() -> None:
    rows, lattice_events, sprinkled_events = run_experiment()
    summary_path = write_summary(rows)
    cone_path = save_cone_plot(lattice_events, sprinkled_events)
    growth_path = save_growth_plot(rows)
    print_direction_summary()
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote cone figure: {cone_path}")
    print(f"Wrote growth figure: {growth_path}")


if __name__ == "__main__":
    main()
