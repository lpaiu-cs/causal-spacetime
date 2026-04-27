"""Reflection degeneracy of single-observer radar distance."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.oriented_radar import signed_position_from_two_distances
from causal_spacetime_lab.radar import stationary_radar_coordinates_1p1

BEACON_SEPARATION = 0.35
OUTPUT_DATA = Path("outputs/data/single_observer_reflection_degeneracy.csv")
OUTPUT_FIGURE = Path("outputs/figures/single_observer_reflection_degeneracy.png")


def deterministic_events() -> np.ndarray:
    """Return mirrored deterministic 1+1D events for the degeneracy check."""

    rows: list[tuple[float, float]] = []
    for t, x in [(-0.45, 0.10), (-0.15, 0.25), (0.20, 0.40), (0.45, 0.20)]:
        rows.append((t, x))
        rows.append((t, -x))
    return np.asarray(rows, dtype=float)


def run_experiment() -> list[dict[str, float]]:
    """Compute single-observer and oriented two-chain coordinates."""

    rows: list[dict[str, float]] = []
    for index, event in enumerate(deterministic_events()):
        radar = stationary_radar_coordinates_1p1(event)
        beacon_distance = abs(float(event[1]) - BEACON_SEPARATION)
        oriented_x = signed_position_from_two_distances(
            radar.radar_distance,
            beacon_distance,
            BEACON_SEPARATION,
        )
        rows.append(
            {
                "event_index": float(index),
                "t": float(event[0]),
                "true_x": float(event[1]),
                "single_observer_radar_time": radar.radar_time,
                "single_observer_radar_distance": radar.radar_distance,
                "beacon_distance": beacon_distance,
                "two_chain_signed_position": oriented_x,
            }
        )
    return rows


def write_summary(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write deterministic reflection-degeneracy rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(rows: list[dict[str, float]], output_path: Path = OUTPUT_FIGURE) -> Path:
    """Save reflection-degeneracy diagnostic plot."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    true_x = np.asarray([row["true_x"] for row in rows])
    single_distance = np.asarray(
        [row["single_observer_radar_distance"] for row in rows]
    )
    oriented_x = np.asarray([row["two_chain_signed_position"] for row in rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.scatter(true_x, single_distance, marker="o", label="Single-observer distance")
    ax.scatter(true_x, oriented_x, marker="s", label="Two-chain signed coordinate")
    lower = float(min(np.min(true_x), np.min(oriented_x)))
    upper = float(max(np.max(true_x), np.max(oriented_x)))
    ax.plot([lower, upper], [lower, upper], color="black", linewidth=1.0)
    ax.axhline(0.0, color="0.7", linewidth=0.8)
    ax.set_xlabel("True signed x")
    ax.set_ylabel("Reconstructed value")
    ax.set_title("Single-observer radar has reflection degeneracy")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    rows = run_experiment()
    data_path = write_summary(rows)
    figure_path = save_plot(rows)
    print(f"Wrote summary: {data_path}")
    print(f"Wrote figure: {figure_path}")
    print(
        "A single observer gives unsigned distance; the beacon separation is "
        "additional orientation structure."
    )


if __name__ == "__main__":
    main()
