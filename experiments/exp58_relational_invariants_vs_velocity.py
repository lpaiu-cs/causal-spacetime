"""Relational invariants versus transport-dependent velocity statements."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.persistence import finite_difference_velocity_by_object
from causal_spacetime_lab.relational_evolution import (
    apply_per_slice_affine_position_gauge,
    compare_histories_order_disagreement,
    relational_change_rate_between_slices,
    relational_shape_history,
)
from causal_spacetime_lab.relational_scenarios import (
    generate_shear_or_reordering_history_1d,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for relational-invariants versus velocity diagnostics."""

    slice_count: int = 12
    object_count: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Relational invariants vs velocity.")
    parser.add_argument("--slice-count", type=int, default=12)
    parser.add_argument("--object-count", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_count=args.object_count,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _flatten_history(
    history: dict[int, dict[int, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positions: list[float] = []
    slices: list[int] = []
    objects: list[int] = []
    for slice_label, object_positions in sorted(history.items()):
        for object_id, position in sorted(object_positions.items()):
            positions.append(position)
            slices.append(slice_label)
            objects.append(object_id)
    return (
        np.asarray(positions, dtype=float),
        np.asarray(slices, dtype=int),
        np.asarray(objects, dtype=int),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run relational invariants versus velocity diagnostics."""

    slice_labels = np.arange(config.slice_count)
    base_positions = np.linspace(-0.4, 0.4, config.object_count)
    displacement = np.linspace(-0.55, 0.55, config.slice_count)
    history = generate_shear_or_reordering_history_1d(
        slice_labels,
        base_positions,
        moving_object_id=config.object_count // 2,
        displacement_by_slice=displacement,
    )
    signature_history = relational_shape_history(history)
    change_rows = relational_change_rate_between_slices(signature_history)
    relational_change = float(
        np.mean([row["order_disagreement"] for row in change_rows])
    )
    gauged = apply_per_slice_affine_position_gauge(history, seed=config.seed)
    gauge_disagreement = compare_histories_order_disagreement(
        signature_history,
        relational_shape_history(gauged),
    )
    positions, slices, objects = _flatten_history(history)
    time_by_slice = {int(label): float(label) for label in slice_labels}
    transports = {
        "identity": positions,
        "drifting": positions + 0.08 * slices,
        "reflected": -positions,
    }
    rows: list[dict[str, float | str]] = []
    for name, transported_positions in transports.items():
        velocities = finite_difference_velocity_by_object(
            transported_positions,
            slices,
            objects,
            time_by_slice,
        )
        rows.append(
            {
                "transport": name,
                "relational_change_rate": relational_change,
                "gauge_history_disagreement": gauge_disagreement,
                "mean_velocity": float(np.mean(list(velocities.values()))),
                "velocity_spread": float(np.std(list(velocities.values()))),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write relational invariants versus velocity rows."""

    output_path = output_dir / "data" / "relational_invariants_vs_velocity.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save relational invariants versus velocity figure."""

    output_path = (
        output_dir / "figures" / "relational_change_vs_velocity_transport.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    labels = [str(row["transport"]) for row in rows]
    ax.plot(
        labels,
        [float(row["mean_velocity"]) for row in rows],
        marker="o",
        label="mean velocity",
    )
    ax.axhline(
        float(rows[0]["relational_change_rate"]),
        color="black",
        linestyle="--",
        label="relational change rate",
    )
    ax.set_ylabel("Diagnostic value")
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
    print(f"Wrote relational invariants versus velocity data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
