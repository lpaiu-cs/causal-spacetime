"""Gauge invariance of relational pair-distance order histories."""

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
    relational_shape_history,
)
from causal_spacetime_lab.relational_scenarios import (
    generate_shear_or_reordering_history_1d,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for relational gauge-invariance diagnostics."""

    slice_count: int = 12
    object_count: int = 6
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Relational gauge invariance.")
    parser.add_argument("--slice-count", type=int, default=12)
    parser.add_argument("--object-count", type=int, default=6)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_count=args.object_count,
        repetitions=args.repetitions,
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run relational gauge-invariance diagnostics."""

    slice_labels = np.arange(config.slice_count)
    base_positions = np.linspace(-0.45, 0.45, config.object_count)
    displacement = np.linspace(-0.4, 0.5, config.slice_count)
    base_history = generate_shear_or_reordering_history_1d(
        slice_labels,
        base_positions,
        moving_object_id=config.object_count // 2,
        displacement_by_slice=displacement,
    )
    base_signature_history = relational_shape_history(base_history)
    positions, slices, objects = _flatten_history(base_history)
    time_by_slice = {int(label): float(label) for label in slice_labels}
    base_velocities = finite_difference_velocity_by_object(
        positions,
        slices,
        objects,
        time_by_slice,
    )
    rows: list[dict[str, float]] = []
    for repetition in range(config.repetitions):
        gauged = apply_per_slice_affine_position_gauge(
            base_history,
            seed=config.seed + repetition,
        )
        gauged_signature_history = relational_shape_history(gauged)
        disagreement = compare_histories_order_disagreement(
            base_signature_history,
            gauged_signature_history,
        )
        g_positions, g_slices, g_objects = _flatten_history(gauged)
        gauged_velocities = finite_difference_velocity_by_object(
            g_positions,
            g_slices,
            g_objects,
            time_by_slice,
        )
        velocity_changes = [
            abs(gauged_velocities[object_id] - base_velocities[object_id])
            for object_id in sorted(set(base_velocities) & set(gauged_velocities))
        ]
        rows.append(
            {
                "repetition": float(repetition),
                "history_order_disagreement": disagreement,
                "mean_velocity_change": float(np.mean(velocity_changes))
                if velocity_changes
                else float("nan"),
                "max_velocity_change": float(np.max(velocity_changes))
                if velocity_changes
                else float("nan"),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write relational gauge-invariance rows."""

    output_path = output_dir / "data" / "relational_history_gauge_invariance.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save relational gauge-invariance figure."""

    output_path = output_dir / "figures" / "relational_history_gauge_invariance.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.bar(
        ["history disagreement", "velocity change"],
        [
            float(np.mean([row["history_order_disagreement"] for row in rows])),
            float(np.mean([row["mean_velocity_change"] for row in rows])),
        ],
    )
    ax.set_ylabel("Mean diagnostic value")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_path = save_figure(rows, config.output_dir)
    print(f"Wrote relational history gauge data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
