"""Relational shape history from persistence plus slice-local distance order."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.relational_evolution import (
    relational_change_rate_between_slices,
    relational_shape_history,
)
from causal_spacetime_lab.relational_scenarios import (
    generate_expanding_configuration_history,
    generate_shear_or_reordering_history_1d,
    generate_static_configuration_history,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for relational shape-history diagnostics."""

    slice_count: int = 12
    object_count: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Relational shape histories.")
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


def _mean_change(history: dict[int, object]) -> tuple[float, list[dict[str, float]]]:
    rows = relational_change_rate_between_slices(history)  # type: ignore[arg-type]
    values = [row["order_disagreement"] for row in rows]
    return (float(np.mean(values)) if values else 0.0, rows)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run relational shape-history diagnostics without transport."""

    if config.slice_count < 2 or config.object_count < 3:
        raise ValueError("slice_count >= 2 and object_count >= 3 are required")
    slice_labels = np.arange(config.slice_count)
    base_positions = np.linspace(-0.4, 0.4, config.object_count)
    expansion = np.linspace(0.7, 1.6, config.slice_count)
    displacement = np.linspace(-0.7, 0.7, config.slice_count)
    histories = {
        "static": generate_static_configuration_history(slice_labels, base_positions),
        "uniform_expansion": generate_expanding_configuration_history(
            slice_labels,
            base_positions,
            expansion,
        ),
        "reordering": generate_shear_or_reordering_history_1d(
            slice_labels,
            base_positions,
            moving_object_id=config.object_count // 2,
            displacement_by_slice=displacement,
        ),
    }
    rows: list[dict[str, float | str]] = []
    for name, positions_by_slice in histories.items():
        history = relational_shape_history(positions_by_slice)
        mean_change, change_rows = _mean_change(history)
        rows.append(
            {
                "scenario": name,
                "slice_count": float(config.slice_count),
                "object_count": float(config.object_count),
                "mean_change_rate": mean_change,
                "max_change_rate": float(
                    max((row["order_disagreement"] for row in change_rows), default=0.0)
                ),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write relational shape-history rows."""

    output_path = output_dir / "data" / "relational_shape_history_without_transport.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save relational shape-change figure."""

    output_path = output_dir / "figures" / "relational_shape_change_rates.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.bar(
        [str(row["scenario"]) for row in rows],
        [float(row["mean_change_rate"]) for row in rows],
    )
    ax.set_ylabel("Mean pair-order change rate")
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
    print(f"Wrote relational shape history data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
