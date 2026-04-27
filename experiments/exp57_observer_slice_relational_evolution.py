"""Persistence-dependent relational histories on observer-derived slices."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.observer_slice_pipeline import (
    reconstruct_stationary_oriented_slices_1p1,
)
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
    """Configuration for observer-slice relational evolution."""

    T: float = 2.0
    slice_count: int = 10
    object_count: int = 5
    tick_count: int = 128
    bin_width: int = 4
    beacon_separation: float = 0.15
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Observer-slice relational history.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--slice-count", type=int, default=10)
    parser.add_argument("--object-count", type=int, default=5)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument("--bin-width", type=int, default=4)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        slice_count=args.slice_count,
        object_count=args.object_count,
        tick_count=args.tick_count,
        bin_width=args.bin_width,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _events_from_history(
    slice_times: np.ndarray,
    history: dict[int, dict[int, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    events: list[tuple[float, float]] = []
    true_slice_ids: list[int] = []
    object_ids: list[int] = []
    for true_slice_id, time in enumerate(slice_times):
        for object_id, position in sorted(history[int(true_slice_id)].items()):
            events.append((float(time), float(position)))
            true_slice_ids.append(true_slice_id)
            object_ids.append(object_id)
    return (
        np.asarray(events, dtype=float),
        np.asarray(true_slice_ids, dtype=int),
        np.asarray(object_ids, dtype=int),
    )


def _positions_by_derived_slice(
    values: np.ndarray,
    derived_labels: np.ndarray,
    object_ids: np.ndarray,
    accessible: np.ndarray,
    min_objects: int = 3,
) -> dict[int, dict[int, float]]:
    buckets: dict[int, dict[int, list[float]]] = {}
    for value, label, object_id, ok in zip(
        values,
        derived_labels,
        object_ids,
        accessible,
        strict=True,
    ):
        if not ok or label < 0 or not np.isfinite(value):
            continue
        buckets.setdefault(int(label), {}).setdefault(int(object_id), []).append(
            float(value)
        )
    output: dict[int, dict[int, float]] = {}
    for label, objects in buckets.items():
        averaged = {
            object_id: float(np.mean(position_values))
            for object_id, position_values in objects.items()
        }
        if len(averaged) >= min_objects:
            output[label] = averaged
    return output


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run observer-derived slice relational history diagnostics."""

    slice_times = np.linspace(-0.45, 0.45, config.slice_count)
    base_positions = np.linspace(-0.25, 0.25, config.object_count)
    displacement = np.linspace(-0.16, 0.16, config.slice_count)
    true_history = generate_shear_or_reordering_history_1d(
        np.arange(config.slice_count),
        base_positions,
        moving_object_id=config.object_count // 2,
        displacement_by_slice=displacement,
    )
    events, _, object_ids = _events_from_history(slice_times, true_history)
    result = reconstruct_stationary_oriented_slices_1p1(
        events,
        config.T,
        config.tick_count,
        config.beacon_separation,
        bin_width=config.bin_width,
    )
    accessible = result["accessible"] & np.isfinite(result["reconstructed_X"])
    reconstructed_history = _positions_by_derived_slice(
        result["reconstructed_X"],
        result["slice_labels"],
        object_ids,
        accessible,
    )
    validation_history = _positions_by_derived_slice(
        events[:, 1],
        result["slice_labels"],
        object_ids,
        accessible,
    )
    reconstructed_signatures = relational_shape_history(reconstructed_history)
    validation_signatures = relational_shape_history(validation_history)
    order_error = compare_histories_order_disagreement(
        reconstructed_signatures,
        validation_signatures,
    )
    gauged_history = apply_per_slice_affine_position_gauge(
        reconstructed_history,
        seed=config.seed,
    )
    gauge_error = compare_histories_order_disagreement(
        reconstructed_signatures,
        relational_shape_history(gauged_history),
    )
    return [
        {
            "slice_count": float(config.slice_count),
            "object_count": float(config.object_count),
            "observer_slice_count": float(len(reconstructed_history)),
            "accessible_event_count": float(np.count_nonzero(accessible)),
            "relational_history_order_error": order_error,
            "gauge_history_disagreement": gauge_error,
        }
    ]


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write observer-slice relational evolution rows."""

    output_path = output_dir / "data" / "observer_slice_relational_evolution.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save observer-slice relational evolution figure."""

    output_path = (
        output_dir / "figures" / "observer_slice_relational_evolution_error.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    row = rows[0]
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.bar(
        ["validation error", "gauge disagreement"],
        [row["relational_history_order_error"], row["gauge_history_disagreement"]],
    )
    ax.set_ylabel("Order-history disagreement")
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
    print(f"Wrote observer-slice relational evolution data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
