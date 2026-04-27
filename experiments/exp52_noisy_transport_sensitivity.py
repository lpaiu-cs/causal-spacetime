"""Sensitivity of derived cross-slice quantities to noisy anchor transport."""

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
from causal_spacetime_lab.ordinal_embedding import procrustes_align
from causal_spacetime_lab.persistence import (
    finite_difference_velocity_by_object,
    generate_persistent_object_events_1p1,
)
from causal_spacetime_lab.slice_anchors import (
    anchor_indices_by_slice,
    anchor_reference_positions_by_slice,
)
from causal_spacetime_lab.slice_transport import align_slices_by_anchor_points
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for noisy transport sensitivity."""

    T: float = 2.0
    n_events: int = 800
    tick_count: int = 128
    bin_width: int = 4
    anchor_positions: tuple[float, ...] = (-0.35, 0.0, 0.35)
    noise_levels: tuple[float, ...] = (0.0, 0.01, 0.03, 0.05, 0.10)
    repetitions: int = 5
    seed: int = 0
    beacon_separation: float = 0.15
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str], name: str) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Noisy transport sensitivity.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=800)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument("--bin-width", type=int, default=4)
    parser.add_argument(
        "--anchor-positions",
        nargs="+",
        default=["-0.35", "0.0", "0.35"],
    )
    parser.add_argument(
        "--noise-levels",
        nargs="+",
        default=["0.0", "0.01", "0.03", "0.05", "0.10"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_count=args.tick_count,
        bin_width=args.bin_width,
        anchor_positions=_parse_float_values(args.anchor_positions, "anchor"),
        noise_levels=_parse_float_values(args.noise_levels, "noise"),
        repetitions=args.repetitions,
        seed=args.seed,
        beacon_separation=args.beacon_separation,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run noisy anchor transport sensitivity diagnostics."""

    rows: list[dict[str, float]] = []
    anchors = np.asarray(config.anchor_positions, dtype=float)
    true_velocities = np.asarray([-0.12, 0.0, 0.11])
    initial_positions = np.asarray([-0.25, 0.0, 0.25])
    rng = np.random.default_rng(config.seed)
    for repetition in range(config.repetitions):
        support = sprinkle_1p1_causal_diamond(
            config.n_events,
            T=config.T,
            seed=config.seed + repetition,
        )
        result = reconstruct_stationary_oriented_slices_1p1(
            support,
            config.T,
            config.tick_count,
            config.beacon_separation,
            bin_width=config.bin_width,
        )
        slice_values = sorted(
            set(int(value) for value in result["slice_labels"] if value >= 0)
        )
        if not slice_values:
            continue
        anchor_coords = np.tile(anchors, len(slice_values))
        anchor_labels = np.repeat(slice_values, anchors.size)
        anchor_ids = np.tile(np.arange(anchors.size, dtype=int), len(slice_values))
        combined_coords = np.concatenate((result["reconstructed_X"], anchor_coords))
        combined_labels = np.concatenate((result["slice_labels"], anchor_labels))
        anchor_global_indices = np.arange(
            config.n_events,
            combined_coords.size,
            dtype=int,
        )
        by_slice = anchor_indices_by_slice(anchor_global_indices, anchor_labels)
        base_references = anchor_reference_positions_by_slice(
            anchor_ids,
            anchors,
            anchor_labels,
        )
        slice_times = np.asarray(slice_values, dtype=float)
        object_events, object_slices, object_ids = (
            generate_persistent_object_events_1p1(
                slice_times,
                initial_positions,
                true_velocities,
            )
        )
        object_slice_labels = np.asarray([slice_values[int(s)] for s in object_slices])
        true_velocity_by_id = {
            index: value for index, value in enumerate(true_velocities)
        }
        for noise_level in config.noise_levels:
            noisy_refs = {
                label: values + rng.normal(scale=noise_level, size=values.shape)
                for label, values in base_references.items()
            }
            aligned = align_slices_by_anchor_points(
                combined_coords,
                combined_labels,
                by_slice,
                noisy_refs,
            )
            support_aligned = aligned[: config.n_events]
            finite = np.isfinite(support_aligned)
            if np.count_nonzero(finite) >= 3:
                _, diagnostics = procrustes_align(
                    support_aligned[finite, None],
                    support[finite, 1, None],
                )
                global_rmse = diagnostics["rmse"]
            else:
                global_rmse = float("nan")
            object_refs = {
                label: values + rng.normal(scale=noise_level, size=values.shape)
                for label, values in base_references.items()
            }
            object_anchor_coords = np.concatenate((object_events[:, 1], anchor_coords))
            object_anchor_labels = np.concatenate((object_slice_labels, anchor_labels))
            object_anchor_indices = np.arange(
                object_events.shape[0],
                object_anchor_coords.size,
            )
            object_by_slice = anchor_indices_by_slice(
                object_anchor_indices,
                anchor_labels,
            )
            object_aligned = align_slices_by_anchor_points(
                object_anchor_coords,
                object_anchor_labels,
                object_by_slice,
                object_refs,
            )[: object_events.shape[0]]
            time_by_slice = {label: float(i) for i, label in enumerate(slice_values)}
            estimated = finite_difference_velocity_by_object(
                object_aligned,
                object_slice_labels,
                object_ids,
                time_by_slice,
            )
            errors = [
                abs(value - true_velocity_by_id[object_id])
                for object_id, value in estimated.items()
            ]
            rows.append(
                {
                    "noise_level": float(noise_level),
                    "repetition": float(repetition),
                    "global_rmse": global_rmse,
                    "velocity_instability": float(np.mean(errors))
                    if errors
                    else float("nan"),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write noisy transport rows."""

    output_path = output_dir / "data" / "noisy_transport_sensitivity.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(rows: list[dict[str, float]], key: str, ylabel: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    noise = sorted({row["noise_level"] for row in rows})
    values = [
        float(np.nanmean([row[key] for row in rows if row["noise_level"] == level]))
        for level in noise
    ]
    ax.plot(noise, values, marker="o")
    ax.set_xlabel("Anchor noise level")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save noisy transport figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "global_rmse",
            "Global RMSE",
            figure_dir / "noisy_transport_global_rmse.png",
        ),
        _plot(
            rows,
            "velocity_instability",
            "Velocity instability",
            figure_dir / "noisy_transport_velocity_instability.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote noisy transport sensitivity data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
