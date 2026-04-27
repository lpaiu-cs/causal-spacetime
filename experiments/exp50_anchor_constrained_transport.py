"""Anchor-constrained cross-slice transport diagnostics."""

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
from causal_spacetime_lab.ordinal import pair_distance_values_1d
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    procrustes_align,
)
from causal_spacetime_lab.slice_anchors import (
    anchor_indices_by_slice,
    anchor_reference_positions_by_slice,
)
from causal_spacetime_lab.slice_local_embedding import (
    assemble_slice_embeddings_with_nan,
    fit_slice_local_ordinal_embeddings,
)
from causal_spacetime_lab.slice_transport import (
    align_slices_by_anchor_points,
    apply_random_slice_gauge_transforms,
)
from causal_spacetime_lab.sliced_distance_order import (
    pair_slice_labels,
    quadruplet_constraints_from_sliced_pair_distances,
)
from causal_spacetime_lab.spatial_slices import same_slice_unordered_pairs
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for anchor-constrained transport."""

    T: float = 2.0
    n_events: int = 800
    tick_count: int = 128
    bin_width: int = 4
    anchor_positions: tuple[float, ...] = (-0.35, 0.0, 0.35)
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
    constraint_count: int = 5000
    steps: int = 600
    restarts: int = 2
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

    parser = argparse.ArgumentParser(description="Anchor-constrained transport.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=800)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument("--bin-width", type=int, default=4)
    parser.add_argument(
        "--anchor-positions",
        nargs="+",
        default=["-0.35", "0.0", "0.35"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--constraint-count", type=int, default=5000)
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_count=args.tick_count,
        bin_width=args.bin_width,
        anchor_positions=_parse_float_values(args.anchor_positions, "anchor"),
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        constraint_count=args.constraint_count,
        steps=args.steps,
        restarts=args.restarts,
        output_dir=args.output_dir,
    )


def _evaluate(coords: np.ndarray, true_x: np.ndarray) -> tuple[float, float]:
    finite = np.isfinite(coords[: true_x.size])
    if np.count_nonzero(finite) < 3:
        return float("nan"), float("nan")
    aligned, diagnostics = procrustes_align(
        coords[: true_x.size][finite, None],
        true_x[finite, None],
    )
    return (
        diagnostics["rmse"],
        embedding_distance_order_error(aligned, true_x[finite, None], seed=3),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run anchor-constrained transport diagnostics."""

    rows: list[dict[str, float]] = []
    anchors = np.asarray(config.anchor_positions, dtype=float)
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
        anchor_coords = np.tile(anchors, len(slice_values))
        anchor_labels = np.repeat(slice_values, anchors.size)
        anchor_ids = np.tile(np.arange(anchors.size, dtype=int), len(slice_values))
        combined_coords = np.concatenate((result["reconstructed_X"], anchor_coords))
        combined_labels = np.concatenate((result["slice_labels"], anchor_labels))
        pairs = same_slice_unordered_pairs(
            combined_labels,
            max_pairs_per_slice=500,
            seed=config.seed + repetition,
        )
        distances = pair_distance_values_1d(combined_coords, pairs)
        pair_labels = pair_slice_labels(pairs, combined_labels)
        constraints = quadruplet_constraints_from_sliced_pair_distances(
            pairs,
            distances,
            pair_labels,
            config.constraint_count,
            seed=config.seed + 10_000 + repetition,
        )
        embeddings = fit_slice_local_ordinal_embeddings(
            combined_coords.size,
            combined_labels,
            constraints,
            min_points_per_slice=4,
            min_constraints_per_slice=20,
            steps=config.steps,
            restarts=config.restarts,
            learning_rate=0.1,
            seed=config.seed + repetition,
        )
        assembled = assemble_slice_embeddings_with_nan(
            combined_coords.size,
            embeddings,
        )[:, 0]
        anchor_global_indices = np.arange(
            config.n_events,
            combined_coords.size,
            dtype=int,
        )
        by_slice = anchor_indices_by_slice(anchor_global_indices, anchor_labels)
        references = anchor_reference_positions_by_slice(
            anchor_ids,
            anchors,
            anchor_labels,
        )
        anchored = align_slices_by_anchor_points(
            assembled,
            combined_labels,
            by_slice,
            references,
        )
        random_gauge = apply_random_slice_gauge_transforms(
            assembled,
            combined_labels,
            seed=config.seed + repetition,
        )
        no_transport_rmse, no_transport_order = _evaluate(assembled, support[:, 1])
        anchored_rmse, anchored_order = _evaluate(anchored, support[:, 1])
        random_rmse, random_order = _evaluate(random_gauge, support[:, 1])
        rows.append(
            {
                "repetition": float(repetition),
                "slice_count": float(len(slice_values)),
                "constraint_count": float(constraints.shape[0]),
                "embedded_slice_count": float(len(embeddings)),
                "no_transport_rmse": no_transport_rmse,
                "anchor_transport_rmse": anchored_rmse,
                "random_gauge_rmse": random_rmse,
                "no_transport_order_error": no_transport_order,
                "anchor_transport_order_error": anchored_order,
                "random_gauge_order_error": random_order,
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write anchor transport rows."""

    output_path = output_dir / "data" / "anchor_constrained_transport.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save anchor transport figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    path1 = figure_dir / "anchor_transport_global_rmse.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    keys = ["no_transport_rmse", "anchor_transport_rmse", "random_gauge_rmse"]
    means = [float(np.nanmean([row[key] for row in rows])) for key in keys]
    ax.bar(["no transport", "anchor", "random gauge"], means)
    ax.set_ylabel("Global RMSE after alignment")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path1, dpi=200)
    plt.close(fig)

    path2 = figure_dir / "anchor_transport_distance_order_error.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    keys = [
        "no_transport_order_error",
        "anchor_transport_order_error",
        "random_gauge_order_error",
    ]
    means = [float(np.nanmean([row[key] for row in rows])) for key in keys]
    ax.bar(["no transport", "anchor", "random gauge"], means)
    ax.set_ylabel("Distance-order error")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path2, dpi=200)
    plt.close(fig)
    return path1, path2


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote anchor-constrained transport data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
