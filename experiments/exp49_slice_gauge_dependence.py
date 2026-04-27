"""Gauge dependence of cross-slice spatial judgments."""

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
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    procrustes_align,
)
from causal_spacetime_lab.slice_transport import (
    SliceTransportRule,
    apply_random_slice_gauge_transforms,
    same_position_under_transport,
)
from causal_spacetime_lab.sliced_distance_order import (
    sliced_pair_distance_order_inversion_rate,
)
from causal_spacetime_lab.spatial_slices import same_slice_unordered_pairs
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for slice gauge-dependence diagnostics."""

    T: float = 2.0
    n_events: int = 800
    tick_count: int = 128
    bin_width: int = 4
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Slice gauge-dependence checks.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=800)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument("--bin-width", type=int, default=4)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_count=args.tick_count,
        bin_width=args.bin_width,
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _all_pairs(indices: np.ndarray, count: int, seed: int) -> np.ndarray:
    if indices.size < 2:
        return np.empty((0, 2), dtype=int)
    left, right = np.triu_indices(indices.size, k=1)
    pairs = np.column_stack((indices[left], indices[right])).astype(int)
    if pairs.shape[0] > count:
        rng = np.random.default_rng(seed)
        pairs = pairs[rng.choice(pairs.shape[0], size=count, replace=False)]
    return pairs


def _same_slice_local_error(
    baseline: np.ndarray,
    transformed: np.ndarray,
    labels: np.ndarray,
    seed: int,
) -> float:
    values: list[float] = []
    for label in sorted(set(int(value) for value in labels if value >= 0)):
        pairs = same_slice_unordered_pairs(
            np.where(labels == label, label, -1),
            max_pairs_per_slice=300,
            seed=seed + label,
        )
        if pairs.shape[0] >= 2:
            values.append(
                sliced_pair_distance_order_inversion_rate(
                    baseline,
                    transformed,
                    pairs,
                )
            )
    return float(np.nanmean(values)) if values else float("nan")


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run slice gauge-dependence diagnostics."""

    rows: list[dict[str, float]] = []
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
        labels = result["slice_labels"]
        baseline = result["reconstructed_X"]
        gauged = apply_random_slice_gauge_transforms(
            baseline,
            labels,
            seed=config.seed + repetition,
        )
        accessible = np.flatnonzero(result["accessible"] & np.isfinite(gauged))
        pairs = same_slice_unordered_pairs(
            labels,
            max_pairs_per_slice=300,
            seed=repetition,
        )
        all_pairs = _all_pairs(accessible, max(pairs.shape[0], 1), seed=repetition)
        baseline_aligned, base_diag = procrustes_align(
            baseline[accessible, None],
            support[accessible, 1, None],
        )
        gauged_aligned, gauge_diag = procrustes_align(
            gauged[accessible, None],
            support[accessible, 1, None],
        )
        label_values = sorted(set(int(value) for value in labels if value >= 0))
        identity_rule = SliceTransportRule(
            name="identity",
            description="identity per-slice rule",
            slice_labels=labels,
            scale_by_slice={label: 1.0 for label in label_values},
            shift_by_slice={label: 0.0 for label in label_values},
            reflection_by_slice={label: 1.0 for label in label_values},
        )
        shifted_rule = SliceTransportRule(
            name="shifted",
            description="per-slice shifted rule",
            slice_labels=labels,
            scale_by_slice={label: 1.0 for label in label_values},
            shift_by_slice={label: 0.1 * label for label in label_values},
            reflection_by_slice={label: 1.0 for label in label_values},
        )
        if len(label_values) >= 2:
            same_identity = same_position_under_transport(
                0.0,
                label_values[0],
                0.0,
                label_values[1],
                identity_rule,
            )
            same_shifted = same_position_under_transport(
                0.0,
                label_values[0],
                0.0,
                label_values[1],
                shifted_rule,
            )
        else:
            same_identity = same_shifted = same_position_under_transport(
                0.0,
                -1,
                0.0,
                -2,
                identity_rule,
            )
        rows.append(
            {
                "repetition": float(repetition),
                "same_slice_order_error_after_gauge": _same_slice_local_error(
                    baseline,
                    gauged,
                    labels,
                    seed=config.seed + repetition,
                ),
                "baseline_global_rmse": base_diag["rmse"],
                "gauged_global_rmse": gauge_diag["rmse"],
                "baseline_all_pairs_order_error": embedding_distance_order_error(
                    baseline_aligned,
                    support[accessible, 1, None],
                    seed=config.seed + repetition,
                )
                if accessible.size >= 2
                else float("nan"),
                "gauged_all_pairs_order_error": (
                    sliced_pair_distance_order_inversion_rate(
                        support[:, 1],
                        gauged,
                        all_pairs,
                    )
                ),
                "same_position_identity": float(bool(same_identity.value)),
                "same_position_shifted": float(bool(same_shifted.value)),
                "same_slice_pair_count": float(pairs.shape[0]),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write gauge-dependence rows."""

    output_path = output_dir / "data" / "slice_gauge_dependence.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save gauge-dependence figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    path1 = figure_dir / "slice_gauge_same_slice_vs_global_error.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    keys = [
        "same_slice_order_error_after_gauge",
        "baseline_all_pairs_order_error",
        "gauged_all_pairs_order_error",
    ]
    ax.bar(
        ["same-slice", "global baseline", "global gauged"],
        [float(np.nanmean([row[key] for row in rows])) for key in keys],
    )
    ax.set_ylabel("Order error")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path1, dpi=200)
    plt.close(fig)

    path2 = figure_dir / "slice_gauge_cross_slice_judgments.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.bar(
        ["identity", "shifted"],
        [
            float(np.nanmean([row["same_position_identity"] for row in rows])),
            float(np.nanmean([row["same_position_shifted"] for row in rows])),
        ],
    )
    ax.set_ylabel("Same-position true fraction")
    ax.set_ylim(-0.05, 1.05)
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
    print(f"Wrote slice gauge-dependence data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
