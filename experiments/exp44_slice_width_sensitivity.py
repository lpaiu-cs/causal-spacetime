"""Sensitivity of same-slice distance order to radar-time bin width."""

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
from causal_spacetime_lab.sliced_distance_order import (
    sliced_pair_distance_order_inversion_rate,
)
from causal_spacetime_lab.spatial_slices import (
    same_slice_unordered_pairs,
    slice_sizes,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for slice-width sensitivity."""

    T: float = 2.0
    n_events: int = 800
    tick_count: int = 128
    bin_width_values: tuple[int, ...] = (1, 2, 4, 8, 16)
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
    max_pairs_per_slice: int = 300
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str], name: str) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Slice-width sensitivity.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=800)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument(
        "--bin-width-values",
        nargs="+",
        default=["1", "2", "4", "8", "16"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-pairs-per-slice", type=int, default=300)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_count=args.tick_count,
        bin_width_values=_parse_int_values(args.bin_width_values, "bin width"),
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        max_pairs_per_slice=args.max_pairs_per_slice,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run slice-width sensitivity diagnostics."""

    rows: list[dict[str, float]] = []
    for repetition in range(config.repetitions):
        support = sprinkle_1p1_causal_diamond(
            config.n_events,
            T=config.T,
            seed=config.seed + repetition,
        )
        for bin_width in config.bin_width_values:
            result = reconstruct_stationary_oriented_slices_1p1(
                support,
                config.T,
                config.tick_count,
                config.beacon_separation,
                bin_width=bin_width,
            )
            labels = result["slice_labels"]
            sizes = slice_sizes(labels)
            pairs = same_slice_unordered_pairs(
                labels,
                max_pairs_per_slice=config.max_pairs_per_slice,
                seed=config.seed + 10_000 * bin_width + repetition,
            )
            error = sliced_pair_distance_order_inversion_rate(
                support[:, 1],
                result["reconstructed_X"],
                pairs,
            )
            size_values = list(sizes.values())
            rows.append(
                {
                    "bin_width": float(bin_width),
                    "repetition": float(repetition),
                    "slice_count": float(len(sizes)),
                    "mean_slice_size": float(
                        np.mean(size_values) if size_values else np.nan
                    ),
                    "median_slice_size": float(
                        np.median(size_values) if size_values else np.nan
                    ),
                    "pair_count": float(pairs.shape[0]),
                    "same_slice_inversion_rate": error,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write slice-width sensitivity rows."""

    output_path = output_dir / "data" / "slice_width_sensitivity.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(rows: list[dict[str, float]], key: str, ylabel: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    widths = sorted({row["bin_width"] for row in rows})
    values = [
        float(np.nanmean([row[key] for row in rows if row["bin_width"] == width]))
        for width in widths
    ]
    ax.plot(widths, values, marker="o")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Radar-time bin width")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save slice-width sensitivity figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "pair_count",
            "Same-slice pair count",
            figure_dir / "slice_width_pair_count.png",
        ),
        _plot(
            rows,
            "same_slice_inversion_rate",
            "Same-slice inversion rate",
            figure_dir / "slice_width_distance_order_error.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote slice-width sensitivity data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
