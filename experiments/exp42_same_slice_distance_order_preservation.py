"""Same-slice spatial distance-order preservation diagnostics."""

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
    """Configuration for same-slice distance-order diagnostics."""

    T: float = 2.0
    n_values: tuple[int, ...] = (300, 600, 1200)
    tick_values: tuple[int, ...] = (32, 64, 128)
    bin_width_values: tuple[int, ...] = (1, 2, 4)
    repetitions: int = 5
    seed: int = 0
    beacon_separation: float = 0.15
    max_pairs_per_slice: int = 200
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

    parser = argparse.ArgumentParser(description="Same-slice distance-order checks.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--n-values", nargs="+", default=["300", "600", "1200"])
    parser.add_argument("--tick-values", nargs="+", default=["32", "64", "128"])
    parser.add_argument("--bin-width-values", nargs="+", default=["1", "2", "4"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--max-pairs-per-slice", type=int, default=200)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        bin_width_values=_parse_int_values(args.bin_width_values, "bin width"),
        repetitions=args.repetitions,
        seed=args.seed,
        beacon_separation=args.beacon_separation,
        max_pairs_per_slice=args.max_pairs_per_slice,
        output_dir=args.output_dir,
    )


def _sample_unordered_pairs(
    indices: np.ndarray,
    max_pairs: int,
    seed: int,
) -> np.ndarray:
    if indices.size < 2 or max_pairs < 1:
        return np.empty((0, 2), dtype=int)
    left, right = np.triu_indices(indices.size, k=1)
    pairs = np.column_stack((indices[left], indices[right])).astype(int)
    if pairs.shape[0] > max_pairs:
        rng = np.random.default_rng(seed)
        selected = rng.choice(pairs.shape[0], size=max_pairs, replace=False)
        pairs = pairs[selected]
    return pairs


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run same-slice distance-order preservation diagnostics."""

    rows: list[dict[str, float]] = []
    for n_events in config.n_values:
        for tick_count in config.tick_values:
            for bin_width in config.bin_width_values:
                for repetition in range(config.repetitions):
                    support = sprinkle_1p1_causal_diamond(
                        n_events,
                        T=config.T,
                        seed=config.seed + 100_000 * n_events + repetition,
                    )
                    result = reconstruct_stationary_oriented_slices_1p1(
                        support,
                        config.T,
                        tick_count,
                        config.beacon_separation,
                        bin_width=bin_width,
                    )
                    labels = result["slice_labels"]
                    pairs = same_slice_unordered_pairs(
                        labels,
                        max_pairs_per_slice=config.max_pairs_per_slice,
                        seed=config.seed + 10_000 * tick_count + repetition,
                    )
                    sizes = slice_sizes(labels)
                    accessible_indices = np.flatnonzero(result["accessible"])
                    baseline_count = max(pairs.shape[0], 1)
                    all_pairs = _sample_unordered_pairs(
                        accessible_indices,
                        baseline_count,
                        seed=config.seed + 20_000 * tick_count + repetition,
                    )
                    same_slice_error = sliced_pair_distance_order_inversion_rate(
                        support[:, 1],
                        result["reconstructed_X"],
                        pairs,
                    )
                    all_pairs_error = sliced_pair_distance_order_inversion_rate(
                        support[:, 1],
                        result["reconstructed_X"],
                        all_pairs,
                    )
                    rows.append(
                        {
                            "N": float(n_events),
                            "tick_count": float(tick_count),
                            "bin_width": float(bin_width),
                            "repetition": float(repetition),
                            "slice_count": float(len(sizes)),
                            "mean_slice_size": float(
                                np.mean(list(sizes.values())) if sizes else np.nan
                            ),
                            "pair_count": float(pairs.shape[0]),
                            "all_pairs_count": float(all_pairs.shape[0]),
                            "same_slice_inversion_rate": same_slice_error,
                            "all_pairs_inversion_rate": all_pairs_error,
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write same-slice distance-order rows."""

    output_path = output_dir / "data" / "same_slice_distance_order_preservation.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(
    rows: list[dict[str, float]],
    key: str,
    ylabel: str,
    output_path: Path,
) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for bin_width in sorted({row["bin_width"] for row in rows}):
        subset = [row for row in rows if row["bin_width"] == bin_width]
        ticks = sorted({row["tick_count"] for row in subset})
        values = [
            float(np.nanmean([row[key] for row in subset if row["tick_count"] == tick]))
            for tick in ticks
        ]
        ax.plot(ticks, values, marker="o", label=f"bin width={int(bin_width)}")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save same-slice distance-order figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    path1 = _plot(
        rows,
        "same_slice_inversion_rate",
        "Same-slice inversion rate",
        figure_dir / "same_slice_distance_order_inversion_vs_ticks.png",
    )
    path2 = figure_dir / "same_slice_vs_all_pairs_inversion.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    labels = ("same_slice_inversion_rate", "all_pairs_inversion_rate")
    ax.bar(
        ["same-slice", "all accessible"],
        [float(np.nanmean([row[label] for row in rows])) for label in labels],
    )
    ax.set_ylabel("Distance-order inversion rate")
    ax.set_title("Same-Slice vs All-Pairs Baseline")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path2, dpi=200)
    plt.close(fig)

    path3 = figure_dir / "same_slice_pair_count_vs_bin_width.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    widths = sorted({row["bin_width"] for row in rows})
    values = [
        float(np.nanmean([row["pair_count"] for row in rows if row["bin_width"] == w]))
        for w in widths
    ]
    ax.plot(widths, values, marker="o")
    ax.set_xlabel("Radar-time bin width")
    ax.set_ylabel("Same-slice pair count")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path3, dpi=200)
    plt.close(fig)
    return path1, path2, path3


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote same-slice distance-order data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
