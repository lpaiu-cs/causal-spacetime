"""Same-slice constraint graph decomposition diagnostics."""

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
from causal_spacetime_lab.slice_constraint_graph import (
    constraint_cross_slice_fraction,
    slice_component_summary,
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
    """Configuration for sliced constraint graph decomposition."""

    T: float = 2.0
    n_events: int = 600
    tick_count: int = 128
    bin_width_values: tuple[int, ...] = (1, 2, 4, 8)
    constraint_count: int = 3000
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
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

    parser = argparse.ArgumentParser(description="Sliced constraint graph diagnostics.")
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=600)
    parser.add_argument("--tick-count", type=int, default=128)
    parser.add_argument("--bin-width-values", nargs="+", default=["1", "2", "4", "8"])
    parser.add_argument("--constraint-count", type=int, default=3000)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_count=args.tick_count,
        bin_width_values=_parse_int_values(args.bin_width_values, "bin width"),
        constraint_count=args.constraint_count,
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _constraints_for_result(
    result: dict[str, np.ndarray],
    constraint_count: int,
    seed: int,
) -> np.ndarray:
    pairs = same_slice_unordered_pairs(
        result["slice_labels"],
        max_pairs_per_slice=400,
        seed=seed,
    )
    if pairs.shape[0] < 2:
        return np.empty((0, 4), dtype=int)
    distances = pair_distance_values_1d(result["reconstructed_X"], pairs)
    pair_labels = pair_slice_labels(pairs, result["slice_labels"])
    return quadruplet_constraints_from_sliced_pair_distances(
        pairs,
        distances,
        pair_labels,
        constraint_count,
        seed=seed + 1,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run constraint graph decomposition diagnostics."""

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
            constraints = _constraints_for_result(
                result,
                config.constraint_count,
                seed=config.seed + 10_000 * bin_width + repetition,
            )
            summary = slice_component_summary(
                result["slice_labels"],
                constraints,
                config.n_events,
            )
            rows.append(
                {
                    "bin_width": float(bin_width),
                    "repetition": float(repetition),
                    "constraint_count": float(constraints.shape[0]),
                    "cross_slice_fraction": constraint_cross_slice_fraction(
                        constraints,
                        result["slice_labels"],
                    )
                    if constraints.size
                    else float("nan"),
                    **summary,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write sliced graph decomposition rows."""

    output_path = output_dir / "data" / "sliced_constraint_graph_decomposition.csv"
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
    """Save graph decomposition figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "component_count",
            "Constraint component count",
            figure_dir / "sliced_constraint_components_vs_bin_width.png",
        ),
        _plot(
            rows,
            "largest_component_size",
            "Largest component size",
            figure_dir / "sliced_constraint_largest_component_vs_bin_width.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote sliced constraint graph data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
