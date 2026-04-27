"""Exact sanity checks for spatial slice utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.sliced_distance_order import (
    pair_slice_labels,
    quadruplet_constraints_from_sliced_pair_distances,
    sliced_pair_distance_order_inversion_rate,
)
from causal_spacetime_lab.spatial_slices import (
    assign_slices_from_radar_time_rank,
    filter_pairs_by_same_slice,
    radar_distance_rank_from_tick_brackets,
    radar_time_rank_from_tick_brackets,
    same_slice_unordered_pairs,
)

DEFAULT_OUTPUT = Path("outputs/data/spatial_slice_exact_sanity.csv")


def _row(check: str, passed: bool, value: float) -> dict[str, float | str]:
    return {"check": check, "passed": float(passed), "value": float(value)}


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic spatial-slice checks."""

    predecessor = np.asarray([1, 2, -1, 4, 4])
    successor = np.asarray([3, 6, -1, 7, 6])
    accessible = np.asarray([True, True, False, True, True])
    time_rank = radar_time_rank_from_tick_brackets(predecessor, successor, accessible)
    distance_rank = radar_distance_rank_from_tick_brackets(
        predecessor,
        successor,
        accessible,
    )
    labels = assign_slices_from_radar_time_rank(time_rank, accessible, bin_width=4)
    pairs = same_slice_unordered_pairs(labels)
    filtered = filter_pairs_by_same_slice(np.asarray([[0, 1], [1, 3], [1, 2]]), labels)
    pair_labels = pair_slice_labels(np.asarray([[0, 1], [1, 3], [1, 2]]), labels)
    coords = np.asarray([0.0, 1.0, 3.0, 4.0, 6.0])
    translated = coords + 10.0
    inversion = sliced_pair_distance_order_inversion_rate(coords, translated, pairs)
    pair_distances = np.asarray([1.0, 3.0, 2.0, 5.0])
    constraint_pairs = np.asarray([[0, 1], [0, 3], [1, 3], [2, 3]])
    constraint_labels = np.asarray([0, 0, 1, 1])
    constraints = quadruplet_constraints_from_sliced_pair_distances(
        constraint_pairs,
        pair_distances,
        constraint_labels,
        4,
        seed=0,
    )

    return [
        _row(
            "radar_time_rank",
            np.array_equal(time_rank, [4, 8, -1, 11, 10]),
            time_rank[0],
        ),
        _row(
            "radar_distance_rank",
            np.array_equal(distance_rank, [2, 4, -1, 3, 2]),
            distance_rank[0],
        ),
        _row("assign_slices", np.array_equal(labels, [1, 2, -1, 2, 2]), labels[0]),
        _row("same_slice_pairs", pairs.shape[1] == 2, pairs.shape[0]),
        _row(
            "filter_same_slice",
            np.array_equal(filtered, [[1, 3]]),
            filtered.shape[0],
        ),
        _row(
            "pair_slice_labels",
            np.array_equal(pair_labels, [-1, 2, -1]),
            pair_labels[1],
        ),
        _row("translated_distance_order", inversion == 0.0, inversion),
        _row(
            "sliced_constraints_valid",
            constraints.ndim == 2 and constraints.shape[1] == 4,
            constraints.shape[0],
        ),
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write exact sanity rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote spatial slice exact sanity results: {output_path}")


if __name__ == "__main__":
    main()
