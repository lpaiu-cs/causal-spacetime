"""Exact sanity checks for ordinal and metric-representation utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.metric_representation import (
    distance_order_constraints_from_values,
    euclidean_embedding_diagnostics,
    has_order_cycle,
    topological_rank_representation,
    triangle_inequality_violations,
    unordered_pair_indices,
)
from causal_spacetime_lab.ordinal import (
    order_inversion_rate,
    order_matrix_from_values,
    pair_distance_order_inversion_rate,
    pair_distance_values_1d,
)
from causal_spacetime_lab.radar_order import (
    earliest_successor_tick_index,
    latest_predecessor_tick_index,
    radar_return_order_from_successor_ticks,
    radar_tick_brackets_from_order,
)

OUTPUT_DATA = Path("outputs/data/ordinal_exact_sanity.csv")


def _row(
    check: str,
    passed: bool,
    observed: float,
    expected: float,
) -> dict[str, float | str]:
    return {
        "check": check,
        "passed": float(passed),
        "observed": observed,
        "expected": expected,
        "absolute_error": abs(observed - expected),
    }


def run_experiment() -> list[dict[str, float | str]]:
    """Return exact ordinal sanity checks."""

    rows: list[dict[str, float | str]] = []
    values = np.asarray([2.0, 1.0, 3.0])
    strict = order_matrix_from_values(values, strict=True)
    nonstrict = order_matrix_from_values(values, strict=False)
    rows.append(
        _row("strict_order_count", int(np.count_nonzero(strict)) == 3, 3.0, 3.0)
    )
    rows.append(
        _row(
            "nonstrict_order_count",
            int(np.count_nonzero(nonstrict)) == 6,
            6.0,
            6.0,
        )
    )
    base = np.asarray([1.0, 2.0, 4.0])
    rows.append(
        _row(
            "identical_inversion",
            order_inversion_rate(base, base) == 0.0,
            0.0,
            0.0,
        )
    )
    scaled_inversion = order_inversion_rate(base, 3.0 * base)
    rows.append(
        _row(
            "positive_scaling_inversion",
            scaled_inversion == 0.0,
            scaled_inversion,
            0.0,
        )
    )
    reversed_inversion = order_inversion_rate(base, base[::-1])
    rows.append(
        _row(
            "reversed_inversion",
            reversed_inversion == 1.0,
            reversed_inversion,
            1.0,
        )
    )

    coords = np.asarray([0.0, 2.0, 5.0])
    pairs = np.asarray([[0, 1], [1, 2], [0, 2]], dtype=int)
    distances = pair_distance_values_1d(coords, pairs)
    rows.append(
        _row(
            "pair_distance_sum",
            np.sum(distances) == 10.0,
            float(np.sum(distances)),
            10.0,
        )
    )
    distance_inversion = pair_distance_order_inversion_rate(coords, coords + 1.0, pairs)
    rows.append(
        _row(
            "pair_distance_order_translation",
            distance_inversion == 0.0,
            distance_inversion,
            0.0,
        )
    )

    causal_matrix = np.zeros((6, 6), dtype=bool)
    causal_matrix[1, 4] = True
    causal_matrix[2, 4] = True
    causal_matrix[4, 3] = True
    causal_matrix[0, 5] = True
    causal_matrix[5, 3] = True
    observers = np.asarray([0, 1, 2, 3], dtype=int)
    rows.append(
        _row(
            "latest_predecessor_tick",
            latest_predecessor_tick_index(causal_matrix, observers, 4) == 2,
            2.0,
            2.0,
        )
    )
    rows.append(
        _row(
            "earliest_successor_tick",
            earliest_successor_tick_index(causal_matrix, observers, 4) == 3,
            3.0,
            3.0,
        )
    )
    _, successors, accessible = radar_tick_brackets_from_order(
        causal_matrix,
        observers,
        np.asarray([4, 5], dtype=int),
    )
    order = radar_return_order_from_successor_ticks(successors, accessible)
    rows.append(
        _row(
            "radar_bracket_accessible_count",
            np.count_nonzero(accessible) == 2,
            2.0,
            2.0,
        )
    )
    rows.append(
        _row(
            "radar_return_order_count",
            np.count_nonzero(order) == 0,
            0.0,
            0.0,
        )
    )

    unordered = unordered_pair_indices(3)
    rows.append(_row("unordered_pair_count", unordered.shape[0] == 3, 3.0, 3.0))
    constraints = distance_order_constraints_from_values(np.asarray([1.0, 2.0, 3.0]))
    rows.append(
        _row("has_no_order_cycle", not has_order_cycle(3, constraints), 0.0, 0.0)
    )
    ranks = topological_rank_representation(3, constraints)
    rows.append(
        _row(
            "rank_representation_span",
            np.max(ranks) == 2.0,
            float(np.max(ranks)),
            2.0,
        )
    )
    rows.append(
        _row(
            "triangle_violation_count",
            triangle_inequality_violations(np.asarray([1.0, 1.0, 3.0]), 3) == 1,
            1.0,
            1.0,
        )
    )
    matrix = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=float)
    diagnostics = euclidean_embedding_diagnostics(matrix)
    rows.append(
        _row(
            "euclidean_candidate",
            bool(diagnostics["is_euclidean_candidate"]),
            1.0,
            1.0,
        )
    )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write exact ordinal sanity rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote ordinal exact sanity results: {output_path}")


if __name__ == "__main__":
    main()
