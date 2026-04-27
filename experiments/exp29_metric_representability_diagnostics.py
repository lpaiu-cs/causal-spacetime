"""Finite diagnostics for distance-order and metric representability."""

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
from causal_spacetime_lab.ordinal import pair_distance_values_1d

OUTPUT_DATA = Path("outputs/data/metric_representability_diagnostics.csv")


def _distance_matrix_from_values(values: np.ndarray, n_points: int) -> np.ndarray:
    matrix = np.zeros((n_points, n_points), dtype=float)
    pairs = unordered_pair_indices(n_points)
    matrix[pairs[:, 0], pairs[:, 1]] = values
    matrix[pairs[:, 1], pairs[:, 0]] = values
    return matrix


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic representability diagnostics."""

    rows: list[dict[str, float | str]] = []

    coords = np.asarray([0.0, 1.0, 3.0], dtype=float)
    pairs = unordered_pair_indices(coords.size)
    distances = pair_distance_values_1d(coords, pairs)
    constraints = distance_order_constraints_from_values(distances)
    ranks = topological_rank_representation(distances.size, constraints)
    diagnostics = euclidean_embedding_diagnostics(
        _distance_matrix_from_values(distances, coords.size)
    )
    rows.append(
        {
            "case": "consistent_1d_metric",
            "has_order_cycle": float(has_order_cycle(distances.size, constraints)),
            "triangle_violations": float(
                triangle_inequality_violations(distances, coords.size)
            ),
            "min_gram_eigenvalue": diagnostics["min_gram_eigenvalue"],
            "positive_eigenvalue_count": diagnostics["positive_eigenvalue_count"],
            "is_euclidean_candidate": float(diagnostics["is_euclidean_candidate"]),
            "rank_span": float(np.max(ranks) - np.min(ranks)),
        }
    )

    cyclic_constraints = np.asarray([[0, 1], [1, 2], [2, 0]], dtype=int)
    rows.append(
        {
            "case": "contradictory_order_cycle",
            "has_order_cycle": float(has_order_cycle(3, cyclic_constraints)),
            "triangle_violations": float("nan"),
            "min_gram_eigenvalue": float("nan"),
            "positive_eigenvalue_count": float("nan"),
            "is_euclidean_candidate": 0.0,
            "rank_span": float("nan"),
        }
    )

    nonmetric_distances = np.asarray([1.0, 3.0, 1.0], dtype=float)
    nonmetric_matrix = _distance_matrix_from_values(nonmetric_distances, 3)
    nonmetric_diagnostics = euclidean_embedding_diagnostics(nonmetric_matrix)
    rows.append(
        {
            "case": "non_euclidean_candidate",
            "has_order_cycle": float(
                has_order_cycle(
                    nonmetric_distances.size,
                    distance_order_constraints_from_values(nonmetric_distances),
                )
            ),
            "triangle_violations": float(
                triangle_inequality_violations(nonmetric_distances, 3)
            ),
            "min_gram_eigenvalue": nonmetric_diagnostics["min_gram_eigenvalue"],
            "positive_eigenvalue_count": nonmetric_diagnostics[
                "positive_eigenvalue_count"
            ],
            "is_euclidean_candidate": float(
                nonmetric_diagnostics["is_euclidean_candidate"]
            ),
            "rank_span": float("nan"),
        }
    )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write metric representability diagnostics."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote metric representability diagnostics: {output_path}")


if __name__ == "__main__":
    main()
