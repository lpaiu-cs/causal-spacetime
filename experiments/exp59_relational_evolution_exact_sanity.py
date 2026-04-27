"""Exact sanity checks for relational spatial evolution utilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.predicate_definability import predicate_requirement_table
from causal_spacetime_lab.relational_evolution import (
    apply_per_slice_affine_position_gauge,
    compare_histories_order_disagreement,
    pair_distance_order_signature,
    pair_order_comparison_matrix,
    relational_change_rate_between_slices,
    relational_shape_history,
    signature_order_disagreement,
    unordered_object_pairs,
)
from causal_spacetime_lab.relational_scenarios import (
    generate_expanding_configuration_history,
    generate_shear_or_reordering_history_1d,
    generate_static_configuration_history,
)

DEFAULT_OUTPUT = Path("outputs/data/relational_evolution_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Relational evolution sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "relational_evolution_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact sanity checks."""

    pairs = unordered_object_pairs(np.asarray([2, 0, 1]))
    signature = pair_distance_order_signature({0: 0.0, 1: 1.0, 2: 3.0})
    matrix = pair_order_comparison_matrix(signature)
    static = generate_static_configuration_history(
        np.asarray([0, 1, 2]),
        np.asarray([0.0, 1.0, 3.0]),
    )
    expansion = generate_expanding_configuration_history(
        np.asarray([0, 1, 2]),
        np.asarray([0.0, 1.0, 3.0]),
        np.asarray([1.0, 2.0, 3.0]),
    )
    reordering = generate_shear_or_reordering_history_1d(
        np.asarray([0, 1, 2]),
        np.asarray([0.0, 1.0, 3.0]),
        moving_object_id=1,
        displacement_by_slice=np.asarray([0.0, 1.5, 3.0]),
    )
    static_history = relational_shape_history(static)
    expansion_history = relational_shape_history(expansion)
    reorder_history = relational_shape_history(reordering)
    static_change = relational_change_rate_between_slices(static_history)
    reorder_change = relational_change_rate_between_slices(reorder_history)
    gauged = apply_per_slice_affine_position_gauge(reordering, seed=3)
    table = predicate_requirement_table()
    table_names = {str(row["name"]) for row in table}
    return [
        _row(
            "unordered_object_pairs",
            np.array_equal(pairs, [[0, 1], [0, 2], [1, 2]]),
            pairs.shape[0],
        ),
        _row(
            "signature_ranks",
            np.array_equal(signature.pair_order_ranks, [0.0, 2.0, 1.0]),
            signature.pair_order_ranks[0],
        ),
        _row(
            "comparison_matrix",
            matrix[0, 1] == -1 and matrix[1, 0] == 1,
            matrix[0, 1],
        ),
        _row(
            "identical_disagreement",
            signature_order_disagreement(signature, signature) == 0.0,
            0.0,
        ),
        _row(
            "uniform_scaling_preserves_order",
            compare_histories_order_disagreement(
                static_history,
                expansion_history,
            )
            == 0.0,
            0.0,
        ),
        _row(
            "static_change_zero",
            all(row["order_disagreement"] == 0.0 for row in static_change),
            0.0,
        ),
        _row(
            "reordering_change_nonzero",
            any(row["order_disagreement"] > 0.0 for row in reorder_change),
            reorder_change[-1]["order_disagreement"],
        ),
        _row(
            "affine_gauge_preserves_history",
            compare_histories_order_disagreement(
                reorder_history,
                relational_shape_history(gauged),
            )
            == 0.0,
            0.0,
        ),
        _row(
            "predicate_table_entries",
            {"coordinate_velocity", "pair_distance_order_history"} <= table_names,
            len(table_names),
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
    output_path = parse_args()
    rows = run_experiment()
    output_path = write_outputs(rows, output_path)
    print(f"Wrote relational evolution exact sanity data: {output_path}")


if __name__ == "__main__":
    main()
