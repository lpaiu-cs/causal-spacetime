"""Exact sanity checks for representation-stability utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.embedding_stability import (
    pairwise_order_stability,
    pairwise_procrustes_stability,
)
from causal_spacetime_lab.ordinal_embedding import (
    sample_quadruplet_constraints_from_coords,
)
from causal_spacetime_lab.representation_validation import (
    evaluate_embedding_on_constraints,
    flip_constraint_fraction,
    random_quadruplet_constraints,
    shuffle_constraint_sides,
    split_constraints,
)

DEFAULT_OUTPUT = Path("outputs/data/representation_stability_exact_sanity.csv")


def _row(check: str, passed: bool, value: float) -> dict[str, float | str]:
    return {"check": check, "passed": float(passed), "value": float(value)}


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic exact checks."""

    coords = np.asarray([[0.0], [1.0], [3.0], [6.0], [10.0]])
    constraints = sample_quadruplet_constraints_from_coords(coords, 40, seed=1)
    train, test = split_constraints(constraints, train_fraction=0.75, seed=2)
    shuffled = shuffle_constraint_sides(constraints, seed=3)
    random_constraints = random_quadruplet_constraints(5, 20, seed=4)
    no_flip = flip_constraint_fraction(constraints, 0.0, seed=5)
    full_flip = flip_constraint_fraction(constraints, 1.0, seed=5)
    exact_eval = evaluate_embedding_on_constraints(coords, constraints)
    identical_embeddings = [coords, coords.copy(), coords.copy()]
    procrustes = pairwise_procrustes_stability(identical_embeddings)
    order = pairwise_order_stability(identical_embeddings, seed=6)

    return [
        _row(
            "split_constraints_preserves_rows",
            train.shape[0] + test.shape[0] == constraints.shape[0],
            train.shape[0] + test.shape[0],
        ),
        _row(
            "split_constraints_fraction",
            abs(train.shape[0] / constraints.shape[0] - 0.75) < 0.1,
            train.shape[0] / constraints.shape[0],
        ),
        _row(
            "shuffle_preserves_shape",
            shuffled.shape == constraints.shape,
            shuffled.shape[0],
        ),
        _row(
            "random_constraints_valid",
            bool(
                np.all(random_constraints[:, 0] != random_constraints[:, 1])
                and np.all(random_constraints[:, 2] != random_constraints[:, 3])
            ),
            random_constraints.shape[0],
        ),
        _row(
            "flip_zero_returns_original",
            np.array_equal(no_flip, constraints),
            float(np.mean(no_flip == constraints)),
        ),
        _row(
            "flip_one_swaps_all",
            bool(
                np.array_equal(full_flip[:, 0:2], constraints[:, 2:4])
                and np.array_equal(full_flip[:, 2:4], constraints[:, 0:2])
            ),
            full_flip.shape[0],
        ),
        _row(
            "exact_constraints_have_zero_violation",
            exact_eval["violation_rate"] == 0.0,
            exact_eval["violation_rate"],
        ),
        _row(
            "identical_procrustes_stability_zero",
            procrustes["mean_procrustes_rmse"] < 1e-12,
            procrustes["mean_procrustes_rmse"],
        ),
        _row(
            "identical_order_stability_zero",
            order["mean_order_disagreement"] < 1e-12,
            order["mean_order_disagreement"],
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
    print(f"Wrote representation stability exact sanity results: {output_path}")


if __name__ == "__main__":
    main()
