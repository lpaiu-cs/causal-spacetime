"""Exact sanity checks for persistence matching utilities."""

from __future__ import annotations

import argparse
import csv
from math import factorial
from pathlib import Path

import numpy as np

from causal_spacetime_lab.identity_constraints import (
    filter_permutations_by_fixed_points,
    generate_partial_identity_constraints,
)
from causal_spacetime_lab.persistence_definability import (
    object_identity_without_persistence,
)
from causal_spacetime_lab.persistence_history import (
    build_identity_tracks_from_adjacent_matchings,
    infer_adjacent_matchings,
    relational_history_from_inferred_tracks,
)
from causal_spacetime_lab.persistence_matching import (
    all_permutations,
    apply_permutation_to_positions,
    matching_accuracy,
    matching_ambiguity_gap,
    signature_disagreement_for_permutation,
)
from causal_spacetime_lab.persistence_scenarios import (
    generate_unlabeled_static_history,
)

DEFAULT_OUTPUT = Path("outputs/data/persistence_matching_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Persistence matching sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "persistence_matching_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact persistence matching checks."""

    permutations = all_permutations(4)
    positions = np.asarray([10.0, 20.0, 30.0])
    permuted = apply_permutation_to_positions(positions, np.asarray([2, 0, 1]))
    matches = [
        {"cost": 0.0, "permutation": np.asarray([0, 1, 2])},
        {"cost": 0.25, "permutation": np.asarray([1, 0, 2])},
    ]
    filtered = filter_permutations_by_fixed_points(
        all_permutations(3),
        {0: 2},
    )
    fixed = generate_partial_identity_constraints(np.asarray([2, 0, 1]), 2 / 3, seed=1)
    unlabeled, _tracks = generate_unlabeled_static_history(
        2,
        np.asarray([0.0, 1.0, 3.0]),
        random_per_slice_permutation=False,
    )
    inferred = infer_adjacent_matchings(unlabeled, top_k=1)
    best = {
        edge: np.asarray(rows[0]["permutation"], dtype=int)
        for edge, rows in inferred.items()
    }
    built_tracks = build_identity_tracks_from_adjacent_matchings(
        best,
        np.asarray([0, 1, 2]),
    )
    history = relational_history_from_inferred_tracks(unlabeled, built_tracks)
    undefined = object_identity_without_persistence()
    return [
        _row(
            "all_permutations_count",
            permutations.shape[0] == factorial(4),
            permutations.shape[0],
        ),
        _row(
            "apply_permutation_to_positions",
            np.array_equal(permuted, [30.0, 10.0, 20.0]),
            permuted[0],
        ),
        _row(
            "signature_disagreement_identity",
            signature_disagreement_for_permutation(
                positions,
                positions,
                np.asarray([0, 1, 2]),
            )
            == 0.0,
            0.0,
        ),
        _row("matching_ambiguity_gap", matching_ambiguity_gap(matches) == 0.25, 0.25),
        _row(
            "matching_accuracy",
            matching_accuracy(np.asarray([0, 2, 1]), np.asarray([0, 1, 1]))
            == 2 / 3,
            2 / 3,
        ),
        _row("filter_fixed_points", np.all(filtered[:, 0] == 2), filtered.shape[0]),
        _row("partial_identity_constraints", len(fixed) == 2, len(fixed)),
        _row("infer_adjacent_matchings", bool(inferred), len(inferred)),
        _row("build_identity_tracks", 1 in built_tracks, len(built_tracks)),
        _row("relational_history_from_tracks", set(history) == {0, 1}, len(history)),
        _row("object_identity_undefined", not undefined.defined, undefined.reason),
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
    print(f"Wrote persistence matching exact sanity data: {output_path}")


if __name__ == "__main__":
    main()
