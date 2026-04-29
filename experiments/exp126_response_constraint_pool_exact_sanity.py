"""Exact sanity checks for response-comparison constraint pools."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
    evaluate_constraint_pool_on_dissimilarity,
    filter_constraint_pool_by_margin,
    merge_constraint_pools,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseDissimilarity,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def _dissimilarity(values: list[float]) -> PairwiseResponseDissimilarity:
    return PairwiseResponseDissimilarity(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        pair_indices=np.asarray(
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            dtype=int,
        ),
        dissimilarity_values=np.asarray(values, dtype=float),
        valid_pair_mask=np.isfinite(np.asarray(values, dtype=float)),
        protocol_name="gap",
        method="rank_gap_mean",
    )


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact response-comparison pool sanity checks."""

    baseline = _dissimilarity([0.1, 0.2, 0.4, 0.5, 0.7, 0.9])
    inverted = _dissimilarity([0.9, 0.7, 0.5, 0.4, 0.2, 0.1])
    pool = build_constraint_pool_from_dissimilarity(
        baseline,
        max_constraints=100,
        min_margin=0.1,
        seed=0,
    )
    identical = evaluate_constraint_pool_on_dissimilarity(pool, baseline)
    reversed_report = evaluate_constraint_pool_on_dissimilarity(pool, inverted)
    filtered = filter_constraint_pool_by_margin(pool, 0.5)
    merged = merge_constraint_pools([filtered, filtered], source_label="merged")
    checks = [
        ("pool_nonempty", pool.constraints.shape[0] > 0, pool.constraints.shape[0]),
        (
            "identical_agreement",
            identical["agreement_fraction"] == 1.0,
            identical["agreement_fraction"],
        ),
        (
            "inverted_disagreement",
            reversed_report["inversion_fraction"] == 1.0,
            reversed_report["inversion_fraction"],
        ),
        (
            "margin_filter_reduces_or_keeps",
            filtered.constraints.shape[0] <= pool.constraints.shape[0],
            filtered.constraints.shape[0],
        ),
        (
            "merge_deduplicates",
            merged.constraints.shape[0] == filtered.constraints.shape[0],
            merged.constraints.shape[0],
        ),
    ]
    return [
        {"check": name, "passed": float(passed), "value": str(value)}
        for name, passed, value in checks
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact sanity CSV."""

    path = output_dir / "data" / "response_constraint_pool_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote response constraint pool exact sanity: {output_path}")


if __name__ == "__main__":
    main()
