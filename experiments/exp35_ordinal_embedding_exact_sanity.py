"""Exact sanity checks for ordinal embedding and calibration utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.calibration import (
    affine_fit_1d,
    equal_step_deviation,
    ratio_error_under_transform,
)
from causal_spacetime_lab.ordinal_embedding import (
    procrustes_align,
    quadruplet_hinge_loss,
    quadruplet_violation_rate,
    sample_quadruplet_constraints_from_coords,
    squared_distance_matrix,
)

OUTPUT_DATA = Path("outputs/data/ordinal_embedding_exact_sanity.csv")


def _row(
    check: str,
    observed: float,
    expected: float,
    passed: bool,
) -> dict[str, float | str]:
    return {
        "check": check,
        "observed": observed,
        "expected": expected,
        "absolute_error": abs(observed - expected),
        "passed": float(passed),
    }


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic sanity checks."""

    rows: list[dict[str, float | str]] = []
    points = np.asarray([[0.0, 0.0], [3.0, 4.0], [1.0, 0.0]], dtype=float)
    distances = squared_distance_matrix(points)
    rows.append(_row("squared_distance_3_4", float(distances[0, 1]), 25.0, True))

    constraints = sample_quadruplet_constraints_from_coords(points, 10, seed=3)
    violation = quadruplet_violation_rate(points, constraints)
    rows.append(_row("exact_violation_rate", violation, 0.0, violation == 0.0))
    loss = quadruplet_hinge_loss(points, constraints, margin=0.0)
    rows.append(_row("exact_hinge_loss", loss, 0.0, loss < 1e-12))

    theta = np.pi / 5.0
    rotation = np.asarray(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )
    transformed = 2.5 * points @ rotation + np.asarray([4.0, -2.0])
    _, procrustes = procrustes_align(transformed, points)
    rows.append(
        _row("procrustes_rmse", procrustes["rmse"], 0.0, procrustes["rmse"] < 1e-12)
    )

    source = np.arange(5.0)
    target = 2.0 * source + 3.0
    scale, offset, rmse = affine_fit_1d(source, target)
    rows.append(_row("affine_fit_scale", scale, 2.0, abs(scale - 2.0) < 1e-12))
    rows.append(_row("affine_fit_offset", offset, 3.0, abs(offset - 3.0) < 1e-12))
    rows.append(_row("affine_fit_rmse", rmse, 0.0, rmse < 1e-12))

    step_deviation = equal_step_deviation(source)
    rows.append(
        _row(
            "equal_step_deviation",
            step_deviation,
            0.0,
            step_deviation < 1e-12,
        )
    )
    intervals = np.asarray([[0, 2, 0, 4], [1, 3, 0, 2]], dtype=int)
    affine_ratio_error = ratio_error_under_transform(source, target, intervals)
    square_ratio_error = ratio_error_under_transform(source, source * source, intervals)
    rows.append(
        _row(
            "affine_ratio_error",
            affine_ratio_error,
            0.0,
            affine_ratio_error < 1e-12,
        )
    )
    rows.append(
        _row(
            "square_ratio_error_positive",
            square_ratio_error,
            0.0,
            square_ratio_error > 0.0,
        )
    )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
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
    print(f"Wrote ordinal embedding exact sanity results: {output_path}")


if __name__ == "__main__":
    main()
