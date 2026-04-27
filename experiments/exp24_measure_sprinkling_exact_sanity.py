"""Deterministic sanity checks for measure-sprinkling utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.conformal import constant_profile, flat_profile
from causal_spacetime_lab.measure_sprinkling import (
    conformal_full_diamond_volume_1p1,
    conformal_time_bin_masses_1p1,
    coordinate_time_bin_volumes_1p1,
    normalize_profile_shape,
)

OUTPUT_DATA = Path("outputs/data/measure_sprinkling_exact_sanity.csv")


def run_experiment(T: float = 2.0, scale: float = 1.5) -> list[dict[str, float | str]]:
    """Return exact sanity checks for measure-sprinkling utilities."""

    edges = np.linspace(-T / 2.0, T / 2.0, 9)
    flat_expected = T * T / 2.0
    flat_volume = conformal_full_diamond_volume_1p1(T, flat_profile())
    constant_volume = conformal_full_diamond_volume_1p1(T, constant_profile(scale))
    coordinate_bins = coordinate_time_bin_volumes_1p1(T, edges)
    conformal_bins = conformal_time_bin_masses_1p1(T, constant_profile(scale), edges)
    normalized = normalize_profile_shape(np.asarray([2.0, 4.0, 6.0]))
    return [
        {
            "check": "flat_full_diamond_volume",
            "expected": flat_expected,
            "observed": flat_volume,
            "absolute_error": abs(flat_volume - flat_expected),
        },
        {
            "check": "constant_full_diamond_volume",
            "expected": scale * scale * flat_expected,
            "observed": constant_volume,
            "absolute_error": abs(constant_volume - scale * scale * flat_expected),
        },
        {
            "check": "coordinate_bin_volume_sum",
            "expected": flat_expected,
            "observed": float(np.sum(coordinate_bins)),
            "absolute_error": abs(float(np.sum(coordinate_bins)) - flat_expected),
        },
        {
            "check": "conformal_bin_mass_sum",
            "expected": constant_volume,
            "observed": float(np.sum(conformal_bins)),
            "absolute_error": abs(float(np.sum(conformal_bins)) - constant_volume),
        },
        {
            "check": "normalize_profile_shape_mean",
            "expected": 1.0,
            "observed": float(np.nanmean(normalized)),
            "absolute_error": abs(float(np.nanmean(normalized)) - 1.0),
        },
    ]


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
    print(f"Wrote measure-sprinkling exact sanity results: {output_path}")


if __name__ == "__main__":
    main()
