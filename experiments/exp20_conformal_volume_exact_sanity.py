"""Deterministic sanity checks for conformal volume utilities."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.conformal import (
    central_observer_proper_time_1p1,
    constant_profile,
    flat_profile,
    integrate_conformal_interval_volume_1p1,
)

OUTPUT_DATA = Path("outputs/data/conformal_volume_exact_sanity.csv")


def run_experiment(T: float = 2.0, scale: float = 1.5) -> list[dict[str, float | str]]:
    """Return deterministic conformal volume sanity rows."""

    p = np.asarray([-T / 2.0, 0.0])
    q = np.asarray([T / 2.0, 0.0])
    flat_volume_expected = T * T / 2.0
    flat_volume = integrate_conformal_interval_volume_1p1(p, q, flat_profile())
    scaled_profile = constant_profile(scale)
    scaled_volume = integrate_conformal_interval_volume_1p1(p, q, scaled_profile)
    scaled_time = central_observer_proper_time_1p1(
        -T / 2.0,
        T / 2.0,
        scaled_profile,
    )
    return [
        {
            "check": "flat_full_diamond_volume",
            "expected": flat_volume_expected,
            "observed": flat_volume,
            "absolute_error": abs(flat_volume - flat_volume_expected),
        },
        {
            "check": "constant_scale_full_diamond_volume",
            "expected": scale * scale * flat_volume_expected,
            "observed": scaled_volume,
            "absolute_error": abs(scaled_volume - scale * scale * flat_volume_expected),
        },
        {
            "check": "constant_scale_central_proper_time",
            "expected": scale * T,
            "observed": scaled_time,
            "absolute_error": abs(scaled_time - scale * T),
        },
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write exact conformal sanity rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote exact conformal sanity results: {output_path}")


if __name__ == "__main__":
    main()
