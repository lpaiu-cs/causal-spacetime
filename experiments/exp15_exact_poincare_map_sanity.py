"""Exact-coordinate sanity checks for affine Lorentz/Poincare map fitting."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    affine_lab_to_rest_1p1,
)
from causal_spacetime_lab.poincare_maps import (
    compose_affine_lorentz_maps_1p1,
    expected_relative_beta_1p1,
    fit_affine_lorentz_beta_grid_1p1,
)

OUTPUT_DATA = Path("outputs/data/exact_poincare_map_sanity.csv")


def default_protocols() -> tuple[ObserverProtocolSpec, ...]:
    """Return three deterministic affine observer protocols."""

    return (
        ObserverProtocolSpec("A_lab", 0.0, 0.0, 0.0, 0.15),
        ObserverProtocolSpec("B_moving_pos", 0.3, 0.05, -0.05, 0.15),
        ObserverProtocolSpec("C_moving_neg", -0.25, -0.04, 0.04, 0.15),
    )


def _random_lab_events(seed: int = 0, count: int = 200) -> np.ndarray:
    rng = np.random.default_rng(seed)
    t = rng.uniform(-0.5, 0.5, size=count)
    x = rng.uniform(-0.45, 0.45, size=count)
    return np.column_stack((t, x)).astype(float, copy=False)


def run_experiment() -> list[dict[str, float | str]]:
    """Fit exact affine maps between hidden protocol coordinates."""

    protocols = default_protocols()
    lab_events = _random_lab_events()
    coords = {
        spec.name: affine_lab_to_rest_1p1(lab_events, spec.beta, spec.origin_lab)
        for spec in protocols
    }
    specs = {spec.name: spec for spec in protocols}
    transitions = [("A_lab", "B_moving_pos"), ("B_moving_pos", "C_moving_neg")]
    transitions.append(("A_lab", "C_moving_neg"))
    rows: list[dict[str, float | str]] = []
    fit_results: dict[tuple[str, str], tuple[float, np.ndarray, float]] = {}

    for source_name, target_name in transitions:
        source = specs[source_name]
        target = specs[target_name]
        expected_beta = expected_relative_beta_1p1(source.beta, target.beta)
        fitted_beta, translation, residual = fit_affine_lorentz_beta_grid_1p1(
            coords[source_name],
            coords[target_name],
            beta_min=-0.95,
            beta_max=0.95,
            num_grid=1901,
        )
        fit_results[(source_name, target_name)] = (
            fitted_beta,
            translation,
            residual,
        )
        rows.append(
            {
                "kind": "transition",
                "source_protocol": source_name,
                "target_protocol": target_name,
                "expected_beta": expected_beta,
                "fitted_beta": fitted_beta,
                "fitted_beta_error": fitted_beta - expected_beta,
                "translation_t": float(translation[0]),
                "translation_x": float(translation[1]),
                "rmse": residual,
                "beta_composition_error": float("nan"),
                "translation_composition_error_norm": float("nan"),
            }
        )

    beta_ab, trans_ab, _ = fit_results[("A_lab", "B_moving_pos")]
    beta_bc, trans_bc, _ = fit_results[("B_moving_pos", "C_moving_neg")]
    beta_ac, trans_ac, _ = fit_results[("A_lab", "C_moving_neg")]
    beta_composed, trans_composed = compose_affine_lorentz_maps_1p1(
        beta_ab,
        trans_ab,
        beta_bc,
        trans_bc,
    )
    rows.append(
        {
            "kind": "loop",
            "source_protocol": "A_lab",
            "target_protocol": "C_moving_neg",
            "expected_beta": expected_relative_beta_1p1(0.0, -0.25),
            "fitted_beta": beta_ac,
            "fitted_beta_error": beta_ac - expected_relative_beta_1p1(0.0, -0.25),
            "translation_t": float(trans_ac[0]),
            "translation_x": float(trans_ac[1]),
            "rmse": 0.0,
            "beta_composition_error": beta_composed - beta_ac,
            "translation_composition_error_norm": float(
                np.linalg.norm(trans_composed - trans_ac)
            ),
        }
    )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write exact-coordinate Poincare sanity results."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    print(f"Wrote exact Poincare sanity results: {output_path}")


if __name__ == "__main__":
    main()
