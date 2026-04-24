"""Timelike proper-time reconstruction in a 1+1D causal diamond."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.chains import longest_chain_length
from causal_spacetime_lab.intervals import alexandrov_interval_indices
from causal_spacetime_lab.metrics import (
    causal_diamond_volume_1p1,
    estimate_tau_from_interval_count,
    minkowski_tau_1p1,
)
from causal_spacetime_lab.plotting import save_timelike_reconstruction_plot
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

EVENT_COUNTS = (200, 500, 1000, 2000)
T = 2.0
BASE_SEED = 20260424
OUTPUT_DATA = Path("outputs/data/timelike_reconstruction_summary.csv")
OUTPUT_FIGURE = Path("outputs/figures/timelike_reconstruction_error.png")


def run_experiment() -> list[dict[str, float]]:
    """Run the milestone timelike reconstruction experiment."""

    p = np.array([-T / 2.0, 0.0])
    q = np.array([T / 2.0, 0.0])
    true_tau = minkowski_tau_1p1(p, q)
    diamond_volume = causal_diamond_volume_1p1(T)

    rows: list[dict[str, float]] = []

    for offset, n_events in enumerate(EVENT_COUNTS):
        interior_events = sprinkle_1p1_causal_diamond(
            n_events,
            T=T,
            seed=BASE_SEED + offset,
        )
        events = np.vstack((p, interior_events, q))
        start_index = 0
        end_index = events.shape[0] - 1

        causal_matrix = causal_matrix_1p1(events)
        chain_length = longest_chain_length(
            causal_matrix,
            start=start_index,
            end=end_index,
            event_times=events[:, 0],
        )
        interval_indices = alexandrov_interval_indices(
            causal_matrix,
            start_index,
            end_index,
        )
        interval_count = int(interval_indices.size)

        rho = n_events / diamond_volume
        interval_tau_estimate = estimate_tau_from_interval_count(interval_count, rho)

        chain_interior_length = max(chain_length - 2, 0)
        chain_tau_estimate = chain_interior_length / np.sqrt(2.0 * rho)

        rows.append(
            {
                "n_events": float(n_events),
                "T": T,
                "diamond_volume": diamond_volume,
                "rho": rho,
                "true_tau": true_tau,
                "longest_chain_length": float(chain_length),
                "chain_tau_estimate": float(chain_tau_estimate),
                "chain_tau_abs_error": abs(float(chain_tau_estimate) - true_tau),
                "interval_count": float(interval_count),
                "interval_tau_estimate": float(interval_tau_estimate),
                "interval_tau_abs_error": abs(float(interval_tau_estimate) - true_tau),
            }
        )

    return rows


def write_summary(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write experiment rows to CSV."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    rows = run_experiment()
    summary_path = write_summary(rows)
    figure_path = save_timelike_reconstruction_plot(rows, OUTPUT_FIGURE)

    print(f"Wrote summary: {summary_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
