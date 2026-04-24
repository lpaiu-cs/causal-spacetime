"""Lorentz length-contraction experiment from lab-simultaneous event selection."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.lorentz import gamma, measured_rod_length_lab

L0 = 1.0
BETAS = np.linspace(0.0, 0.95, 40)
OUTPUT_DATA = Path("outputs/data/lorentz_length_contraction_summary.csv")
OUTPUT_FIGURE = Path("outputs/figures/lorentz_length_contraction.png")


def run_experiment() -> list[dict[str, float]]:
    """Compute lab-measured rod length over a range of relative velocities."""

    rows: list[dict[str, float]] = []
    for beta in BETAS:
        measured_length = measured_rod_length_lab(L0, float(beta))
        expected_length = L0 / float(gamma(beta))
        rows.append(
            {
                "beta": float(beta),
                "gamma": float(gamma(beta)),
                "L0": L0,
                "measured_length_lab": measured_length,
                "expected_length": expected_length,
                "abs_error": abs(measured_length - expected_length),
            }
        )
    return rows


def write_summary(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write experiment rows to CSV."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(rows: list[dict[str, float]], output_path: Path = OUTPUT_FIGURE) -> Path:
    """Save measured length versus beta."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    betas = np.asarray([row["beta"] for row in rows])
    measured_lengths = np.asarray([row["measured_length_lab"] for row in rows])
    expected_lengths = np.asarray([row["expected_length"] for row in rows])

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(betas, measured_lengths, marker="o", markersize=3, label="Measured")
    ax.plot(betas, expected_lengths, linestyle="--", label=r"$L_0 / \gamma$")
    ax.set_xlabel(r"Relative velocity $\beta$")
    ax.set_ylabel("Lab-measured rod length")
    ax.set_title("Length contraction from lab-simultaneous endpoint events")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    rows = run_experiment()
    summary_path = write_summary(rows)
    figure_path = save_plot(rows)
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
