"""Scale and monotone-invariance checks for ordinal representations."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal import order_inversion_rate

OUTPUT_DATA = Path("outputs/data/metric_representation_scale_invariance.csv")
OUTPUT_FIGURE = Path("outputs/figures/metric_representation_order_preservation.png")


def _reference_ratio(values: np.ndarray) -> float:
    numerator = values[3] - values[1]
    denominator = values[4] - values[0]
    return float(numerator / denominator)


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic order-preservation transformations."""

    values = np.asarray([0.4, 0.9, 1.7, 2.8, 4.2], dtype=float)
    original_ratio = _reference_ratio(values)
    transforms: list[tuple[str, np.ndarray]] = []
    for alpha in (0.5, 1.0, 2.0, 10.0):
        transforms.append((f"scale_{alpha:g}", alpha * values))
    transforms.extend(
        [
            ("monotone_square", values * values),
            ("nonmonotone_sine", np.sin(values)),
        ]
    )

    rows: list[dict[str, float | str]] = []
    for name, transformed in transforms:
        inversion = order_inversion_rate(values, transformed, ignore_ties=True)
        transformed_ratio = _reference_ratio(transformed)
        ratio_error = abs(transformed_ratio - original_ratio)
        rows.append(
            {
                "transformation": name,
                "order_inversion_rate": float(inversion),
                "order_preserved": float(inversion == 0.0),
                "original_ratio": original_ratio,
                "transformed_ratio": transformed_ratio,
                "ratio_error": float(ratio_error),
                "ratio_preserved": float(ratio_error < 1e-12),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write scale-invariance rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_plot(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_FIGURE,
) -> Path:
    """Save order-preservation and ratio-error plot."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    labels = [str(row["transformation"]) for row in rows]
    inversion = np.asarray([row["order_inversion_rate"] for row in rows], dtype=float)
    ratio_error = np.asarray([row["ratio_error"] for row in rows], dtype=float)
    x = np.arange(len(labels))
    fig, ax1 = plt.subplots(figsize=(8.0, 4.8))
    ax1.bar(x - 0.18, inversion, width=0.36, label="order inversion")
    ax1.set_ylabel("Order inversion/error rate")
    ax1.set_ylim(0.0, 1.05)
    ax2 = ax1.twinx()
    ax2.bar(x + 0.18, ratio_error, width=0.36, color="tab:orange", label="ratio error")
    ax2.set_ylabel("Reference ratio error")
    ax1.set_xticks(x, labels, rotation=30, ha="right")
    ax1.set_title("Ordinal preservation is weaker than ratio preservation")
    ax1.grid(True, axis="y", alpha=0.3)
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    figure_path = save_plot(rows)
    print(f"Wrote metric representation scale-invariance data: {output_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
