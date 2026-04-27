"""Ratio stability from equal-step calibration constraints."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal import order_inversion_rate

OUTPUT_DATA = Path("outputs/data/ratio_stability_from_calibration.csv")
OUTPUT_FIGURE = Path("outputs/figures/ratio_stability_from_calibration.png")


def _interval_ratio(values: np.ndarray) -> float:
    return float((values[5] - values[1]) / (values[3] - values[1]))


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic calibration checks."""

    ticks = np.arange(0.0, 8.0, dtype=float)
    original_ratio = _interval_ratio(ticks)
    transforms = {
        "scale_2t": 2.0 * ticks,
        "translate_t_plus_3": ticks + 3.0,
        "square_t2": ticks * ticks,
        "log_1_plus_t": np.log1p(ticks),
    }
    rows: list[dict[str, float | str]] = []
    for name, transformed in transforms.items():
        steps = np.diff(transformed)
        step_std = float(np.std(steps))
        step_mean = float(np.mean(steps))
        ratio = _interval_ratio(transformed)
        ratio_error = abs(ratio - original_ratio)
        inversion = order_inversion_rate(ticks, transformed, ignore_ties=True)
        rows.append(
            {
                "transformation": name,
                "order_inversion_rate": float(inversion),
                "order_preserved": float(inversion == 0.0),
                "mean_step": step_mean,
                "step_std": step_std,
                "equal_step_preserved": float(step_std < 1e-12),
                "original_ratio": original_ratio,
                "transformed_ratio": ratio,
                "ratio_error": float(ratio_error),
                "ratio_preserved": float(ratio_error < 1e-12),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write ratio calibration rows."""

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
    """Save calibration diagnostic plot."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    labels = [str(row["transformation"]) for row in rows]
    step_std = np.asarray([row["step_std"] for row in rows], dtype=float)
    ratio_error = np.asarray([row["ratio_error"] for row in rows], dtype=float)
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.bar(x - 0.18, step_std, width=0.36, label="step std")
    ax.bar(x + 0.18, ratio_error, width=0.36, label="ratio error")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylabel("Deviation")
    ax.set_title("Calibration restrictions stabilize ratios")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    rows = run_experiment()
    output_path = write_outputs(rows)
    figure_path = save_plot(rows)
    print(f"Wrote ratio stability data: {output_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
