"""Persistence-dependent, transport-relative velocity diagnostics."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.cross_slice import velocity_without_transport
from causal_spacetime_lab.persistence import generate_persistent_object_events_1p1
from causal_spacetime_lab.slice_transport import (
    SliceTransportRule,
    velocity_under_transport,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Persistence-dependent velocity diagnostics.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return args.output_dir


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic protocol-relative velocity diagnostics."""

    slice_times = np.linspace(-0.8, 0.8, 12)
    initial_positions = np.linspace(-0.35, 0.35, 5)
    velocities = np.asarray([-0.15, -0.05, 0.0, 0.08, 0.16])
    events, slice_ids, object_ids = generate_persistent_object_events_1p1(
        slice_times,
        initial_positions,
        velocities,
    )
    labels = sorted(set(int(value) for value in slice_ids))
    identity_rule = SliceTransportRule(
        name="identity_transport",
        description="identity slice transport",
        slice_labels=slice_ids,
        scale_by_slice={label: 1.0 for label in labels},
        shift_by_slice={label: 0.0 for label in labels},
        reflection_by_slice={label: 1.0 for label in labels},
    )
    drifting_rule = SliceTransportRule(
        name="drifting_transport",
        description="slice-dependent shifted transport",
        slice_labels=slice_ids,
        scale_by_slice={label: 1.0 for label in labels},
        shift_by_slice={label: 0.05 * label for label in labels},
        reflection_by_slice={label: 1.0 for label in labels},
    )
    time_by_slice = {index: float(value) for index, value in enumerate(slice_times)}
    rows: list[dict[str, float | str]] = []
    undefined = velocity_without_transport()
    rows.append(
        {
            "transport": "none",
            "object_id": -1.0,
            "defined": float(undefined.defined),
            "velocity_mean": float("nan"),
            "velocity_std": float("nan"),
            "true_velocity": float("nan"),
        }
    )
    for rule in (identity_rule, drifting_rule):
        estimates = velocity_under_transport(
            events[:, 1],
            slice_ids,
            object_ids,
            time_by_slice,
            rule,
        )
        for object_id, summary in estimates.items():
            rows.append(
                {
                    "transport": rule.name,
                    "object_id": float(object_id),
                    "defined": 1.0,
                    "velocity_mean": summary["velocity_mean"],
                    "velocity_std": summary["velocity_std"],
                    "true_velocity": float(velocities[object_id]),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write persistence velocity rows."""

    output_path = output_dir / "data" / "persistence_dependent_velocity.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save protocol-relative velocity plot."""

    output_path = output_dir / "figures" / "persistence_velocity_by_transport.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for transport in ("identity_transport", "drifting_transport"):
        subset = [row for row in rows if row["transport"] == transport]
        ax.plot(
            [float(row["object_id"]) for row in subset],
            [float(row["velocity_mean"]) for row in subset],
            marker="o",
            label=transport,
        )
    true_rows = [row for row in rows if row["transport"] == "identity_transport"]
    ax.plot(
        [float(row["object_id"]) for row in true_rows],
        [float(row["true_velocity"]) for row in true_rows],
        linestyle="--",
        color="black",
        label="true validation velocity",
    )
    ax.set_xlabel("Object id")
    ax.set_ylabel("Velocity under chosen transport")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    output_dir = parse_args()
    rows = run_experiment()
    data_path = write_outputs(rows, output_dir)
    figure_path = save_figure(rows, output_dir)
    print(f"Wrote persistence-dependent velocity data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
