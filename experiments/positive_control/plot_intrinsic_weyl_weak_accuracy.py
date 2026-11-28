"""Render only an audited completed campaign; emit its table and resource totals."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render(folder, audit_path, output):
    summary_path = folder / "summary.json"
    summary = json.loads(summary_path.read_text())
    audit = json.loads(audit_path.read_text())
    assert audit["audit"] == "VERIFIED"
    assert (
        hashlib.sha256(summary_path.read_bytes()).hexdigest() == audit["summary_sha256"]
    )
    assert summary["precision_gate"] == audit["precision_gate"]
    rows = sorted(summary["results"], key=lambda row: row["case"])
    assert [row["case"] for row in rows] == list(range(8))
    assert all("failure" not in row and "precision" in row for row in rows)
    plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.8), width_ratios=[3, 1.15])
    for axis, group, weak in ((axes[0], rows[2:], True), (axes[1], rows[:2], False)):
        x = np.arange(len(group))
        scale = 200 if weak else 1
        error = np.array([row["precision"]["error"] for row in group]) * scale
        radius = np.array([row["precision"]["radius"] for row in group]) * scale
        target = 25 if weak else 0.125
        axis.bar(x, error, color="#174d70", width=0.65, label="Point error")
        axis.bar(
            x,
            radius,
            bottom=error,
            color="#8bc5e3",
            width=0.65,
            label="Adjusted confidence radius",
        )
        axis.axhline(
            target,
            color="#ab3f32",
            linewidth=1.5,
            linestyle="--",
            label="Frozen accuracy limit",
        )
        for xx, total in zip(x, error + radius, strict=True):
            axis.text(
                xx,
                total + target * 0.035,
                f"{total:.2f}%" if weak else f"{total:.4f}",
                ha="center",
                fontsize=10,
            )
        labels = []
        for row in group:
            source = "Fresh" if row["retained_pilot_case"] is None else "Retained"
            labels.append(
                f"{source}\nw={row['rotation_rate']:+g}\n"
                f"N={row['measurement_n'] / 1e6:.2f}M"
            )
        axis.set_xticks(x, labels)
        axis.set_ylim(0, max(target, float(np.max(error + radius))) * 1.25)
        axis.set_ylabel(
            "Error / truth norm (%)" if weak else "Absolute full-tidal error"
        )
        axis.set_title(
            "Six weak-curvature cases (A=0.25)" if weak else "Two flat controls (A=0)",
            pad=14,
        )
        axis.grid(axis="y", alpha=0.18)
        axis.set_axisbelow(True)
        axis.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.12),
        frameon=False,
        ncol=3,
    )
    fig.suptitle(
        f"Weak-curvature precision: {audit['precision_gate']}", fontsize=17, y=0.97
    )
    fig.text(
        0.5,
        0.055,
        "Nominal 95% simultaneous family confidence; 99.375% per case. "
        "2% assumed relative bias plus fixed numerical slack.\n"
        "Height = point error + confidence radius. "
        "Gaussian/covariance assumptions and tested local plane-wave scope apply.",
        ha="center",
        fontsize=10,
        color="#454545",
    )
    fig.subplots_adjust(left=0.065, right=0.985, top=0.83, bottom=0.32, wspace=0.3)
    output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".svg", ".png"):
        fig.savefig(output.with_suffix(suffix), dpi=180, facecolor="white")
    plt.close(fig)
    metrics = {
        "precision_gate": summary["precision_gate"],
        "measurement_events": summary["sampled_measurement_events"],
        "selection_events": summary["new_selection_events"],
        "coordinator_elapsed_seconds": summary["elapsed_seconds"],
        "generation_wall_seconds": sum(
            row["generation"]["wall_seconds"] for row in rows
        ),
        "fit_and_evaluation_wall_seconds": sum(row["wall_seconds"] for row in rows),
        "process_cpu_seconds": sum(
            row["generation"]["process_cpu_seconds"] + row["process_cpu_seconds"]
            for row in rows
        ),
        "directional_predicates": sum(
            row["generation"]["directional_predicates_classified"] for row in rows
        ),
        "peak_worker_GiB": max(row["peak_working_set_bytes"] for row in rows) / 2**30,
        "peak_generator_GiB": max(
            row["generation"]["parent_lifetime_peak_working_set_bytes"] for row in rows
        )
        / 2**30,
        "archive_GiB": sum(row["archive_bytes"] for row in rows) / 2**30,
        "audit_seconds": audit["wall_seconds"],
        "pure_replay_seconds": audit["post_completion_pure_replay"]["wall_seconds"],
        "rows": audit["rows"],
    }
    output.with_suffix(".metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.folder, args.audit, args.output)
