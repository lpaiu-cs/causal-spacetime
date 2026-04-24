"""Exploratory spacelike-distance proxy comparison in a 1+1D causal set."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.metrics import minkowski_spacelike_distance_1p1
from causal_spacetime_lab.spacelike import (
    alexandrov_interval_count_matrix,
    common_future_overlap_count,
    common_past_overlap_count,
    is_spacelike_pair,
    minimal_enclosing_alexandrov_interval_count,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

T = 4.0
N_EVENTS = 700
N_PAIRS = 350
SEED = 20260424
OUTPUT_DATA = Path("outputs/data/spacelike_distance_proxy_summary.csv")
OUTPUT_FIGURE = Path("outputs/figures/spacelike_distance_proxy_scatter.png")


def sample_spacelike_pairs(
    causal_matrix: np.ndarray,
    rng: np.random.Generator,
    n_pairs: int,
) -> list[tuple[int, int]]:
    """Sample spacelike index pairs from a causal matrix."""

    upper_i, upper_j = np.triu_indices(causal_matrix.shape[0], k=1)
    spacelike_mask = np.asarray(
        [
            is_spacelike_pair(causal_matrix, int(i), int(j))
            for i, j in zip(upper_i, upper_j, strict=True)
        ],
        dtype=bool,
    )
    spacelike_pairs = np.column_stack(
        (upper_i[spacelike_mask], upper_j[spacelike_mask])
    )
    if spacelike_pairs.shape[0] == 0:
        raise RuntimeError("no spacelike pairs were found")

    sample_size = min(n_pairs, spacelike_pairs.shape[0])
    chosen = rng.choice(spacelike_pairs.shape[0], size=sample_size, replace=False)
    return [(int(i), int(j)) for i, j in spacelike_pairs[chosen]]


def run_experiment() -> list[dict[str, float]]:
    """Compute exploratory spacelike proxy counts for sampled event pairs."""

    rng = np.random.default_rng(SEED)
    boundary_events = np.array([[-T / 2.0, 0.0], [T / 2.0, 0.0]])
    interior_events = sprinkle_1p1_causal_diamond(N_EVENTS, T=T, seed=rng)
    events = np.vstack((boundary_events, interior_events))
    causal_matrix = causal_matrix_1p1(events)
    interval_counts = alexandrov_interval_count_matrix(causal_matrix)

    candidate_pairs = sample_spacelike_pairs(causal_matrix[2:, 2:], rng, N_PAIRS)
    sampled_pairs = [(i + 2, j + 2) for i, j in candidate_pairs]

    rows: list[dict[str, float]] = []
    for i, j in sampled_pairs:
        true_distance = minkowski_spacelike_distance_1p1(events[i], events[j])
        enclosing_count = minimal_enclosing_alexandrov_interval_count(
            causal_matrix,
            i,
            j,
            interval_counts,
        )
        if enclosing_count is None:
            continue

        rows.append(
            {
                "i": float(i),
                "j": float(j),
                "t_i": float(events[i, 0]),
                "x_i": float(events[i, 1]),
                "t_j": float(events[j, 0]),
                "x_j": float(events[j, 1]),
                "true_spacelike_distance": true_distance,
                "common_future_count": float(
                    common_future_overlap_count(causal_matrix, i, j)
                ),
                "common_past_count": float(
                    common_past_overlap_count(causal_matrix, i, j)
                ),
                "minimal_enclosing_interval_count": float(enclosing_count),
            }
        )

    return rows


def write_summary(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_DATA,
) -> Path:
    """Write proxy comparison rows to CSV."""

    if not rows:
        raise RuntimeError("no rows to write")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _pearson_label(x: np.ndarray, y: np.ndarray) -> str:
    if x.size < 2 or np.all(y == y[0]):
        return "r unavailable"
    return f"r = {np.corrcoef(x, y)[0, 1]:.3f}"


def save_proxy_plot(
    rows: list[dict[str, float]],
    output_path: Path = OUTPUT_FIGURE,
) -> Path:
    """Save scatter plots comparing proxy counts to true spacelike distance."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    distances = np.asarray([row["true_spacelike_distance"] for row in rows])
    proxy_specs = [
        ("common_future_count", "Common future overlap"),
        ("common_past_count", "Common past overlap"),
        ("minimal_enclosing_interval_count", "Minimal enclosing interval"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharex=True)
    for ax, (column, title) in zip(axes, proxy_specs, strict=True):
        values = np.asarray([row[column] for row in rows])
        ax.scatter(distances, values, s=14, alpha=0.72)
        ax.set_title(f"{title}\n{_pearson_label(distances, values)}")
        ax.set_xlabel("True spacelike distance")
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel("Proxy count")
    fig.suptitle("Exploratory spacelike-distance proxy counts")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    rows = run_experiment()
    summary_path = write_summary(rows)
    figure_path = save_proxy_plot(rows)
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote figure: {figure_path}")
    print(
        "These proxies are exploratory and boundary-dependent; "
        "they are not validated distance estimators."
    )


if __name__ == "__main__":
    main()
