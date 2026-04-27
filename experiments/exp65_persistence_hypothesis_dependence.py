"""Relational histories depend on the selected persistence hypothesis."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.persistence_history import (
    build_identity_tracks_from_adjacent_matchings,
    relational_history_from_inferred_tracks,
)
from causal_spacetime_lab.persistence_matching import (
    best_persistence_matchings_by_relational_order,
)
from causal_spacetime_lab.persistence_scenarios import (
    generate_unlabeled_crossing_history_1d,
)
from causal_spacetime_lab.relational_evolution import (
    compare_histories_order_disagreement,
    relational_change_rate_between_slices,
    relational_shape_history,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for persistence-hypothesis dependence."""

    slice_count: int = 6
    object_count: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Persistence-hypothesis dependence of relational history."
    )
    parser.add_argument("--slice-count", type=int, default=6)
    parser.add_argument("--object-count", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_count=args.object_count,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _initial_ids_from_track(track: dict[int, int]) -> np.ndarray:
    return np.asarray(
        [object_id for object_id, _ in sorted(track.items(), key=lambda item: item[1])],
        dtype=int,
    )


def _mean_change_rate(history: dict[int, dict[int, float]]) -> float:
    signatures = relational_shape_history(history)
    changes = relational_change_rate_between_slices(signatures)
    if not changes:
        return 0.0
    return float(np.mean([row["order_disagreement"] for row in changes]))


def _alternative_tracks(
    true_tracks: dict[int, dict[int, int]],
) -> dict[int, dict[int, int]]:
    labels = sorted(true_tracks)
    object_ids = sorted(next(iter(true_tracks.values())))
    if len(object_ids) < 2:
        return {label: dict(track) for label, track in true_tracks.items()}
    first_id, second_id = object_ids[0], object_ids[-1]
    switch_at = labels[len(labels) // 2]
    alternative: dict[int, dict[int, int]] = {}
    for label, track in true_tracks.items():
        copied = dict(track)
        if label >= switch_at:
            copied[first_id], copied[second_id] = copied[second_id], copied[first_id]
        alternative[int(label)] = copied
    return alternative


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run persistence-hypothesis dependence diagnostics."""

    base_positions = np.linspace(-1.0, 1.0, config.object_count)
    moving = (
        max(0, config.object_count // 2 - 1),
        min(config.object_count - 1, config.object_count // 2 + 1),
    )
    positions, true_tracks = generate_unlabeled_crossing_history_1d(
        config.slice_count,
        base_positions,
        moving,
        seed=config.seed,
    )
    labels = sorted(positions)
    inferred_adjacent = {}
    for left, right in zip(labels[:-1], labels[1:], strict=False):
        matches = best_persistence_matchings_by_relational_order(
            positions[left],
            positions[right],
            top_k=1,
        )
        inferred_adjacent[(left, right)] = np.asarray(
            matches[0]["permutation"],
            dtype=int,
        )
    inferred_tracks = build_identity_tracks_from_adjacent_matchings(
        inferred_adjacent,
        _initial_ids_from_track(true_tracks[labels[0]]),
    )
    alternative_tracks = _alternative_tracks(true_tracks)
    hypotheses = {
        "hidden_true_validation": true_tracks,
        "inferred_relational_continuity": inferred_tracks,
        "deliberately_permuted_alternative": alternative_tracks,
    }
    true_signatures = relational_shape_history(
        relational_history_from_inferred_tracks(positions, true_tracks)
    )
    rows: list[dict[str, float | str]] = []
    for name, tracks in hypotheses.items():
        history = relational_history_from_inferred_tracks(positions, tracks)
        signatures = relational_shape_history(history)
        rows.append(
            {
                "hypothesis": name,
                "mean_relational_change_rate": _mean_change_rate(history),
                "disagreement_vs_hidden_true": compare_histories_order_disagreement(
                    true_signatures,
                    signatures,
                ),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write persistence-hypothesis dependence rows."""

    output_path = output_dir / "data" / "persistence_hypothesis_dependence.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Save persistence-hypothesis change-rate figure."""

    output_path = output_dir / "figures" / "persistence_hypothesis_change_rates.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.7))
    labels = [str(row["hypothesis"]) for row in rows]
    values = [float(row["mean_relational_change_rate"]) for row in rows]
    ax.bar(labels, values)
    ax.set_ylabel("Mean relational change rate")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_path = save_figure(rows, config.output_dir)
    print(f"Wrote persistence-hypothesis dependence data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
