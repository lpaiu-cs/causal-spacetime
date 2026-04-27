"""Crossing histories can destabilize persistence matching hypotheses."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.persistence_history import (
    adjacent_true_permutation_from_tracks,
    build_identity_tracks_from_adjacent_matchings,
    relational_history_from_inferred_tracks,
    track_consistency_error,
)
from causal_spacetime_lab.persistence_matching import (
    best_persistence_matchings_by_relational_order,
    matching_accuracy,
    matching_ambiguity_gap,
)
from causal_spacetime_lab.persistence_scenarios import (
    generate_unlabeled_crossing_history_1d,
)
from causal_spacetime_lab.relational_evolution import (
    compare_histories_order_disagreement,
    relational_shape_history,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for crossing persistence diagnostics."""

    slice_count: int = 10
    object_count: int = 5
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Crossing persistence failure.")
    parser.add_argument("--slice-count", type=int, default=10)
    parser.add_argument("--object-count", type=int, default=5)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_count=args.object_count,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _initial_ids_from_track(track: dict[int, int]) -> np.ndarray:
    return np.asarray(
        [object_id for object_id, _ in sorted(track.items(), key=lambda item: item[1])],
        dtype=int,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run crossing-history persistence diagnostics."""

    rows: list[dict[str, float]] = []
    base_positions = np.linspace(-1.0, 1.0, config.object_count)
    moving = (
        max(0, config.object_count // 2 - 1),
        min(config.object_count - 1, config.object_count // 2 + 1),
    )
    for repetition in range(config.repetitions):
        positions, true_tracks = generate_unlabeled_crossing_history_1d(
            config.slice_count,
            base_positions,
            moving,
            seed=config.seed + repetition,
        )
        inferred_adjacent: dict[tuple[int, int], np.ndarray] = {}
        accuracies: list[float] = []
        gaps: list[float] = []
        labels = sorted(positions)
        for left, right in zip(labels[:-1], labels[1:], strict=False):
            matches = best_persistence_matchings_by_relational_order(
                positions[left],
                positions[right],
                top_k=3,
            )
            best_perm = np.asarray(matches[0]["permutation"], dtype=int)
            true_perm = adjacent_true_permutation_from_tracks(true_tracks, left, right)
            inferred_adjacent[(left, right)] = best_perm
            accuracies.append(matching_accuracy(best_perm, true_perm))
            gaps.append(matching_ambiguity_gap(matches))
        inferred_tracks = build_identity_tracks_from_adjacent_matchings(
            inferred_adjacent,
            _initial_ids_from_track(true_tracks[labels[0]]),
        )
        true_history = relational_shape_history(
            relational_history_from_inferred_tracks(positions, true_tracks)
        )
        inferred_history = relational_shape_history(
            relational_history_from_inferred_tracks(positions, inferred_tracks)
        )
        rows.append(
            {
                "repetition": float(repetition),
                "mean_adjacent_accuracy": float(np.mean(accuracies)),
                "mean_ambiguity_gap": float(np.nanmean(gaps)),
                "track_consistency_error": track_consistency_error(
                    inferred_tracks,
                    true_tracks,
                ),
                "relational_history_error": compare_histories_order_disagreement(
                    true_history,
                    inferred_history,
                ),
                "failure_or_swap": float(np.mean(accuracies) < 1.0),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write crossing persistence diagnostics."""

    output_path = output_dir / "data" / "crossing_persistence_failure.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save crossing persistence track-error figure."""

    output_path = output_dir / "figures" / "crossing_persistence_track_error.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    reps = [row["repetition"] for row in rows]
    ax.plot(
        reps,
        [row["track_consistency_error"] for row in rows],
        marker="o",
        label="track error",
    )
    ax.plot(
        reps,
        [row["relational_history_error"] for row in rows],
        marker="s",
        label="history error",
    )
    ax.set_xlabel("Repetition")
    ax.set_ylabel("Error")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_path = save_figure(rows, config.output_dir)
    print(f"Wrote crossing persistence failure data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
