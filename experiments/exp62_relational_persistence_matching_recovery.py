"""Relational-order persistence matching recovery diagnostics."""

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
    generate_unlabeled_small_motion_history,
)
from causal_spacetime_lab.relational_evolution import (
    compare_histories_order_disagreement,
    relational_shape_history,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for persistence matching recovery."""

    slice_count: int = 8
    object_counts: tuple[int, ...] = (4, 5, 6)
    motion_scales: tuple[float, ...] = (0.01, 0.05, 0.10, 0.20)
    repetitions: int = 10
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Relational persistence matching.")
    parser.add_argument("--slice-count", type=int, default=8)
    parser.add_argument("--object-counts", nargs="+", default=["4", "5", "6"])
    parser.add_argument(
        "--motion-scales",
        nargs="+",
        default=["0.01", "0.05", "0.10", "0.20"],
    )
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_counts=_parse_int_values(args.object_counts),
        motion_scales=_parse_float_values(args.motion_scales),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _initial_ids_from_track(track: dict[int, int]) -> np.ndarray:
    return np.asarray(
        [object_id for object_id, _ in sorted(track.items(), key=lambda item: item[1])],
        dtype=int,
    )


def _run_one(
    positions: dict[int, np.ndarray],
    true_tracks: dict[int, dict[int, int]],
) -> tuple[float, float, float, float]:
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
    track_error = track_consistency_error(inferred_tracks, true_tracks)
    true_history = relational_shape_history(
        relational_history_from_inferred_tracks(positions, true_tracks)
    )
    inferred_history = relational_shape_history(
        relational_history_from_inferred_tracks(positions, inferred_tracks)
    )
    history_error = compare_histories_order_disagreement(true_history, inferred_history)
    return (
        float(np.mean(accuracies)),
        track_error,
        float(np.nanmean(gaps)),
        history_error,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run persistence matching recovery diagnostics."""

    rows: list[dict[str, float]] = []
    for object_count in config.object_counts:
        for motion_scale in config.motion_scales:
            for repetition in range(config.repetitions):
                rng = np.random.default_rng(
                    config.seed + 10000 * object_count + 100 * repetition
                )
                initial = np.sort(rng.uniform(-1.0, 1.0, size=object_count))
                velocities = rng.normal(scale=motion_scale, size=object_count)
                positions, true_tracks = generate_unlabeled_small_motion_history(
                    config.slice_count,
                    initial,
                    velocities,
                    seed=config.seed + repetition,
                )
                accuracy, track_error, gap, history_error = _run_one(
                    positions,
                    true_tracks,
                )
                rows.append(
                    {
                        "object_count": float(object_count),
                        "motion_scale": float(motion_scale),
                        "repetition": float(repetition),
                        "adjacent_matching_accuracy": accuracy,
                        "track_consistency_error": track_error,
                        "ambiguity_gap": gap,
                        "relational_history_error": history_error,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write persistence matching recovery rows."""

    output_path = output_dir / "data" / "relational_persistence_matching_recovery.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save persistence matching recovery figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = (
        figure_dir / "persistence_matching_accuracy_vs_motion.png",
        figure_dir / "persistence_matching_ambiguity_gap_vs_motion.png",
    )
    for path, key, ylabel in [
        (paths[0], "adjacent_matching_accuracy", "Adjacent matching accuracy"),
        (paths[1], "ambiguity_gap", "Ambiguity gap"),
    ]:
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for object_count in sorted({row["object_count"] for row in rows}):
            subset = [row for row in rows if row["object_count"] == object_count]
            scales = sorted({row["motion_scale"] for row in subset})
            values = [
                float(
                    np.nanmean(
                        [row[key] for row in subset if row["motion_scale"] == scale]
                    )
                )
                for scale in scales
            ]
            ax.plot(scales, values, marker="o", label=f"n={int(object_count)}")
        ax.set_xlabel("Motion scale")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote relational persistence matching data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
