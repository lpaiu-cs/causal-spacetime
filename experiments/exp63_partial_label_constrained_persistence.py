"""Partial identity labels restrict persistence-matching ambiguity."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.identity_constraints import (
    best_persistence_matchings_with_fixed_points,
    generate_partial_identity_constraints,
)
from causal_spacetime_lab.persistence_history import (
    adjacent_true_permutation_from_tracks,
    build_identity_tracks_from_adjacent_matchings,
)
from causal_spacetime_lab.persistence_matching import (
    matching_accuracy,
    matching_ambiguity_gap,
)
from causal_spacetime_lab.persistence_scenarios import (
    generate_unlabeled_small_motion_history,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for partial-label persistence diagnostics."""

    slice_count: int = 8
    object_count: int = 6
    known_fractions: tuple[float, ...] = (0.0, 0.2, 0.5, 0.8)
    motion_scale: float = 0.10
    repetitions: int = 10
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Partial-label constrained persistence matching."
    )
    parser.add_argument("--slice-count", type=int, default=8)
    parser.add_argument("--object-count", type=int, default=6)
    parser.add_argument(
        "--known-fractions",
        nargs="+",
        default=["0.0", "0.2", "0.5", "0.8"],
    )
    parser.add_argument("--motion-scale", type=float, default=0.10)
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_count=args.object_count,
        known_fractions=_parse_float_values(args.known_fractions),
        motion_scale=args.motion_scale,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _initial_ids_from_track(track: dict[int, int]) -> np.ndarray:
    return np.asarray(
        [object_id for object_id, _ in sorted(track.items(), key=lambda item: item[1])],
        dtype=int,
    )


def _nanmean_or_nan(values: list[float]) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else float("nan")


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run partial-label constrained persistence diagnostics."""

    rows: list[dict[str, float]] = []
    for known_fraction in config.known_fractions:
        for repetition in range(config.repetitions):
            rng = np.random.default_rng(config.seed + 1000 * repetition)
            initial = np.sort(rng.uniform(-1.0, 1.0, size=config.object_count))
            velocities = rng.normal(scale=config.motion_scale, size=config.object_count)
            positions, true_tracks = generate_unlabeled_small_motion_history(
                config.slice_count,
                initial,
                velocities,
                seed=config.seed + repetition,
            )
            inferred_adjacent: dict[tuple[int, int], np.ndarray] = {}
            accuracies: list[float] = []
            gaps: list[float] = []
            fixed_counts: list[int] = []
            labels = sorted(positions)
            for step, (left, right) in enumerate(
                zip(labels[:-1], labels[1:], strict=False)
            ):
                true_perm = adjacent_true_permutation_from_tracks(
                    true_tracks,
                    left,
                    right,
                )
                fixed = generate_partial_identity_constraints(
                    true_perm,
                    known_fraction,
                    seed=config.seed + 10000 * repetition + step,
                )
                matches = best_persistence_matchings_with_fixed_points(
                    positions[left],
                    positions[right],
                    fixed,
                    top_k=3,
                )
                best_perm = np.asarray(matches[0]["permutation"], dtype=int)
                inferred_adjacent[(left, right)] = best_perm
                accuracies.append(matching_accuracy(best_perm, true_perm))
                gaps.append(matching_ambiguity_gap(matches))
                fixed_counts.append(len(fixed))
            inferred_tracks = build_identity_tracks_from_adjacent_matchings(
                inferred_adjacent,
                _initial_ids_from_track(true_tracks[labels[0]]),
            )
            track_error = _track_error(inferred_tracks, true_tracks)
            rows.append(
                {
                    "known_fraction": float(known_fraction),
                    "repetition": float(repetition),
                    "motion_scale": float(config.motion_scale),
                    "mean_fixed_count": float(np.mean(fixed_counts)),
                    "adjacent_matching_accuracy": float(np.mean(accuracies)),
                    "ambiguity_gap": _nanmean_or_nan(gaps),
                    "track_consistency_error": track_error,
                }
            )
    return rows


def _track_error(
    inferred_tracks: dict[int, dict[int, int]],
    true_tracks: dict[int, dict[int, int]],
) -> float:
    total = 0
    wrong = 0
    for label in sorted(set(inferred_tracks) & set(true_tracks)):
        for object_id, local_index in inferred_tracks[label].items():
            total += 1
            wrong += int(true_tracks[label].get(object_id) != local_index)
    return float(wrong / total) if total else 0.0


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write partial-label diagnostics."""

    output_path = output_dir / "data" / "partial_label_constrained_persistence.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> tuple[Path, Path]:
    """Save partial-label persistence figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = (
        figure_dir / "partial_label_matching_accuracy.png",
        figure_dir / "partial_label_ambiguity_gap.png",
    )
    specs = (
        (paths[0], "adjacent_matching_accuracy", "Adjacent matching accuracy"),
        (paths[1], "ambiguity_gap", "Ambiguity gap"),
    )
    fractions = sorted({row["known_fraction"] for row in rows})
    for path, key, ylabel in specs:
        values = [
            float(
                _nanmean_or_nan(
                    [row[key] for row in rows if row["known_fraction"] == fraction]
                )
            )
            for fraction in fractions
        ]
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(fractions, values, marker="o")
        ax.set_xlabel("Known identity fraction")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote partial-label persistence data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
