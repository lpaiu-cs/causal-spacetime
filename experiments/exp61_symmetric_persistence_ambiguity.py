"""Symmetric configurations create persistence-matching ambiguity."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.persistence_history import (
    adjacent_true_permutation_from_tracks,
)
from causal_spacetime_lab.persistence_matching import (
    all_permutations,
    matching_accuracy,
    matching_ambiguity_gap,
    signature_disagreement_for_permutation,
)
from causal_spacetime_lab.persistence_scenarios import (
    generate_symmetric_configuration_history,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for symmetric persistence ambiguity."""

    slice_count: int = 4
    object_counts: tuple[int, ...] = (4, 5, 6)
    spacing: float = 1.0
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Symmetric persistence ambiguity.")
    parser.add_argument("--slice-count", type=int, default=4)
    parser.add_argument("--object-counts", nargs="+", default=["4", "5", "6"])
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        slice_count=args.slice_count,
        object_counts=_parse_int_values(args.object_counts),
        spacing=args.spacing,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _all_costs(previous: np.ndarray, current: np.ndarray) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for perm in all_permutations(previous.size):
        rows.append(
            {
                "permutation": perm,
                "cost": signature_disagreement_for_permutation(
                    previous,
                    current,
                    perm,
                ),
            }
        )
    rows.sort(key=lambda row: (float(row["cost"]), tuple(row["permutation"])))
    return rows


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run symmetric ambiguity diagnostics."""

    rows: list[dict[str, float]] = []
    for object_count in config.object_counts:
        for repetition in range(config.repetitions):
            positions, true_tracks = generate_symmetric_configuration_history(
                config.slice_count,
                object_count,
                spacing=config.spacing,
                seed=config.seed + 1000 * object_count + repetition,
            )
            for left, right in zip(
                range(config.slice_count - 1),
                range(1, config.slice_count),
                strict=True,
            ):
                costs = _all_costs(positions[left], positions[right])
                best_cost = float(costs[0]["cost"])
                equal_best = sum(
                    abs(float(row["cost"]) - best_cost) <= 1e-12 for row in costs
                )
                true_perm = adjacent_true_permutation_from_tracks(
                    true_tracks,
                    left,
                    right,
                )
                rows.append(
                    {
                        "object_count": float(object_count),
                        "repetition": float(repetition),
                        "slice_a": float(left),
                        "slice_b": float(right),
                        "best_cost": best_cost,
                        "second_best_cost": float(costs[1]["cost"]),
                        "ambiguity_gap": matching_ambiguity_gap(costs[:2]),
                        "equal_best_count": float(equal_best),
                        "best_accuracy_vs_hidden": matching_accuracy(
                            np.asarray(costs[0]["permutation"]),
                            true_perm,
                        ),
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write symmetric ambiguity rows."""

    output_path = output_dir / "data" / "symmetric_persistence_ambiguity.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save ambiguity gap figure."""

    output_path = output_dir / "figures" / "symmetric_persistence_ambiguity_gap.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    counts = sorted({row["object_count"] for row in rows})
    gaps = [
        float(
            np.mean(
                [row["ambiguity_gap"] for row in rows if row["object_count"] == count]
            )
        )
        for count in counts
    ]
    ax.plot(counts, gaps, marker="o")
    ax.set_xlabel("Object count")
    ax.set_ylabel("Ambiguity gap")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_path = save_figure(rows, config.output_dir)
    print(f"Wrote symmetric persistence ambiguity data: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
