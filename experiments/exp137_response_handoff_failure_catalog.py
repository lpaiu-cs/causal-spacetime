"""Catalog failed handoff reasons under weak synthetic profiles."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from pairwise_response_experiment_helpers import synthetic_random_profile

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
    summarize_handoff_manifests,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for handoff failure catalog."""

    target_count: int = 60
    protocol_count_values: tuple[int, ...] = (2, 3, 5)
    reachable_probabilities: tuple[float, ...] = (0.3, 0.5, 0.8)
    unique_rank_counts: tuple[int, ...] = (3, 5)
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

    parser = argparse.ArgumentParser(description="Handoff failure catalog.")
    parser.add_argument("--target-count", type=int, default=60)
    parser.add_argument("--protocol-count-values", nargs="+", default=["2", "3", "5"])
    parser.add_argument(
        "--reachable-probabilities",
        nargs="+",
        default=["0.3", "0.5", "0.8"],
    )
    parser.add_argument("--unique-rank-counts", nargs="+", default=["3", "5"])
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        target_count=args.target_count,
        protocol_count_values=_parse_int_values(args.protocol_count_values),
        reachable_probabilities=_parse_float_values(args.reachable_probabilities),
        unique_rank_counts=_parse_int_values(args.unique_rank_counts),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run failed-handoff catalog diagnostics."""

    rows: list[dict[str, float | str]] = []
    protocol = PairwiseResponseComparisonProtocol("gap_mean", "rank_gap_mean")
    gate = ConstraintValidationGate("failure_catalog")
    for protocol_count in config.protocol_count_values:
        for reachable_probability in config.reachable_probabilities:
            for unique_rank_count in config.unique_rank_counts:
                for repetition in range(config.repetitions):
                    profile = synthetic_random_profile(
                        config.target_count,
                        protocol_count,
                        reachable_probability,
                        unique_rank_count,
                        config.seed
                        + repetition
                        + 10 * protocol_count
                        + int(100 * reachable_probability),
                    )
                    manifest = build_candidate_handoff_manifest(
                        profile,
                        protocol,
                        gate,
                        max_constraints=1000,
                        min_margin=0.05,
                        constraint_seed=config.seed + repetition,
                        bootstrap_count=5,
                        null_repetitions=2,
                        source_label="failure_catalog",
                    )
                    base = summarize_handoff_manifests([manifest])[0]
                    base.update(
                        {
                            "protocol_count": float(protocol_count),
                            "reachable_probability": float(reachable_probability),
                            "unique_rank_count": float(unique_rank_count),
                            "repetition": float(repetition),
                        }
                    )
                    reasons = manifest.handoff_decision.failed_reasons or ["none"]
                    for reason in reasons:
                        rows.append({**base, "failure_reason": reason})
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write failure catalog CSV."""

    path = output_dir / "data" / "response_handoff_failure_catalog.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save failure mode figure."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    counts = Counter(str(row["failure_reason"]) for row in rows)
    labels = list(counts)
    values = [float(counts[label]) for label in labels]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(labels, values)
    ax.set_ylabel("row count")
    ax.tick_params(axis="x", labelrotation=40)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_handoff_failure_mode_counts.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return [path]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote response handoff failure catalog: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
