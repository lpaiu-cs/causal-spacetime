"""Null-baseline separation for response-comparison constraint pools."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import (
    default_pairwise_protocols,
    profile_from_synthetic_config,
)

from causal_spacetime_lab.state_change_response_constraint_null_validation import (
    evaluate_constraint_pool_against_nulls,
    null_separation_by_type,
)
from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    pairwise_response_dissimilarity,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for constraint null-separation validation."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    max_constraints: int = 5000
    min_margin: float = 0.05
    null_repetitions: int = 10
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

    parser = argparse.ArgumentParser(description="Constraint null separation.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument(
        "--emission-positions",
        nargs="+",
        default=["8", "16", "24", "32"],
    )
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument("--max-constraints", type=int, default=5000)
    parser.add_argument("--min-margin", type=float, default=0.05)
    parser.add_argument("--null-repetitions", type=int, default=10)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        max_constraints=args.max_constraints,
        min_margin=args.min_margin,
        null_repetitions=args.null_repetitions,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run null-separation diagnostics."""

    summary_rows: list[dict[str, float | str]] = []
    by_type_rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        profile = profile_from_synthetic_config(
            config.reference_length,
            config.emission_positions,
            config.layer_delay_ranks,
            config.targets_per_layer,
            repetition,
            config.seed,
        )
        for protocol in default_pairwise_protocols():
            pool = build_constraint_pool_from_dissimilarity(
                pairwise_response_dissimilarity(profile, protocol),
                max_constraints=config.max_constraints,
                min_margin=config.min_margin,
                seed=config.seed + repetition,
            )
            base = {
                "repetition": float(repetition),
                "protocol_name": protocol.name,
                "method": protocol.method,
                "constraint_count": float(pool.constraints.shape[0]),
                "min_margin": float(config.min_margin),
            }
            summary_rows.append(
                {
                    **base,
                    **evaluate_constraint_pool_against_nulls(
                        profile,
                        protocol,
                        pool,
                        null_repetitions=config.null_repetitions,
                        seed=config.seed + 1000 * repetition,
                    ),
                }
            )
            for row in null_separation_by_type(
                profile,
                protocol,
                pool,
                null_repetitions=config.null_repetitions,
                seed=config.seed + 1000 * repetition,
            ):
                by_type_rows.append({**base, **row})
    return summary_rows, by_type_rows


def write_outputs(
    summary_rows: list[dict[str, float | str]],
    by_type_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write null-separation CSV files."""

    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    summary_path = data_dir / "response_constraint_null_separation.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    by_type_path = data_dir / "response_constraint_null_separation_by_type.csv"
    with by_type_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(by_type_rows[0]))
        writer.writeheader()
        writer.writerows(by_type_rows)
    return summary_path, by_type_path


def _means(
    rows: list[dict[str, float | str]],
    group_key: str,
    value_key: str,
) -> tuple[list[str], list[float]]:
    labels = sorted({str(row[group_key]) for row in rows})
    values = []
    for label in labels:
        raw = np.asarray(
            [float(row[value_key]) for row in rows if row[group_key] == label],
            dtype=float,
        )
        raw[~np.isfinite(raw)] = np.nan
        values.append(float(np.nanmean(raw)))
    return labels, values


def save_figures(
    summary_rows: list[dict[str, float | str]],
    by_type_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save null-separation figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    labels, values = _means(summary_rows, "protocol_name", "null_z_score")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, values)
    ax.set_ylabel("null z score")
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_constraint_null_z_score.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    labels, values = _means(
        by_type_rows,
        "null_type",
        "mean_null_agreement_fraction",
    )
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, values)
    ax.set_ylabel("mean null agreement")
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_constraint_null_agreement_by_type.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def main() -> None:
    config = parse_args()
    summary_rows, by_type_rows = run_experiment(config)
    data_paths = write_outputs(summary_rows, by_type_rows, config.output_dir)
    figure_paths = save_figures(summary_rows, by_type_rows, config.output_dir)
    print(f"Wrote response constraint null separation: {data_paths[0]}")
    print(f"Wrote response constraint null separation by type: {data_paths[1]}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
