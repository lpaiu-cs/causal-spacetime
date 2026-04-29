"""Threshold sensitivity for predeclared handoff eligibility."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pairwise_response_experiment_helpers import (
    default_pairwise_protocols,
    profile_from_synthetic_config,
)
from response_handoff_experiment_helpers import failure_reason_rows

from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
)
from causal_spacetime_lab.state_change_response_handoff import (
    decide_handoff_eligibility,
)
from causal_spacetime_lab.state_change_response_handoff_pipeline import (
    build_candidate_handoff_manifest,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for handoff threshold sensitivity."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32, 40)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    min_agreement_values: tuple[float, ...] = (0.6, 0.7, 0.8, 0.9)
    min_null_z_values: tuple[float, ...] = (0.0, 1.0, 2.0)
    min_bootstrap_values: tuple[float, ...] = (0.6, 0.7, 0.8)
    repetitions: int = 5
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

    parser = argparse.ArgumentParser(description="Handoff threshold sensitivity.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument(
        "--emission-positions",
        nargs="+",
        default=["8", "16", "24", "32", "40"],
    )
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument(
        "--min-agreement-values",
        nargs="+",
        default=["0.6", "0.7", "0.8", "0.9"],
    )
    parser.add_argument("--min-null-z-values", nargs="+", default=["0.0", "1.0", "2.0"])
    parser.add_argument(
        "--min-bootstrap-values",
        nargs="+",
        default=["0.6", "0.7", "0.8"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        min_agreement_values=_parse_float_values(args.min_agreement_values),
        min_null_z_values=_parse_float_values(args.min_null_z_values),
        min_bootstrap_values=_parse_float_values(args.min_bootstrap_values),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run handoff threshold sensitivity."""

    rows: list[dict[str, float | str]] = []
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
            manifest = build_candidate_handoff_manifest(
                profile,
                protocol,
                ConstraintValidationGate("base", min_constraint_count=1),
                max_constraints=1000,
                min_margin=0.05,
                constraint_seed=config.seed + repetition,
                bootstrap_count=10,
                null_repetitions=3,
                source_label=f"threshold_rep{repetition}",
            )
            summary = manifest.validation_summary
            for min_agreement in config.min_agreement_values:
                for min_null_z in config.min_null_z_values:
                    for min_bootstrap in config.min_bootstrap_values:
                        gate = ConstraintValidationGate(
                            "sweep",
                            min_constraint_count=100,
                            min_agreement_fraction=min_agreement,
                            min_null_z_score=min_null_z,
                            min_bootstrap_confidence=min_bootstrap,
                        )
                        decision = decide_handoff_eligibility(summary, gate)
                        base = {
                            "repetition": float(repetition),
                            "protocol_name": protocol.name,
                            "min_agreement": float(min_agreement),
                            "min_null_z": float(min_null_z),
                            "min_bootstrap": float(min_bootstrap),
                            "eligible": float(decision.eligible),
                            "constraint_count": float(summary.constraint_count),
                            "heldout_agreement": float(
                                summary.heldout_agreement_fraction
                            ),
                            "bootstrap_agreement": float(
                                summary.bootstrap_mean_agreement_fraction
                            ),
                            "null_z_score": float(summary.null_z_score),
                        }
                        rows.extend(failure_reason_rows(base, decision.failed_reasons))
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write threshold sensitivity CSV."""

    path = output_dir / "data" / "response_handoff_threshold_sensitivity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save threshold sensitivity figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    groups = sorted({float(row["min_agreement"]) for row in rows})
    labels = [f"{value:g}" for value in groups]
    pass_rates = [
        float(
            np.mean(
                [
                    float(row["eligible"])
                    for row in rows
                    if float(row["min_agreement"]) == value
                ]
            )
        )
        for value in groups
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels, pass_rates)
    ax.set_xlabel("minimum held-out agreement")
    ax.set_ylabel("handoff pass rate")
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_handoff_pass_rate_by_threshold.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    counts = Counter(str(row["failure_reason"]) for row in rows)
    reason_labels = list(counts)
    reason_values = [float(counts[label]) for label in reason_labels]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(reason_labels, reason_values)
    ax.set_ylabel("row count")
    ax.tick_params(axis="x", labelrotation=40)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_handoff_failure_reasons.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote response handoff threshold sensitivity: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
