"""Integrated validation-gate summary for response-comparison pools."""

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

from causal_spacetime_lab.state_change_response_constraint_coverage import (
    constraint_pool_summary,
)
from causal_spacetime_lab.state_change_response_constraint_null_validation import (
    evaluate_constraint_pool_against_nulls,
)
from causal_spacetime_lab.state_change_response_constraint_pool import (
    build_constraint_pool_from_dissimilarity,
)
from causal_spacetime_lab.state_change_response_constraint_validation import (
    ConstraintValidationGate,
    bootstrap_constraint_stability,
    heldout_protocol_constraint_validation,
    validation_gate_pass_fail,
)
from causal_spacetime_lab.state_change_response_pairwise import (
    pairwise_response_dissimilarity,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for integrated validation-gate summary."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24, 32, 40)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    max_constraints: int = 5000
    min_margin_values: tuple[float, ...] = (0.0, 0.05, 0.10)
    repetitions: int = 5
    bootstrap_count: int = 10
    null_repetitions: int = 5
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

    parser = argparse.ArgumentParser(description="Constraint validation gates.")
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
    parser.add_argument("--max-constraints", type=int, default=5000)
    parser.add_argument(
        "--min-margin-values",
        nargs="+",
        default=["0.0", "0.05", "0.10"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--bootstrap-count", type=int, default=10)
    parser.add_argument("--null-repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        max_constraints=args.max_constraints,
        min_margin_values=_parse_float_values(args.min_margin_values),
        repetitions=args.repetitions,
        bootstrap_count=args.bootstrap_count,
        null_repetitions=args.null_repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run integrated validation-gate diagnostics."""

    gate = ConstraintValidationGate("integrated_default")
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
            dissimilarity = pairwise_response_dissimilarity(profile, protocol)
            for min_margin in config.min_margin_values:
                seed = config.seed + 1000 * repetition + int(min_margin * 1000)
                pool = build_constraint_pool_from_dissimilarity(
                    dissimilarity,
                    max_constraints=config.max_constraints,
                    min_margin=min_margin,
                    seed=seed,
                )
                heldout = heldout_protocol_constraint_validation(
                    profile,
                    protocol,
                    train_fraction=0.6,
                    max_constraints=config.max_constraints,
                    min_margin=min_margin,
                    seed=seed,
                )
                bootstrap = bootstrap_constraint_stability(
                    profile,
                    protocol,
                    pool,
                    bootstrap_count=config.bootstrap_count,
                    seed=seed,
                )
                nulls = evaluate_constraint_pool_against_nulls(
                    profile,
                    protocol,
                    pool,
                    null_repetitions=config.null_repetitions,
                    seed=seed,
                )
                gate_input = {**heldout, **bootstrap, **nulls}
                gate_row = validation_gate_pass_fail(gate_input, gate)
                rows.append(
                    {
                        "repetition": float(repetition),
                        "protocol_name": protocol.name,
                        "method": protocol.method,
                        "margin_threshold": float(min_margin),
                        **constraint_pool_summary(
                            pool,
                            int(profile.target_event_ids.size),
                        ),
                        **{f"heldout_{key}": value for key, value in heldout.items()},
                        **bootstrap,
                        **nulls,
                        **gate_row,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write integrated validation-gate CSV."""

    path = output_dir / "data" / "response_constraint_validation_gate_summary.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save integrated validation-gate figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    labels = sorted({str(row["protocol_name"]) for row in rows})
    pass_rates = [
        float(
            np.mean(
                [float(row["passed"]) for row in rows if row["protocol_name"] == label]
            )
        )
        for label in labels
    ]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, pass_rates)
    ax.set_ylabel("gate pass rate")
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_constraint_gate_pass_rate.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    reasons: Counter[str] = Counter()
    for row in rows:
        for reason in str(row["failed_reasons"]).split(";"):
            if reason:
                reasons[reason] += 1
    reason_labels = list(reasons) or ["none"]
    reason_values = [float(reasons[label]) for label in reason_labels] or [0.0]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(reason_labels, reason_values)
    ax.set_ylabel("failure count")
    ax.tick_params(axis="x", labelrotation=40)
    ax.grid(True, axis="y", alpha=0.3)
    path = figure_dir / "response_constraint_gate_failure_reasons.png"
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
    print(f"Wrote response constraint validation gate summary: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
