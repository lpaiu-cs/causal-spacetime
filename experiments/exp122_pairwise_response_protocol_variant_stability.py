"""Pairwise response comparison stability under protocol variants."""

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

from causal_spacetime_lab.state_change_response_pairwise import (
    pairwise_response_dissimilarity,
    pairwise_response_order_inversion_rate,
    response_pair_comparison_constraints,
)
from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    pairwise_protocol_admissibility_report,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for pairwise response variant stability."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    shortcut_probability_values: tuple[float, ...] = (0.0, 0.3)
    reference_strides: tuple[int, ...] = (1, 2, 4)
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

    parser = argparse.ArgumentParser(description="Pairwise response variant stability.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument("--emission-positions", nargs="+", default=["8", "16", "24"])
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument(
        "--shortcut-probability-values",
        nargs="+",
        default=["0.0", "0.3"],
    )
    parser.add_argument("--reference-strides", nargs="+", default=["1", "2", "4"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        shortcut_probability_values=_parse_float_values(
            args.shortcut_probability_values
        ),
        reference_strides=_parse_int_values(args.reference_strides),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _shortcut_profile(
    profile: EchoResponseProfile,
    probability: float,
    seed: int,
) -> EchoResponseProfile:
    rng = np.random.default_rng(seed)
    delays = profile.delay_rank_matrix.copy()
    mask = profile.reachable_matrix & (rng.random(delays.shape) < probability)
    reductions = rng.integers(1, 4, size=delays.shape)
    delays[mask] = np.maximum(1, delays[mask] - reductions[mask])
    return EchoResponseProfile(
        profile.target_event_ids.copy(),
        list(profile.protocol_labels),
        delays,
        profile.reachable_matrix.copy(),
    )


def _stride_profile(profile: EchoResponseProfile, stride: int) -> EchoResponseProfile:
    delays = profile.delay_rank_matrix.copy()
    mask = profile.reachable_matrix
    delays[mask] = np.ceil(delays[mask] / stride).astype(int)
    return EchoResponseProfile(
        profile.target_event_ids.copy(),
        list(profile.protocol_labels),
        delays,
        profile.reachable_matrix.copy(),
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run pairwise response variant stability diagnostics."""

    rows: list[dict[str, float | str]] = []
    protocols = default_pairwise_protocols()
    for repetition in range(config.repetitions):
        baseline = profile_from_synthetic_config(
            config.reference_length,
            config.emission_positions,
            config.layer_delay_ranks,
            config.targets_per_layer,
            repetition,
            config.seed,
        )
        variants: list[tuple[str, float, EchoResponseProfile]] = []
        for probability in config.shortcut_probability_values:
            variants.append(
                (
                    "shortcut",
                    float(probability),
                    _shortcut_profile(
                        baseline,
                        probability,
                        config.seed + repetition + int(probability * 1000),
                    ),
                )
            )
        for stride in config.reference_strides:
            variants.append(
                ("reference_stride", float(stride), _stride_profile(baseline, stride))
            )
        for protocol in protocols:
            baseline_dissimilarity = pairwise_response_dissimilarity(baseline, protocol)
            for variant_kind, variant_value, variant in variants:
                dissimilarity = pairwise_response_dissimilarity(variant, protocol)
                report = pairwise_protocol_admissibility_report(variant, protocol)
                constraints = response_pair_comparison_constraints(
                    dissimilarity,
                    num_constraints=500,
                    seed=config.seed + repetition,
                )
                rows.append(
                    {
                        "variant_kind": variant_kind,
                        "variant_value": variant_value,
                        "protocol_name": protocol.name,
                        "repetition": float(repetition),
                        "order_inversion_rate": pairwise_response_order_inversion_rate(
                            baseline_dissimilarity,
                            dissimilarity,
                        ),
                        "constraint_count": float(constraints.shape[0]),
                        **report,
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "pairwise_response_protocol_variant_stability.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_by_variant(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    groups = sorted(
        {
            (str(row["variant_kind"]), str(row["variant_value"]))
            for row in rows
            if row["protocol_name"] == "gap_mean"
        }
    )
    labels = [
        f"{variant_kind}={variant_value}" for variant_kind, variant_value in groups
    ]
    values = []
    for variant_kind, variant_value in groups:
        matching_values = [
            float(row[key])
            for row in rows
            if row["protocol_name"] == "gap_mean"
            and str(row["variant_kind"]) == variant_kind
            and str(row["variant_value"]) == variant_value
        ]
        values.append(float(np.nanmean(matching_values)) if matching_values else np.nan)
    return labels, values


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save variant stability figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "order_inversion_rate",
            "pairwise_response_variant_inversion.png",
            "order inversion",
        ),
        (
            "valid_pair_fraction",
            "pairwise_response_variant_valid_pairs.png",
            "valid pair fraction",
        ),
    ]:
        labels, values = _mean_by_variant(rows, key)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(labels, values)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", labelrotation=35)
        ax.grid(True, axis="y", alpha=0.3)
        path = figure_dir / filename
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
    print(f"Wrote pairwise response protocol variant stability: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
