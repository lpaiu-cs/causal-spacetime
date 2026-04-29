"""Response-profile stability under protocol variants."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from exp115_multi_reference_response_profile_diagnostics import (
    ExperimentConfig as ProfileConfig,
)
from exp115_multi_reference_response_profile_diagnostics import (
    build_synthetic_protocol_signatures,
)

from causal_spacetime_lab.state_change_response_profiles import (
    compare_response_profiles,
    response_profile_from_signatures,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for response-profile stability diagnostics."""

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

    parser = argparse.ArgumentParser(description="Response-profile stability.")
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


def _profile_config(config: ExperimentConfig) -> ProfileConfig:
    return ProfileConfig(
        reference_length=config.reference_length,
        emission_positions=config.emission_positions,
        layer_delay_ranks=config.layer_delay_ranks,
        targets_per_layer=config.targets_per_layer,
        repetitions=1,
        seed=config.seed,
        output_dir=config.output_dir,
    )


def _signature_from_variant(
    signature: EchoResponseSignature,
    delays: np.ndarray,
) -> EchoResponseSignature:
    reachable = delays >= 0
    return EchoResponseSignature(
        target_event_ids=signature.target_event_ids,
        delay_ranks=delays.astype(int),
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delays.astype(int), reachable),
        label=signature.label,
    )


def _shortcut_variant(
    signatures: list[EchoResponseSignature],
    probability: float,
    seed: int,
) -> list[EchoResponseSignature]:
    rng = np.random.default_rng(seed)
    variants: list[EchoResponseSignature] = []
    for signature in signatures:
        delays = signature.delay_ranks.copy()
        mask = signature.reachable_mask & (rng.random(delays.size) < probability)
        reductions = rng.integers(1, 4, size=delays.size)
        delays[mask] = np.maximum(1, delays[mask] - reductions[mask])
        variants.append(_signature_from_variant(signature, delays))
    return variants


def _stride_variant(
    signatures: list[EchoResponseSignature],
    stride: int,
) -> list[EchoResponseSignature]:
    variants: list[EchoResponseSignature] = []
    for signature in signatures:
        delays = signature.delay_ranks.copy()
        mask = signature.reachable_mask
        delays[mask] = np.ceil(delays[mask] / stride).astype(int)
        variants.append(_signature_from_variant(signature, delays))
    return variants


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run response-profile stability diagnostics."""

    rows: list[dict[str, float | str]] = []
    profile_config = _profile_config(config)
    for repetition in range(config.repetitions):
        baseline_signatures = build_synthetic_protocol_signatures(
            profile_config,
            repetition,
        )
        baseline = response_profile_from_signatures(baseline_signatures)
        for probability in config.shortcut_probability_values:
            variant = response_profile_from_signatures(
                _shortcut_variant(
                    baseline_signatures,
                    probability,
                    config.seed + repetition + int(probability * 1000),
                )
            )
            comparison = compare_response_profiles(baseline, variant)
            rows.append(
                {
                    "variant_kind": "shortcut",
                    "variant_value": float(probability),
                    "repetition": float(repetition),
                    **comparison,
                }
            )
        for stride in config.reference_strides:
            variant = response_profile_from_signatures(
                _stride_variant(baseline_signatures, stride)
            )
            comparison = compare_response_profiles(baseline, variant)
            rows.append(
                {
                    "variant_kind": "reference_stride",
                    "variant_value": float(stride),
                    "repetition": float(repetition),
                    **comparison,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = (
        output_dir
        / "data"
        / "response_profile_stability_under_protocol_variants.csv"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_rows(
    rows: list[dict[str, float | str]],
    variant_kind: str,
    key: str,
) -> tuple[list[float], list[float]]:
    xs = sorted(
        {
            float(row["variant_value"])
            for row in rows
            if row["variant_kind"] == variant_kind
        }
    )
    means = [
        float(
            np.nanmean(
                [
                    float(row[key])
                    for row in rows
                    if row["variant_kind"] == variant_kind
                    and float(row["variant_value"]) == x
                ]
            )
        )
        for x in xs
    ]
    return xs, means


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "profile_entry_agreement_fraction",
            "response_profile_entry_agreement.png",
            "profile entry agreement",
        ),
        (
            "profile_pair_separation_agreement",
            "response_profile_pair_separation_agreement.png",
            "pair separation agreement",
        ),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for variant_kind in sorted({str(row["variant_kind"]) for row in rows}):
            xs, ys = _mean_rows(rows, variant_kind, key)
            ax.plot(xs, ys, marker="o", label=variant_kind)
        ax.set_xlabel("variant value")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend()
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
    print(f"Wrote response profile stability: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
