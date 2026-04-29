"""Diagnostics for multi-reference response profiles."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_response_profiles import (
    EchoResponseProfile,
    response_profile_equivalence_classes,
    response_profile_from_signatures,
    response_profile_reachable_fraction,
    response_profile_separation_fraction,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    response_order_sign_matrix,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for response-profile diagnostics."""

    reference_length: int = 96
    emission_positions: tuple[int, ...] = (8, 16, 24)
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
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

    parser = argparse.ArgumentParser(description="Multi-reference profile diagnostics.")
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument("--emission-positions", nargs="+", default=["8", "16", "24"])
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_positions=_parse_int_values(args.emission_positions),
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _signature_from_delays(
    target_event_ids: np.ndarray,
    delays: np.ndarray,
    label: str,
) -> EchoResponseSignature:
    reachable = delays >= 0
    return EchoResponseSignature(
        target_event_ids=target_event_ids.astype(int),
        delay_ranks=delays.astype(int),
        reachable_mask=reachable,
        order_sign_matrix=response_order_sign_matrix(delays.astype(int), reachable),
        label=label,
    )


def build_synthetic_protocol_signatures(
    config: ExperimentConfig,
    repetition: int,
) -> list[EchoResponseSignature]:
    """Build shared-target synthetic protocol signatures.

    These protocol columns are response-profile diagnostics only. They are not
    metric coordinates or calibrated reference protocols.
    """

    rng = np.random.default_rng(config.seed + repetition)
    target_count = len(config.layer_delay_ranks) * config.targets_per_layer
    target_ids = np.arange(target_count, dtype=int)
    base_delays = np.repeat(config.layer_delay_ranks, config.targets_per_layer)
    within_layer = np.tile(
        np.arange(config.targets_per_layer),
        len(config.layer_delay_ranks),
    )
    signatures: list[EchoResponseSignature] = []
    for column, emission in enumerate(config.emission_positions):
        delays = base_delays.astype(int) * 10
        if column > 0:
            offsets = (within_layer + column) % (column + 2)
            delays = delays + offsets
        missing_count = int(max(0, min(target_count // 10, column)))
        if missing_count:
            missing = rng.choice(target_count, size=missing_count, replace=False)
            delays = delays.copy()
            delays[missing] = -1
        signatures.append(
            _signature_from_delays(target_ids, delays, f"emission_{emission}")
        )
    return signatures


def _profile_summary(
    profile: EchoResponseProfile,
    profile_kind: str,
    repetition: int,
) -> dict[str, float | str]:
    classes = response_profile_equivalence_classes(profile)
    largest = max((group.size for group in classes), default=0)
    return {
        "profile_kind": profile_kind,
        "repetition": float(repetition),
        "protocol_count": float(len(profile.protocol_labels)),
        "target_count": float(profile.target_event_ids.size),
        "reachable_fraction": response_profile_reachable_fraction(profile),
        "pair_separation_fraction": response_profile_separation_fraction(profile),
        "equivalence_class_count": float(len(classes)),
        "largest_equivalence_class_size": float(largest),
    }


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run multi-reference profile diagnostics."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        signatures = build_synthetic_protocol_signatures(config, repetition)
        single_profile = response_profile_from_signatures(signatures[:1])
        multi_profile = response_profile_from_signatures(signatures)
        rows.append(_profile_summary(single_profile, "single_protocol", repetition))
        rows.append(_profile_summary(multi_profile, "multi_protocol", repetition))
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "multi_reference_response_profile_diagnostics.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_by_kind(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    kinds = sorted({str(row["profile_kind"]) for row in rows})
    values = [
        float(np.mean([float(row[key]) for row in rows if row["profile_kind"] == kind]))
        for kind in kinds
    ]
    return kinds, values


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "pair_separation_fraction",
            "multi_reference_profile_separation.png",
            "pair separation fraction",
        ),
        (
            "equivalence_class_count",
            "multi_reference_equivalence_classes.png",
            "equivalence class count",
        ),
    ]:
        labels, values = _mean_by_kind(rows, key)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(labels, values)
        ax.set_ylabel(ylabel)
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
    print(f"Wrote multi-reference response profile diagnostics: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
