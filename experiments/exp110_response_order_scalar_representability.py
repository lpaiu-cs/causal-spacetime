"""Scalar ordinal representability diagnostics for response-order cores."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_motifs import insert_multiple_echo_motifs
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_response_layers import (
    EchoResponseLayerSpec,
    build_layered_echo_motif_specs,
)
from causal_spacetime_lab.state_change_response_representability import (
    scalar_representability_report,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    consensus_response_order_matrix,
    stable_response_order_core,
)
from causal_spacetime_lab.state_change_response_signature_stability import (
    build_protocol_variant_signatures,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for scalar ordinal representability diagnostics."""

    reference_length: int = 96
    emission_position: int = 8
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13, 21)
    targets_per_layer: int = 20
    shortcut_probabilities: tuple[float, ...] = (0.0, 0.2, 0.5, 0.8)
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

    parser = argparse.ArgumentParser(
        description="Response-order scalar representability."
    )
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument("--emission-position", type=int, default=8)
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13", "21"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=20)
    parser.add_argument(
        "--shortcut-probabilities",
        nargs="+",
        default=["0.0", "0.2", "0.5", "0.8"],
    )
    parser.add_argument("--reference-strides", nargs="+", default=["1", "2", "4"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_position=args.emission_position,
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        shortcut_probabilities=_parse_float_values(args.shortcut_probabilities),
        reference_strides=_parse_int_values(args.reference_strides),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _layer_specs(config: ExperimentConfig) -> list[EchoResponseLayerSpec]:
    return [
        EchoResponseLayerSpec(delay, config.targets_per_layer)
        for delay in config.layer_delay_ranks
    ]


def _variant_specs(config: ExperimentConfig, seed: int) -> list[dict[str, object]]:
    specs: list[dict[str, object]] = [{"kind": "baseline", "label": "baseline"}]
    for probability in config.shortcut_probabilities:
        specs.append(
            {
                "kind": "shortcut_injection",
                "probability": probability,
                "label": f"shortcut_{probability:g}",
                "seed": seed + int(probability * 1000),
            }
        )
    for stride in config.reference_strides:
        specs.append(
            {
                "kind": "reference_subsampling",
                "stride": stride,
                "label": f"reference_stride_{stride}",
            }
        )
    return specs


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run scalar ordinal representability diagnostics."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        seed = config.seed + repetition
        network, reference = generate_reference_backbone_network(
            config.reference_length
        )
        specs = build_layered_echo_motif_specs(
            reference,
            config.emission_position,
            _layer_specs(config),
            seed=seed,
        )
        network, motifs = insert_multiple_echo_motifs(network, reference, specs)
        signatures = build_protocol_variant_signatures(
            network,
            reference,
            motifs,
            _variant_specs(config, seed + 1000),
        )
        core = stable_response_order_core(signatures)
        consensus = consensus_response_order_matrix(signatures)
        for report_type, matrix in [
            ("stable_core", np.asarray(core["stable_order_sign_matrix"], dtype=int)),
            ("consensus", consensus),
        ]:
            report = scalar_representability_report(matrix)
            rows.append(
                {
                    "repetition": float(repetition),
                    "report_type": report_type,
                    "stable_pair_fraction": float(core["stable_pair_fraction"]),
                    **report,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "response_order_scalar_representability.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_by_type(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    types = sorted({str(row["report_type"]) for row in rows})
    means: list[float] = []
    for report_type in types:
        values = [
            float(row[key])
            for row in rows
            if row["report_type"] == report_type and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return types, means


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for key, filename, ylabel in [
        (
            "directed_3cycle_count",
            "response_order_scalar_cycle_count.png",
            "directed 3-cycle count",
        ),
        (
            "scalar_representable",
            "response_order_scalar_representable_fraction.png",
            "scalar representable fraction",
        ),
    ]:
        labels, values = _mean_by_type(rows, key)
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
    print(f"Wrote response-order scalar representability: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
