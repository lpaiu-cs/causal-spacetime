"""Delay-rank spacing and reference-subsampling tie-resolution diagnostics."""

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
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_response_layers import (
    EchoResponseLayerSpec,
    build_layered_echo_motif_specs,
    planted_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_response_signature import (
    echo_response_signature_from_motifs,
    signature_strict_pair_fraction,
    signature_tie_fraction,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
)
from causal_spacetime_lab.state_change_response_signature_stability import (
    signature_after_reference_subsampling,
)

DEFAULT_OUTPUT_DIR = Path("outputs")
LAYER_SETS = {
    "compact": (3, 4, 5, 6),
    "medium": (3, 5, 8, 13),
    "wide": (10, 20, 30, 40),
}


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for delay-spacing tie-resolution diagnostics."""

    reference_length_values: tuple[int, ...] = (64, 128)
    emission_position: int = 8
    layer_sets: tuple[str, ...] = ("compact", "medium", "wide")
    targets_per_layer: int = 20
    reference_strides: tuple[int, ...] = (1, 2, 4)
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

    parser = argparse.ArgumentParser(description="Echo delay spacing diagnostics.")
    parser.add_argument("--reference-length-values", nargs="+", default=["64", "128"])
    parser.add_argument("--emission-position", type=int, default=8)
    parser.add_argument(
        "--layer-sets",
        nargs="+",
        choices=sorted(LAYER_SETS),
        default=["compact", "medium", "wide"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=20)
    parser.add_argument("--reference-strides", nargs="+", default=["1", "2", "4"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length_values=_parse_int_values(args.reference_length_values),
        emission_position=args.emission_position,
        layer_sets=tuple(args.layer_sets),
        targets_per_layer=args.targets_per_layer,
        reference_strides=_parse_int_values(args.reference_strides),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _layer_specs(layer_set: str, targets_per_layer: int) -> list[EchoResponseLayerSpec]:
    return [
        EchoResponseLayerSpec(delay, targets_per_layer)
        for delay in LAYER_SETS[layer_set]
    ]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run delay spacing and tie-resolution diagnostics."""

    rows: list[dict[str, float | str]] = []
    for reference_length in config.reference_length_values:
        for layer_set in config.layer_sets:
            max_delay = max(LAYER_SETS[layer_set])
            if config.emission_position + max_delay >= reference_length:
                continue
            for stride in config.reference_strides:
                for repetition in range(config.repetitions):
                    network, reference = generate_reference_backbone_network(
                        reference_length
                    )
                    specs = build_layered_echo_motif_specs(
                        reference,
                        config.emission_position,
                        _layer_specs(layer_set, config.targets_per_layer),
                        seed=config.seed + repetition,
                    )
                    network, motifs = insert_multiple_echo_motifs(
                        network,
                        reference,
                        specs,
                    )
                    closure = transitive_closure_dag(
                        immediate_trigger_adjacency(network)
                    )
                    planted = planted_response_signature_from_motifs(motifs)
                    fine = echo_response_signature_from_motifs(
                        closure,
                        reference,
                        motifs,
                        label="fine",
                    )
                    coarse = signature_after_reference_subsampling(
                        network,
                        reference,
                        motifs,
                        stride=stride,
                        label=f"stride_{stride}",
                    )
                    planted_comparison = compare_response_signatures(planted, coarse)
                    fine_comparison = compare_response_signatures(fine, coarse)
                    rows.append(
                        {
                            "reference_length": float(reference_length),
                            "layer_set": layer_set,
                            "stride": float(stride),
                            "repetition": float(repetition),
                            "tie_fraction": signature_tie_fraction(coarse),
                            "strict_pair_fraction": signature_strict_pair_fraction(
                                coarse
                            ),
                            "pair_agreement_with_planted": planted_comparison[
                                "pair_agreement_fraction"
                            ],
                            "order_inversion_vs_fine": fine_comparison[
                                "pair_disagreement_fraction"
                            ],
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "echo_delay_spacing_tie_resolution.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _mean_rows(
    rows: list[dict[str, float | str]],
    layer_set: str,
    key: str,
) -> tuple[list[float], list[float]]:
    strides = sorted({float(row["stride"]) for row in rows})
    means: list[float] = []
    for stride in strides:
        values = [
            float(row[key])
            for row in rows
            if row["layer_set"] == layer_set
            and float(row["stride"]) == stride
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return strides, means


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    layer_sets = sorted({str(row["layer_set"]) for row in rows})
    for key, filename, ylabel in [
        ("tie_fraction", "echo_delay_spacing_tie_fraction.png", "tie fraction"),
        (
            "pair_agreement_with_planted",
            "echo_delay_spacing_pair_agreement.png",
            "pair agreement with planted",
        ),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for layer_set in layer_sets:
            xs, ys = _mean_rows(rows, layer_set, key)
            ax.plot(xs, ys, marker="o", label=layer_set)
        ax.set_xlabel("reference stride")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend(title="layer set")
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
    print(f"Wrote echo delay spacing tie resolution: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
