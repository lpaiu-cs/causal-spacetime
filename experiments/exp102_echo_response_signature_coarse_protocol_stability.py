"""Response-signature stability across coarse protocol variants."""

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
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    response_order_cycle_count,
    stable_response_order_core,
)
from causal_spacetime_lab.state_change_response_signature_stability import (
    build_protocol_variant_signatures,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for coarse protocol signature stability."""

    reference_length: int = 96
    emission_position: int = 8
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13, 21)
    targets_per_layer: int = 20
    event_keep_probabilities: tuple[float, ...] = (1.0, 0.5)
    reference_strides: tuple[int, ...] = (1, 2, 4)
    edge_removal_probabilities: tuple[float, ...] = (0.0, 0.10)
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
        description="Echo response signature coarse protocol stability."
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
        "--event-keep-probabilities",
        nargs="+",
        default=["1.0", "0.5"],
    )
    parser.add_argument("--reference-strides", nargs="+", default=["1", "2", "4"])
    parser.add_argument(
        "--edge-removal-probabilities",
        nargs="+",
        default=["0.0", "0.10"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_position=args.emission_position,
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        event_keep_probabilities=_parse_float_values(args.event_keep_probabilities),
        reference_strides=_parse_int_values(args.reference_strides),
        edge_removal_probabilities=_parse_float_values(
            args.edge_removal_probabilities
        ),
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
    for keep in config.event_keep_probabilities:
        specs.append(
            {
                "kind": "closure_event_thinning",
                "keep_probability": keep,
                "label": f"event_keep_{keep:g}",
                "seed": seed + int(keep * 1000),
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
    for removal in config.edge_removal_probabilities:
        for protect in (True, False):
            specs.append(
                {
                    "kind": "edge_thinning",
                    "removal_probability": removal,
                    "protect_motif_edges": protect,
                    "label": f"edge_{removal:g}_protect_{int(protect)}",
                    "seed": seed + int(removal * 1000) + int(protect),
                }
            )
    return specs


def _variant_kind(label: str) -> str:
    if label.startswith("event_keep"):
        return "closure_event_thinning"
    if label.startswith("reference_stride"):
        return "reference_subsampling"
    if label.startswith("edge"):
        return "edge_thinning"
    return "baseline"


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run response-signature stability across shared protocol variants."""

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
        baseline = signatures[0]
        core = stable_response_order_core(signatures)
        stable_fraction = float(core["stable_pair_fraction"])
        for signature in signatures:
            comparison = compare_response_signatures(baseline, signature)
            rows.append(
                {
                    "repetition": float(repetition),
                    "variant_label": signature.label,
                    "variant_kind": _variant_kind(signature.label),
                    "pair_agreement_fraction": comparison[
                        "pair_agreement_fraction"
                    ],
                    "pair_disagreement_fraction": comparison[
                        "pair_disagreement_fraction"
                    ],
                    "pair_tie_changed_fraction": comparison[
                        "pair_tie_changed_fraction"
                    ],
                    "baseline_tie_fraction": comparison["baseline_tie_fraction"],
                    "variant_tie_fraction": comparison["variant_tie_fraction"],
                    "stable_pair_fraction": stable_fraction,
                    "cycle_count": float(
                        response_order_cycle_count(signature.order_sign_matrix)
                    ),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "echo_response_signature_coarse_protocol_stability.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _group_means(
    rows: list[dict[str, float | str]],
    value: str,
) -> tuple[list[str], list[float]]:
    groups = sorted({str(row["variant_kind"]) for row in rows})
    means: list[float] = []
    for group in groups:
        values = np.asarray(
            [float(row[value]) for row in rows if row["variant_kind"] == group],
            dtype=float,
        )
        values = values[np.isfinite(values)]
        means.append(float(np.mean(values)) if values.size else float("nan"))
    return groups, means


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for value, filename, ylabel in [
        (
            "pair_agreement_fraction",
            "echo_response_signature_pair_agreement_by_variant.png",
            "pair agreement with baseline",
        ),
        (
            "stable_pair_fraction",
            "echo_response_signature_stable_core.png",
            "stable pair fraction",
        ),
    ]:
        labels, means = _group_means(rows, value)
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(labels, means)
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.tick_params(axis="x", rotation=25)
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
    print(f"Wrote echo response signature coarse protocol stability: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
