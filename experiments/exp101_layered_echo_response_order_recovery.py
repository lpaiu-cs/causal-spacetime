"""Clean layered echo-response order recovery diagnostics."""

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
    signature_reachable_fraction,
    signature_strict_pair_fraction,
    signature_tie_fraction,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    response_order_cycle_count,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for layered response-order recovery."""

    reference_length: int = 96
    emission_position: int = 8
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer_values: tuple[int, ...] = (5, 20, 50)
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

    parser = argparse.ArgumentParser(
        description="Layered echo-response order recovery."
    )
    parser.add_argument("--reference-length", type=int, default=96)
    parser.add_argument("--emission-position", type=int, default=8)
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument(
        "--targets-per-layer-values",
        nargs="+",
        default=["5", "20", "50"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        emission_position=args.emission_position,
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer_values=_parse_int_values(args.targets_per_layer_values),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _layer_specs(
    delay_ranks: tuple[int, ...],
    targets_per_layer: int,
) -> list[EchoResponseLayerSpec]:
    return [
        EchoResponseLayerSpec(delay, targets_per_layer, layer_label=f"layer_{delay}")
        for delay in delay_ranks
    ]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run clean layered response-order recovery."""

    rows: list[dict[str, float]] = []
    for targets_per_layer in config.targets_per_layer_values:
        for repetition in range(config.repetitions):
            network, reference = generate_reference_backbone_network(
                config.reference_length
            )
            specs = build_layered_echo_motif_specs(
                reference,
                config.emission_position,
                _layer_specs(config.layer_delay_ranks, targets_per_layer),
                seed=config.seed + repetition,
            )
            network, motifs = insert_multiple_echo_motifs(network, reference, specs)
            closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            planted = planted_response_signature_from_motifs(motifs)
            recovered = echo_response_signature_from_motifs(
                closure,
                reference,
                motifs,
                label="recovered",
            )
            comparison = compare_response_signatures(planted, recovered)
            rows.append(
                {
                    "targets_per_layer": float(targets_per_layer),
                    "repetition": float(repetition),
                    "target_count": float(len(motifs)),
                    "reachable_fraction": signature_reachable_fraction(recovered),
                    "pair_agreement_fraction": comparison["pair_agreement_fraction"],
                    "pair_disagreement_fraction": comparison[
                        "pair_disagreement_fraction"
                    ],
                    "tie_fraction": signature_tie_fraction(recovered),
                    "strict_pair_fraction": signature_strict_pair_fraction(recovered),
                    "cycle_count": float(
                        response_order_cycle_count(recovered.order_sign_matrix)
                    ),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "layered_echo_response_order_recovery.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _group_means(
    rows: list[dict[str, float]],
    key: str,
    value: str,
) -> tuple[np.ndarray, np.ndarray]:
    groups = sorted({row[key] for row in rows})
    means = [
        float(np.nanmean([row[value] for row in rows if row[key] == group]))
        for group in groups
    ]
    return np.asarray(groups, dtype=float), np.asarray(means, dtype=float)


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for value, filename, ylabel in [
        (
            "pair_agreement_fraction",
            "layered_echo_pair_agreement.png",
            "pair agreement",
        ),
        ("tie_fraction", "layered_echo_tie_fraction.png", "tie fraction"),
    ]:
        x_values, y_values = _group_means(rows, "targets_per_layer", value)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(x_values, y_values, marker="o")
        ax.set_xlabel("targets per planted layer")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
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
    print(f"Wrote layered echo response order recovery: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()

