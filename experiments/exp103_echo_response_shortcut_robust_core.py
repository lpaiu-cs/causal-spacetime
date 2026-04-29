"""Shortcut robustness diagnostics for echo-response order signatures."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo_interference import (
    return_spectrum_report_for_motif,
    summarize_return_spectrum_reports,
)
from causal_spacetime_lab.state_change_echo_motifs import insert_multiple_echo_motifs
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    inject_shortcut_returns,
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
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    stable_response_order_core,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for shortcut robust-core diagnostics."""

    reference_length: int = 96
    emission_position: int = 8
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13, 21)
    targets_per_layer: int = 20
    shortcut_probabilities: tuple[float, ...] = (0.0, 0.1, 0.3, 0.6)
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
        description="Echo response shortcut robust core."
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
        default=["0.0", "0.1", "0.3", "0.6"],
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
        shortcut_probabilities=_parse_float_values(args.shortcut_probabilities),
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _layer_specs(config: ExperimentConfig) -> list[EchoResponseLayerSpec]:
    return [
        EchoResponseLayerSpec(delay, config.targets_per_layer)
        for delay in config.layer_delay_ranks
    ]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run shortcut robustness diagnostics."""

    rows: list[dict[str, float]] = []
    for probability in config.shortcut_probabilities:
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
            clean_closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            planted = planted_response_signature_from_motifs(motifs)
            clean = echo_response_signature_from_motifs(
                clean_closure,
                reference,
                motifs,
                label="clean",
            )
            perturbed_network, _ = inject_shortcut_returns(
                network,
                reference,
                motifs,
                ShortcutInjectionSpec(probability=probability),
                seed=seed + 1000,
            )
            perturbed_closure = transitive_closure_dag(
                immediate_trigger_adjacency(perturbed_network)
            )
            shortcut = echo_response_signature_from_motifs(
                perturbed_closure,
                reference,
                motifs,
                label="shortcut",
            )
            reports = [
                return_spectrum_report_for_motif(perturbed_closure, reference, motif)
                for motif in motifs
            ]
            summary = summarize_return_spectrum_reports(reports)
            clean_vs_planted = compare_response_signatures(planted, clean)
            shortcut_vs_planted = compare_response_signatures(planted, shortcut)
            shortcut_vs_clean = compare_response_signatures(clean, shortcut)
            core = stable_response_order_core([clean, shortcut])
            rows.append(
                {
                    "shortcut_probability": float(probability),
                    "repetition": float(repetition),
                    "clean_pair_agreement_with_planted": clean_vs_planted[
                        "pair_agreement_fraction"
                    ],
                    "shortcut_pair_agreement_with_planted": shortcut_vs_planted[
                        "pair_agreement_fraction"
                    ],
                    "shortcut_pair_agreement_with_clean": shortcut_vs_clean[
                        "pair_agreement_fraction"
                    ],
                    "shortcut_pair_disagreement_with_clean": shortcut_vs_clean[
                        "pair_disagreement_fraction"
                    ],
                    "tie_changed_fraction": shortcut_vs_clean[
                        "pair_tie_changed_fraction"
                    ],
                    "stable_pair_fraction": float(core["stable_pair_fraction"]),
                    "shortcut_fraction": summary["shortcut_fraction"],
                    "mean_shortcut_depth": summary["mean_shortcut_depth"],
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "echo_response_shortcut_robust_core.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _group_means(
    rows: list[dict[str, float]],
    value: str,
) -> tuple[np.ndarray, np.ndarray]:
    groups = sorted({row["shortcut_probability"] for row in rows})
    means: list[float] = []
    for group in groups:
        values = np.asarray(
            [
                row[value]
                for row in rows
                if row["shortcut_probability"] == group
            ],
            dtype=float,
        )
        values = values[np.isfinite(values)]
        means.append(float(np.mean(values)) if values.size else float("nan"))
    return np.asarray(groups), np.asarray(means)


def save_figures(rows: list[dict[str, float]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for value, filename, ylabel in [
        (
            "shortcut_pair_agreement_with_planted",
            "echo_response_shortcut_pair_agreement.png",
            "pair agreement with planted",
        ),
        (
            "stable_pair_fraction",
            "echo_response_shortcut_stable_core.png",
            "stable pair fraction",
        ),
    ]:
        x_values, y_values = _group_means(rows, value)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(x_values, y_values, marker="o")
        ax.set_xlabel("shortcut probability")
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
    print(f"Wrote echo response shortcut robust core: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
