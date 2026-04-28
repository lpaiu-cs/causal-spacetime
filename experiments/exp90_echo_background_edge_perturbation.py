"""Generic acyclic background edge perturbations for echo motifs."""

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
from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_multiple_echo_motifs,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    add_random_acyclic_background_edges,
    classify_added_edge_effects_on_motifs,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for generic background edge perturbation diagnostics."""

    reference_length: int = 64
    motif_count: int = 60
    delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    edge_probabilities: tuple[float, ...] = (0.0, 0.001, 0.003, 0.01, 0.03)
    repetitions: int = 5
    seed: int = 0
    max_edges: int = 500
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo background edge perturbation.")
    parser.add_argument("--reference-length", type=int, default=64)
    parser.add_argument("--motif-count", type=int, default=60)
    parser.add_argument("--delay-ranks", nargs="+", default=["3", "5", "8", "13"])
    parser.add_argument(
        "--edge-probabilities",
        nargs="+",
        default=["0.0", "0.001", "0.003", "0.01", "0.03"],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-edges", type=int, default=500)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        reference_length=args.reference_length,
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        edge_probabilities=_parse_float_values(args.edge_probabilities),
        repetitions=args.repetitions,
        seed=args.seed,
        max_edges=args.max_edges,
        output_dir=args.output_dir,
    )


def _random_specs(config: ExperimentConfig, seed: int) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in config.delay_ranks if delay < config.reference_length],
        dtype=int,
    )
    specs: list[EchoMotifSpec] = []
    for _ in range(config.motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, config.reference_length - delay))
        specs.append(EchoMotifSpec(emission, delay))
    return specs


def _effect_fraction(rows: list[dict[str, float | str]], effect: str) -> float:
    if not rows:
        return float("nan")
    return float(np.mean([row["effect"] == effect for row in rows]))


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run generic background edge perturbation diagnostics."""

    rows: list[dict[str, float | str]] = []
    for probability in config.edge_probabilities:
        for repetition in range(config.repetitions):
            seed = config.seed + 1000 * repetition + int(probability * 100000)
            network, reference = generate_reference_backbone_network(
                config.reference_length
            )
            network, motifs = insert_multiple_echo_motifs(
                network,
                reference,
                _random_specs(config, seed),
            )
            before = transitive_closure_dag(immediate_trigger_adjacency(network))
            before_reports = [
                return_spectrum_report_for_motif(before, reference, motif)
                for motif in motifs
            ]
            before_summary = summarize_return_spectrum_reports(before_reports)
            before_order = motif_order_recovery_rate(before_reports)
            perturbed, added_edge_count = add_random_acyclic_background_edges(
                network,
                probability,
                seed=seed + 500,
                max_edges=config.max_edges,
            )
            after = transitive_closure_dag(immediate_trigger_adjacency(perturbed))
            after_reports = [
                return_spectrum_report_for_motif(after, reference, motif)
                for motif in motifs
            ]
            after_summary = summarize_return_spectrum_reports(after_reports)
            after_order = motif_order_recovery_rate(after_reports)
            effects = classify_added_edge_effects_on_motifs(
                before,
                after,
                reference,
                motifs,
            )
            rows.append(
                {
                    "reference_length": float(config.reference_length),
                    "motif_count_setting": float(config.motif_count),
                    "edge_probability": float(probability),
                    "repetition": float(repetition),
                    "added_edge_count": float(added_edge_count),
                    "new_shortcut_fraction": _effect_fraction(
                        effects,
                        "new_shortcut",
                    ),
                    "shortcut_deepened_fraction": _effect_fraction(
                        effects,
                        "shortcut_deepened",
                    ),
                    "unchanged_fraction": _effect_fraction(effects, "unchanged"),
                    "newly_reachable_fraction": _effect_fraction(
                        effects,
                        "newly_reachable",
                    ),
                    "shortcut_removed_fraction": _effect_fraction(
                        effects,
                        "shortcut_removed",
                    ),
                    "became_missing_fraction": _effect_fraction(
                        effects,
                        "became_missing",
                    ),
                    "exact_recovery_before": before_summary[
                        "exact_recovery_fraction"
                    ],
                    "exact_recovery_after": after_summary["exact_recovery_fraction"],
                    "shortcut_fraction_after": after_summary["shortcut_fraction"],
                    "mean_spectrum_size_after": after_summary["mean_spectrum_size"],
                    "order_agreement_before": before_order["order_agreement_rate"],
                    "order_agreement_after": after_order["order_agreement_rate"],
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write background perturbation rows."""

    output_path = output_dir / "data" / "echo_background_edge_perturbation.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_probability(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["edge_probability"]) for row in rows})
    means: list[float] = []
    for probability in probabilities:
        values = [
            float(row[key])
            for row in rows
            if float(row["edge_probability"]) == probability
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return probabilities, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save background perturbation figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    shortcut_path = (
        figure_dir / "echo_background_new_shortcuts_vs_edge_probability.png"
    )
    recovery_path = figure_dir / "echo_background_recovery_before_after.png"
    xs, ys = _mean_by_probability(rows, "new_shortcut_fraction")
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.plot(xs, ys, marker="o")
    ax.set_xlabel("Generic background edge probability")
    ax.set_ylabel("New shortcut fraction")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(shortcut_path, dpi=200)
    plt.close(fig)

    xs, before = _mean_by_probability(rows, "exact_recovery_before")
    _, after = _mean_by_probability(rows, "exact_recovery_after")
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.plot(xs, before, marker="o", label="Before")
    ax.plot(xs, after, marker="s", label="After")
    ax.set_xlabel("Generic background edge probability")
    ax.set_ylabel("Exact recovery fraction")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(recovery_path, dpi=200)
    plt.close(fig)
    return shortcut_path, recovery_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo background edge perturbation: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
