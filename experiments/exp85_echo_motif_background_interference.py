"""Background plus controlled shortcut diagnostics for echo-response motifs.

This Milestone 23 experiment uses random background trigger networks and also
adds controlled shortcut-return edges as a stress test. Milestone 24 separates
that targeted injection from generic acyclic background edge perturbations.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import StateChangeNetwork, TriggerEdge
from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
    motif_recovery_report,
    summarize_motif_recovery,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifRecord,
    EchoMotifSpec,
    insert_multiple_echo_motifs,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_background_state_change_network_with_reference,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_reference_chains import (
    ReferenceChainCandidate,
)
from causal_spacetime_lab.state_change_reference_quality import (
    evaluate_reference_chain_candidate,
    rank_reference_chain_candidates,
)
from causal_spacetime_lab.state_change_validation import trigger_graph_summary

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for echo motif background-interference diagnostics."""

    num_systems: int = 8
    max_events: int = 300
    trigger_probability_values: tuple[float, ...] = (0.05, 0.15, 0.30)
    motif_count: int = 30
    delay_ranks: tuple[int, ...] = (2, 3, 5, 8)
    repetitions: int = 5
    seed: int = 0
    reference_source: str = "highest_utility"
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

    parser = argparse.ArgumentParser(description="Echo motif background interference.")
    parser.add_argument("--num-systems", type=int, default=8)
    parser.add_argument("--max-events", type=int, default=300)
    parser.add_argument(
        "--trigger-probability-values",
        nargs="+",
        default=["0.05", "0.15", "0.30"],
    )
    parser.add_argument("--motif-count", type=int, default=30)
    parser.add_argument("--delay-ranks", nargs="+", default=["2", "3", "5", "8"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--reference-source", default="highest_utility")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability_values=_parse_float_values(
            args.trigger_probability_values
        ),
        motif_count=args.motif_count,
        delay_ranks=_parse_int_values(args.delay_ranks),
        repetitions=args.repetitions,
        seed=args.seed,
        reference_source=args.reference_source,
        output_dir=args.output_dir,
    )


def _random_specs(
    reference_length: int,
    motif_count: int,
    delay_ranks: tuple[int, ...],
    seed: int,
) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in delay_ranks if delay < reference_length]
    )
    if valid_delays.size == 0:
        return []
    specs: list[EchoMotifSpec] = []
    for _ in range(motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, reference_length - delay))
        specs.append(EchoMotifSpec(emission, delay))
    return specs


def _add_shortcut_return_edges(
    network: StateChangeNetwork,
    reference: np.ndarray,
    motifs: list[EchoMotifRecord],
    probability: float,
    seed: int,
) -> tuple[StateChangeNetwork, int]:
    """Add controlled shortcut-return trigger paths for interference tests."""

    rng = np.random.default_rng(seed)
    edges = list(network.trigger_edges)
    added = 0
    for motif in motifs:
        if motif.planted_delay_rank <= 1 or rng.random() > probability:
            continue
        shortcut_position = int(
            rng.integers(
                motif.emission_position + 1,
                motif.planted_return_position,
            )
        )
        trial_edges = edges + [
            TriggerEdge(
                motif.target_event_id,
                int(reference[shortcut_position]),
                "external_trigger",
            )
        ]
        trial = StateChangeNetwork(
            network.events,
            trial_edges,
            network.system_event_ids,
        )
        try:
            topological_order_from_adjacency(immediate_trigger_adjacency(trial))
        except ValueError:
            continue
        edges = trial_edges
        added += 1
    return StateChangeNetwork(network.events, edges, network.system_event_ids), added


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run background-interference motif recovery diagnostics."""

    rows: list[dict[str, float | str]] = []
    for trigger_probability in config.trigger_probability_values:
        for repetition in range(config.repetitions):
            seed = config.seed + int(1000 * trigger_probability) + repetition
            network, reference = (
                generate_background_state_change_network_with_reference(
                    config.num_systems,
                    config.max_events,
                    trigger_probability,
                    reference_source=config.reference_source,
                    seed=seed,
                )
            )
            if reference.size < 3:
                continue
            base_closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            reference_report = evaluate_reference_chain_candidate(
                network,
                base_closure,
                ReferenceChainCandidate(
                    "selected_reference",
                    reference,
                    "manual",
                ),
            )
            specs = _random_specs(
                reference.size,
                config.motif_count,
                config.delay_ranks,
                seed + 2000,
            )
            network, motifs = insert_multiple_echo_motifs(network, reference, specs)
            shortcut_probability = min(0.85, trigger_probability * 1.5)
            network, shortcut_edges_added = _add_shortcut_return_edges(
                network,
                reference,
                motifs,
                shortcut_probability,
                seed + 3000,
            )
            adjacency = immediate_trigger_adjacency(network)
            closure = transitive_closure_dag(adjacency)
            graph_summary = trigger_graph_summary(adjacency, closure)
            motif_rows = [
                motif_recovery_report(closure, reference, motif)
                for motif in motifs
            ]
            summary = summarize_motif_recovery(motif_rows)
            order_summary = motif_order_recovery_rate(motif_rows)
            rows.append(
                {
                    "num_systems": float(config.num_systems),
                    "max_events": float(config.max_events),
                    "trigger_probability": float(trigger_probability),
                    "repetition": float(repetition),
                    "reference_source": config.reference_source,
                    "reference_chain_length": float(reference.size),
                    "reference_chain_utility_score": (
                        rank_reference_chain_candidates([reference_report])[0]["score"]
                    ),
                    "shortcut_edges_added": float(shortcut_edges_added),
                    "shortcut_edge_probability": float(shortcut_probability),
                    "relation_density": graph_summary["relation_density"],
                    **summary,
                    **order_summary,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write background-interference rows."""

    output_path = output_dir / "data" / "echo_motif_background_interference.csv"
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
    probabilities = sorted({float(row["trigger_probability"]) for row in rows})
    means: list[float] = []
    for probability in probabilities:
        values = [
            float(row[key])
            for row in rows
            if float(row["trigger_probability"]) == probability
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return probabilities, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save background-interference figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    shortcut_path = figure_dir / "echo_motif_shortcut_fraction_vs_trigger_density.png"
    recovery_path = figure_dir / "echo_motif_recovery_vs_trigger_density.png"
    for path, key, ylabel in (
        (shortcut_path, "early_shortcut_fraction", "Early shortcut fraction"),
        (recovery_path, "exact_recovery_fraction", "Exact recovery fraction"),
    ):
        xs, ys = _mean_by_probability(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel("External trigger probability")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return shortcut_path, recovery_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo motif background interference: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
