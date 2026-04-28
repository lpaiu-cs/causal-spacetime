"""Trigger-density diagnostics for same-emission echo reachability."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_echo import echo_reachability_from_emission
from causal_spacetime_lab.state_change_echo_diagnostics import (
    echo_order_resolution_summary,
    echo_reachability_summary,
)
from causal_spacetime_lab.state_change_echo_selection import (
    emission_positions_for_reference_chain,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_reference_chains import (
    ReferenceChainCandidate,
    greedy_reference_chain_candidate_from_order,
    local_system_reference_chain_candidates,
    longest_reference_chain_candidate_from_order,
    random_reference_chain_candidate_from_order,
)
from causal_spacetime_lab.state_change_reference_quality import (
    evaluate_reference_chain_candidate,
    rank_reference_chain_candidates,
)
from causal_spacetime_lab.state_change_validation import trigger_graph_summary

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for echo coverage versus trigger density."""

    num_systems: int = 10
    max_events: int = 500
    trigger_probability_values: tuple[float, ...] = (0.05, 0.10, 0.20, 0.35, 0.50)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    emission_count: int = 3
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="State-change echo coverage versus trigger density."
    )
    parser.add_argument("--num-systems", type=int, default=10)
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument(
        "--trigger-probability-values",
        nargs="+",
        default=["0.05", "0.10", "0.20", "0.35", "0.50"],
    )
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--emission-count", type=int, default=3)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability_values=_parse_float_values(
            args.trigger_probability_values
        ),
        max_triggers_per_event=args.max_triggers_per_event,
        repetitions=args.repetitions,
        seed=args.seed,
        emission_count=args.emission_count,
        output_dir=args.output_dir,
    )


def _generate_candidates(
    order_matrix: np.ndarray,
    network,
    seed: int,
) -> list[ReferenceChainCandidate]:
    candidates = local_system_reference_chain_candidates(network)
    candidates.append(greedy_reference_chain_candidate_from_order(order_matrix))
    candidates.append(longest_reference_chain_candidate_from_order(order_matrix))
    for index in range(5):
        candidates.append(
            random_reference_chain_candidate_from_order(
                order_matrix,
                seed=seed + index,
                name=f"random_reference_chain_{index}",
            )
        )
    return candidates


def _selected_candidates(
    candidates: list[ReferenceChainCandidate],
    ranked_rows: list[dict[str, float | str]],
) -> list[tuple[str, ReferenceChainCandidate, dict[str, float | str]]]:
    by_name = {candidate.name: candidate for candidate in candidates}
    selected: list[tuple[str, ReferenceChainCandidate, dict[str, float | str]]] = []
    if ranked_rows:
        selected.append(
            (
                "highest_utility_overall",
                by_name[str(ranked_rows[0]["name"])],
                ranked_rows[0],
            )
        )
    local_rows = [row for row in ranked_rows if row["source"] == "local_system"]
    if local_rows and str(local_rows[0]["name"]) != str(ranked_rows[0]["name"]):
        selected.append(
            (
                "highest_utility_local_system",
                by_name[str(local_rows[0]["name"])],
                local_rows[0],
            )
        )
    return selected


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run echo coverage diagnostics over trigger probabilities."""

    rows: list[dict[str, float | str]] = []
    for trigger_probability in config.trigger_probability_values:
        for repetition in range(config.repetitions):
            seed = config.seed + int(1000 * trigger_probability) + repetition
            network = generate_state_change_network(
                config.num_systems,
                config.max_events,
                trigger_probability=trigger_probability,
                max_triggers_per_event=config.max_triggers_per_event,
                seed=seed,
            )
            adjacency = immediate_trigger_adjacency(network)
            closure = transitive_closure_dag(adjacency)
            graph_summary = trigger_graph_summary(adjacency, closure)
            candidates = _generate_candidates(closure, network, seed + 10000)
            reports = [
                evaluate_reference_chain_candidate(network, closure, candidate)
                for candidate in candidates
            ]
            ranked_rows = rank_reference_chain_candidates(reports)
            for label, candidate, rank_row in _selected_candidates(
                candidates,
                ranked_rows,
            ):
                emissions = emission_positions_for_reference_chain(
                    candidate.chain_event_ids,
                    strategy="interior_quantiles",
                    count=config.emission_count,
                )
                for emission_position in emissions:
                    returns, delays, reachable = echo_reachability_from_emission(
                        closure,
                        candidate.chain_event_ids,
                        int(emission_position),
                    )
                    reach_summary = echo_reachability_summary(
                        returns,
                        delays,
                        reachable,
                        emission_position=int(emission_position),
                        reference_chain_length=candidate.chain_event_ids.size,
                    )
                    resolution = echo_order_resolution_summary(delays, reachable)
                    row: dict[str, float | str] = {
                        "num_systems": float(config.num_systems),
                        "max_events": float(config.max_events),
                        "trigger_probability": float(trigger_probability),
                        "max_triggers_per_event": float(
                            config.max_triggers_per_event
                        ),
                        "repetition": float(repetition),
                        "selection_label": label,
                        "candidate_name": candidate.name,
                        "candidate_source": candidate.source,
                        "utility_score": float(rank_row["score"]),
                        "relation_density": graph_summary["relation_density"],
                        "chain_length": float(candidate.chain_event_ids.size),
                    }
                    row.update(reach_summary)
                    row.update(
                        {
                            "tied_pair_fraction": resolution["tied_pair_fraction"],
                            "strict_order_pair_fraction": resolution[
                                "strict_order_pair_fraction"
                            ],
                        }
                    )
                    rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write echo coverage rows."""

    output_path = (
        output_dir / "data" / "state_change_echo_coverage_vs_trigger_density.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else ["empty"]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_probability(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["trigger_probability"]) for row in rows})
    means = []
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
) -> tuple[Path, ...]:
    """Save echo trigger-density figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "state_change_echo_reachability_vs_trigger_density.png",
            "reachable_fraction",
            "Reachable fraction",
        ),
        (
            "state_change_echo_resolution_vs_trigger_density.png",
            "strict_order_pair_fraction",
            "Strictly ordered pair fraction",
        ),
    )
    paths: list[Path] = []
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        xs, ys = _mean_by_probability(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel("External trigger probability")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)
    return tuple(paths)


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote state-change echo coverage versus trigger density: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
