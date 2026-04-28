"""Rank reference-chain candidates in state-change trigger networks."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_observer_ambiguity import (
    chain_candidate_diversity,
    top_score_gap,
)
from causal_spacetime_lab.state_change_observer_quality import (
    evaluate_chain_candidate,
    rank_chain_candidates,
)
from causal_spacetime_lab.state_change_observers import (
    ObserverChainCandidate,
    greedy_chain_candidate_from_order,
    local_system_chain_candidates,
    longest_chain_candidate_from_order,
    random_chain_candidate_from_order,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for reference-chain candidate ranking."""

    num_systems_values: tuple[int, ...] = (5, 10)
    max_events_values: tuple[int, ...] = (100, 300, 600)
    trigger_probability_values: tuple[float, ...] = (0.15, 0.30)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    random_candidate_count: int = 5
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
        description="Observer-like chain candidate ranking."
    )
    parser.add_argument("--num-systems-values", nargs="+", default=["5", "10"])
    parser.add_argument("--max-events-values", nargs="+", default=["100", "300", "600"])
    parser.add_argument(
        "--trigger-probability-values",
        nargs="+",
        default=["0.15", "0.30"],
    )
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--random-candidate-count", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems_values=_parse_int_values(args.num_systems_values),
        max_events_values=_parse_int_values(args.max_events_values),
        trigger_probability_values=_parse_float_values(
            args.trigger_probability_values
        ),
        max_triggers_per_event=args.max_triggers_per_event,
        repetitions=args.repetitions,
        random_candidate_count=args.random_candidate_count,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def generate_candidates(
    order_matrix: np.ndarray,
    network,
    random_candidate_count: int,
    seed: int,
) -> list[ObserverChainCandidate]:
    """Generate local-system, order-only, and random baseline candidates."""

    candidates = local_system_chain_candidates(network)
    candidates.append(greedy_chain_candidate_from_order(order_matrix))
    candidates.append(longest_chain_candidate_from_order(order_matrix))
    for index in range(random_candidate_count):
        candidates.append(
            random_chain_candidate_from_order(
                order_matrix,
                seed=seed + index,
                name=f"random_order_chain_{index}",
            )
        )
    return candidates


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run candidate ranking diagnostics."""

    ranking_rows: list[dict[str, float | str]] = []
    ambiguity_rows: list[dict[str, float | str]] = []
    for num_systems in config.num_systems_values:
        for max_events in config.max_events_values:
            for trigger_probability in config.trigger_probability_values:
                for repetition in range(config.repetitions):
                    run_seed = (
                        config.seed
                        + 100000 * num_systems
                        + 100 * max_events
                        + int(1000 * trigger_probability)
                        + repetition
                    )
                    network = generate_state_change_network(
                        num_systems,
                        max_events,
                        trigger_probability=trigger_probability,
                        max_triggers_per_event=config.max_triggers_per_event,
                        seed=run_seed,
                    )
                    closure = transitive_closure_dag(
                        immediate_trigger_adjacency(network)
                    )
                    candidates = generate_candidates(
                        closure,
                        network,
                        config.random_candidate_count,
                        run_seed + 7000,
                    )
                    reports = [
                        evaluate_chain_candidate(network, closure, candidate)
                        for candidate in candidates
                    ]
                    ranked = rank_chain_candidates(reports)
                    for row in ranked:
                        output_row = dict(row)
                        output_row.update(
                            {
                                "num_systems": float(num_systems),
                                "max_events": float(max_events),
                                "trigger_probability": float(trigger_probability),
                                "max_triggers_per_event": float(
                                    config.max_triggers_per_event
                                ),
                                "repetition": float(repetition),
                            }
                        )
                        ranking_rows.append(output_row)
                    diversity = chain_candidate_diversity(candidates)
                    ambiguity_row = {
                        "num_systems": float(num_systems),
                        "max_events": float(max_events),
                        "trigger_probability": float(trigger_probability),
                        "max_triggers_per_event": float(
                            config.max_triggers_per_event
                        ),
                        "repetition": float(repetition),
                        "top_score_gap": top_score_gap(ranked),
                        "best_score": float(ranked[0]["score"]),
                        "best_source": str(ranked[0]["source"]),
                        "best_name": str(ranked[0]["name"]),
                    }
                    ambiguity_row.update(diversity)
                    ambiguity_rows.append(ambiguity_row)
    return ranking_rows, ambiguity_rows


def _write_csv(rows: list[dict[str, float | str]], output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def write_outputs(
    ranking_rows: list[dict[str, float | str]],
    ambiguity_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write candidate ranking and ambiguity rows."""

    data_dir = output_dir / "data"
    ranking_path = data_dir / "observer_chain_candidate_ranking.csv"
    ambiguity_path = data_dir / "observer_chain_candidate_ambiguity.csv"
    return _write_csv(ranking_rows, ranking_path), _write_csv(
        ambiguity_rows,
        ambiguity_path,
    )


def _mean_by_source(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    sources = sorted({str(row["source"]) for row in rows})
    means = [
        float(np.mean([float(row[key]) for row in rows if row["source"] == source]))
        for source in sources
    ]
    return sources, means


def save_figures(
    ranking_rows: list[dict[str, float | str]],
    ambiguity_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save candidate ranking figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    score_path = figure_dir / "observer_chain_score_by_source.png"
    sources, means = _mean_by_source(ranking_rows, "score")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(sources, means)
    ax.set_ylabel("Mean heuristic score")
    ax.set_xlabel("Candidate source")
    ax.grid(True, axis="y", alpha=0.3)
    fig.autofmt_xdate(rotation=20)
    fig.tight_layout()
    fig.savefig(score_path, dpi=200)
    plt.close(fig)

    bracket_path = figure_dir / "observer_chain_bracketing_by_source.png"
    sources, means = _mean_by_source(ranking_rows, "bracketed_fraction")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar(sources, means)
    ax.set_ylabel("Mean bracketed fraction")
    ax.set_xlabel("Candidate source")
    ax.grid(True, axis="y", alpha=0.3)
    fig.autofmt_xdate(rotation=20)
    fig.tight_layout()
    fig.savefig(bracket_path, dpi=200)
    plt.close(fig)

    gap_path = figure_dir / "observer_chain_ambiguity_gap.png"
    max_events_values = sorted({row["max_events"] for row in ambiguity_rows})
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for trigger_probability in sorted(
        {row["trigger_probability"] for row in ambiguity_rows}
    ):
        ys = [
            float(
                np.mean(
                    [
                        row["top_score_gap"]
                        for row in ambiguity_rows
                        if row["trigger_probability"] == trigger_probability
                        and row["max_events"] == max_events
                    ]
                )
            )
            for max_events in max_events_values
        ]
        ax.plot(max_events_values, ys, marker="o", label=f"p={trigger_probability:.2f}")
    ax.set_xlabel("Maximum stored events")
    ax.set_ylabel("Mean top-score gap")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gap_path, dpi=200)
    plt.close(fig)
    return score_path, bracket_path, gap_path


def main() -> None:
    config = parse_args()
    ranking_rows, ambiguity_rows = run_experiment(config)
    ranking_path, ambiguity_path = write_outputs(
        ranking_rows,
        ambiguity_rows,
        config.output_dir,
    )
    figure_paths = save_figures(ranking_rows, ambiguity_rows, config.output_dir)
    print(f"Wrote observer-chain candidate ranking: {ranking_path}")
    print(f"Wrote observer-chain candidate ambiguity: {ambiguity_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
