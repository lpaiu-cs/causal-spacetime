"""Order-level bracket diagnostics for selected reference chains."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import (
    StateChangeNetwork,
    generate_state_change_network,
)
from causal_spacetime_lab.state_change_bracket_diagnostics import (
    bracket_coverage_summary,
    rank_slice_summary,
)
from causal_spacetime_lab.state_change_brackets import (
    assign_reference_rank_slices,
    radar_time_rank_from_reference_brackets,
    reference_tick_brackets_from_order,
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

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for reference-chain bracket diagnostics."""

    num_systems_values: tuple[int, ...] = (5, 10)
    max_events_values: tuple[int, ...] = (100, 300, 600)
    trigger_probability_values: tuple[float, ...] = (0.15, 0.30)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    random_candidate_count: int = 5
    seed: int = 0
    bin_width: int = 2
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
        description="State-change reference-chain bracket diagnostics."
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
    parser.add_argument("--bin-width", type=int, default=2)
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
        bin_width=args.bin_width,
        output_dir=args.output_dir,
    )


def generate_reference_candidates(
    order_matrix: np.ndarray,
    network: StateChangeNetwork,
    random_candidate_count: int,
    seed: int,
) -> list[ReferenceChainCandidate]:
    """Generate local-system, order-only, and random reference chains."""

    candidates = local_system_reference_chain_candidates(network)
    candidates.append(greedy_reference_chain_candidate_from_order(order_matrix))
    candidates.append(longest_reference_chain_candidate_from_order(order_matrix))
    for index in range(random_candidate_count):
        candidates.append(
            random_reference_chain_candidate_from_order(
                order_matrix,
                seed=seed + index,
                name=f"random_reference_chain_{index}",
            )
        )
    return candidates


def _select_reference_candidates(
    candidates: list[ReferenceChainCandidate],
    ranked_rows: list[dict[str, float | str]],
) -> list[tuple[str, ReferenceChainCandidate, dict[str, float | str]]]:
    by_name = {candidate.name: candidate for candidate in candidates}
    selected: list[tuple[str, ReferenceChainCandidate, dict[str, float | str]]] = []
    seen: set[str] = set()

    def add(label: str, row: dict[str, float | str]) -> None:
        name = str(row["name"])
        if name in seen:
            return
        seen.add(name)
        selected.append((label, by_name[name], row))

    add("highest_utility_overall", ranked_rows[0])
    for source in sorted({str(row["source"]) for row in ranked_rows}):
        source_rows = [row for row in ranked_rows if row["source"] == source]
        add(f"highest_utility_{source}", source_rows[0])
    random_rows = [row for row in ranked_rows if row["source"] == "random_order"]
    for index, row in enumerate(random_rows[:2]):
        add(f"random_reference_{index}", row)
    return selected


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run reference-chain bracket diagnostics."""

    rows: list[dict[str, float | str]] = []
    for num_systems in config.num_systems_values:
        for max_events in config.max_events_values:
            for trigger_probability in config.trigger_probability_values:
                for repetition in range(config.repetitions):
                    seed = (
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
                        seed=seed,
                    )
                    closure = transitive_closure_dag(
                        immediate_trigger_adjacency(network)
                    )
                    candidates = generate_reference_candidates(
                        closure,
                        network,
                        config.random_candidate_count,
                        seed + 5000,
                    )
                    reports = [
                        evaluate_reference_chain_candidate(network, closure, candidate)
                        for candidate in candidates
                    ]
                    ranked_rows = rank_reference_chain_candidates(reports)
                    for label, candidate, rank_row in _select_reference_candidates(
                        candidates,
                        ranked_rows,
                    ):
                        predecessors, successors, accessible = (
                            reference_tick_brackets_from_order(
                                closure,
                                candidate.chain_event_ids,
                            )
                        )
                        time_ranks = radar_time_rank_from_reference_brackets(
                            predecessors,
                            successors,
                            accessible,
                        )
                        slices = assign_reference_rank_slices(
                            time_ranks,
                            accessible,
                            bin_width=config.bin_width,
                        )
                        row: dict[str, float | str] = {
                            "num_systems": float(num_systems),
                            "max_events": float(max_events),
                            "trigger_probability": float(trigger_probability),
                            "max_triggers_per_event": float(
                                config.max_triggers_per_event
                            ),
                            "repetition": float(repetition),
                            "selection_label": label,
                            "candidate_name": candidate.name,
                            "candidate_source": candidate.source,
                            "utility_score": float(rank_row["score"]),
                            "utility_rank": float(rank_row["rank"]),
                            "bin_width": float(config.bin_width),
                        }
                        row.update(
                            bracket_coverage_summary(
                                predecessors,
                                successors,
                                accessible,
                                reference_chain_length=candidate.chain_event_ids.size,
                            )
                        )
                        row.update(rank_slice_summary(slices))
                        rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write bracket diagnostics."""

    output_path = output_dir / "data" / "state_change_reference_bracket_diagnostics.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_source(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    sources = sorted({str(row["candidate_source"]) for row in rows})
    means = [
        float(
            np.mean(
                [
                    float(row[key])
                    for row in rows
                    if row["candidate_source"] == source
                ]
            )
        )
        for source in sources
    ]
    return sources, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, ...]:
    """Save bracket diagnostic figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "state_change_reference_bracket_accessible_fraction_by_source.png",
            "accessible_fraction",
            "Accessible fraction",
        ),
        (
            "state_change_reference_bracket_width_by_source.png",
            "mean_bracket_width_rank",
            "Mean bracket-width rank",
        ),
        (
            "state_change_reference_rank_slice_count_by_source.png",
            "slice_count",
            "Rank slice count",
        ),
    )
    paths: list[Path] = []
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        sources, means = _mean_by_source(rows, key)
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.bar(sources, means)
        ax.set_xlabel("Reference-chain source")
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.3)
        fig.autofmt_xdate(rotation=20)
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
    print(f"Wrote state-change reference bracket diagnostics: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
