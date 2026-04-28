"""Reference-chain bracket coverage versus trigger density."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
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
    """Configuration for bracket coverage versus trigger density."""

    num_systems: int = 10
    max_events: int = 500
    trigger_probability_values: tuple[float, ...] = (0.05, 0.10, 0.20, 0.35, 0.50)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    bin_width: int = 2
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Reference-chain bracket coverage versus trigger density."
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
    parser.add_argument("--bin-width", type=int, default=2)
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
        bin_width=args.bin_width,
        output_dir=args.output_dir,
    )


def _rank_candidates(network, closure: np.ndarray, seed: int):
    candidates = local_system_reference_chain_candidates(network)
    candidates.append(greedy_reference_chain_candidate_from_order(closure))
    candidates.append(longest_reference_chain_candidate_from_order(closure))
    for index in range(5):
        candidates.append(
            random_reference_chain_candidate_from_order(
                closure,
                seed=seed + index,
                name=f"random_reference_chain_{index}",
            )
        )
    reports = [
        evaluate_reference_chain_candidate(network, closure, candidate)
        for candidate in candidates
    ]
    ranked = rank_reference_chain_candidates(reports)
    return candidates, ranked


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run bracket coverage diagnostics over trigger probabilities."""

    rows: list[dict[str, float | str]] = []
    for trigger_probability in config.trigger_probability_values:
        for repetition in range(config.repetitions):
            seed = config.seed + int(10000 * trigger_probability) + repetition
            network = generate_state_change_network(
                config.num_systems,
                config.max_events,
                trigger_probability=trigger_probability,
                max_triggers_per_event=config.max_triggers_per_event,
                seed=seed,
            )
            adjacency = immediate_trigger_adjacency(network)
            closure = transitive_closure_dag(adjacency)
            candidates, ranked = _rank_candidates(network, closure, seed + 8000)
            by_name = {candidate.name: candidate for candidate in candidates}
            local_rows = [row for row in ranked if row["source"] == "local_system"]
            selected_rows = [("highest_utility_overall", ranked[0])]
            if local_rows:
                selected_rows.append(("highest_utility_local_system", local_rows[0]))
            relation_density = float(
                np.count_nonzero(closure) / (closure.shape[0] * (closure.shape[0] - 1))
            )
            for label, rank_row in selected_rows:
                candidate = by_name[str(rank_row["name"])]
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
                    "num_systems": float(config.num_systems),
                    "max_events": float(config.max_events),
                    "trigger_probability": float(trigger_probability),
                    "max_triggers_per_event": float(config.max_triggers_per_event),
                    "repetition": float(repetition),
                    "selection_label": label,
                    "candidate_source": candidate.source,
                    "utility_score": float(rank_row["score"]),
                    "relation_density": relation_density,
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
    """Write trigger-density bracket coverage rows."""

    output_path = (
        output_dir
        / "data"
        / "state_change_reference_bracket_coverage_vs_trigger_density.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, ...]:
    """Save trigger-density bracket coverage figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "state_change_reference_bracket_accessibility_vs_trigger_density.png",
            "accessible_fraction",
            "Accessible fraction",
        ),
        (
            "state_change_reference_bracket_width_vs_trigger_density.png",
            "mean_bracket_width_rank",
            "Mean bracket-width rank",
        ),
        (
            "state_change_reference_rank_slices_vs_trigger_density.png",
            "slice_count",
            "Rank slice count",
        ),
    )
    paths: list[Path] = []
    probabilities = sorted({float(row["trigger_probability"]) for row in rows})
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for label in sorted({str(row["selection_label"]) for row in rows}):
            means = []
            for probability in probabilities:
                values = [
                    float(row[key])
                    for row in rows
                    if float(row["trigger_probability"]) == probability
                    and row["selection_label"] == label
                ]
                means.append(float(np.mean(values)) if values else float("nan"))
            ax.plot(probabilities, means, marker="o", label=label)
        ax.set_xlabel("External trigger probability")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend()
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
    print(f"Wrote state-change bracket coverage versus trigger density: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
