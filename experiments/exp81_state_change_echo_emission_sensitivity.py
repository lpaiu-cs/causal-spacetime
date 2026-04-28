"""Emission-position sensitivity for same-emission echo order."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_echo import (
    compare_echo_orders,
    echo_reachability_from_emission,
)
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
    """Configuration for echo emission-position sensitivity."""

    num_systems: int = 10
    max_events: int = 500
    trigger_probability: float = 0.25
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    emission_count: int = 5
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="State-change echo emission-position sensitivity."
    )
    parser.add_argument("--num-systems", type=int, default=10)
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument("--trigger-probability", type=float, default=0.25)
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--emission-count", type=int, default=5)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability=args.trigger_probability,
        max_triggers_per_event=args.max_triggers_per_event,
        repetitions=args.repetitions,
        seed=args.seed,
        emission_count=args.emission_count,
        output_dir=args.output_dir,
    )


def _generate_candidates(order_matrix: np.ndarray, network, seed: int):
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


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run emission-position sensitivity diagnostics."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        seed = config.seed + repetition
        network = generate_state_change_network(
            config.num_systems,
            config.max_events,
            trigger_probability=config.trigger_probability,
            max_triggers_per_event=config.max_triggers_per_event,
            seed=seed,
        )
        closure = transitive_closure_dag(immediate_trigger_adjacency(network))
        candidates = _generate_candidates(closure, network, seed + 9000)
        reports = [
            evaluate_reference_chain_candidate(network, closure, candidate)
            for candidate in candidates
        ]
        ranked_rows = rank_reference_chain_candidates(reports)
        by_name = {candidate.name: candidate for candidate in candidates}
        best_row = ranked_rows[0]
        candidate = by_name[str(best_row["name"])]
        emissions = emission_positions_for_reference_chain(
            candidate.chain_event_ids,
            strategy="interior_quantiles",
            count=config.emission_count,
        )
        computed: list[dict[str, object]] = []
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
            relative_position = float(
                emission_position / max(candidate.chain_event_ids.size - 1, 1)
            )
            row: dict[str, float | str] = {
                "row_type": "emission",
                "num_systems": float(config.num_systems),
                "max_events": float(config.max_events),
                "trigger_probability": float(config.trigger_probability),
                "max_triggers_per_event": float(config.max_triggers_per_event),
                "repetition": float(repetition),
                "candidate_name": candidate.name,
                "candidate_source": candidate.source,
                "utility_score": float(best_row["score"]),
                "chain_length": float(candidate.chain_event_ids.size),
                "left_emission_position": float(emission_position),
                "right_emission_position": float("nan"),
                "relative_emission_position": relative_position,
                "pairwise_echo_order_inversion_rate": float("nan"),
                "pairwise_common_reachable_count": float("nan"),
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
            computed.append(
                {
                    "emission_position": float(emission_position),
                    "delay_ranks": delays,
                    "reachable": reachable,
                }
            )
        for i, left in enumerate(computed):
            for right in computed[i + 1 :]:
                comparison = compare_echo_orders(
                    left["delay_ranks"],
                    right["delay_ranks"],
                    left["reachable"],
                    right["reachable"],
                )
                rows.append(
                    {
                        "row_type": "pairwise_emission",
                        "num_systems": float(config.num_systems),
                        "max_events": float(config.max_events),
                        "trigger_probability": float(config.trigger_probability),
                        "max_triggers_per_event": float(config.max_triggers_per_event),
                        "repetition": float(repetition),
                        "candidate_name": candidate.name,
                        "candidate_source": candidate.source,
                        "utility_score": float(best_row["score"]),
                        "chain_length": float(candidate.chain_event_ids.size),
                        "left_emission_position": float(left["emission_position"]),
                        "right_emission_position": float(right["emission_position"]),
                        "relative_emission_position": float("nan"),
                        "target_count": float("nan"),
                        "reachable_count": float("nan"),
                        "reachable_fraction": float("nan"),
                        "emission_position": float("nan"),
                        "reference_chain_length": float(candidate.chain_event_ids.size),
                        "mean_return_position": float("nan"),
                        "median_return_position": float("nan"),
                        "mean_echo_delay_rank": float("nan"),
                        "median_echo_delay_rank": float("nan"),
                        "max_echo_delay_rank": float("nan"),
                        "singleton_delay_fraction": float("nan"),
                        "distinct_delay_rank_count": float("nan"),
                        "tied_pair_fraction": float("nan"),
                        "strict_order_pair_fraction": float("nan"),
                        "pairwise_echo_order_inversion_rate": comparison[
                            "order_inversion_rate"
                        ],
                        "pairwise_common_reachable_count": comparison[
                            "common_reachable_count"
                        ],
                    }
                )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write emission sensitivity rows."""

    output_path = output_dir / "data" / "state_change_echo_emission_sensitivity.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else ["empty"]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_relative_position(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[float], list[float]]:
    emission_rows = [row for row in rows if row["row_type"] == "emission"]
    positions = sorted(
        {float(row["relative_emission_position"]) for row in emission_rows}
    )
    means = []
    for position in positions:
        values = [
            float(row[key])
            for row in emission_rows
            if float(row["relative_emission_position"]) == position
            and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return positions, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, ...]:
    """Save emission sensitivity figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    specs = (
        (
            "state_change_echo_reachability_by_emission_position.png",
            "reachable_fraction",
            "Reachable fraction",
        ),
        (
            "state_change_echo_delay_by_emission_position.png",
            "mean_echo_delay_rank",
            "Mean echo-delay rank",
        ),
    )
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        xs, ys = _mean_by_relative_position(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel("Relative emission position")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)

    pair_path = figure_dir / "state_change_echo_emission_order_disagreement.png"
    pair_rows = [row for row in rows if row["row_type"] == "pairwise_emission"]
    values = [
        float(row["pairwise_echo_order_inversion_rate"])
        for row in pair_rows
        if np.isfinite(float(row["pairwise_echo_order_inversion_rate"]))
    ]
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.bar(["pairwise emissions"], [float(np.mean(values)) if values else float("nan")])
    ax.set_ylabel("Mean echo-order inversion")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(pair_path, dpi=200)
    plt.close(fig)
    paths.append(pair_path)
    return tuple(paths)


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote state-change echo emission sensitivity: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
