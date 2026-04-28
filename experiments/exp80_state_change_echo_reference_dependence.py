"""Reference-protocol dependence for same-emission echo order."""

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

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for echo-order reference dependence."""

    num_systems: int = 10
    max_events: int = 500
    trigger_probability_values: tuple[float, ...] = (0.10, 0.20, 0.35)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    top_k: int = 5
    emission_strategy: str = "middle"
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="State-change echo-order reference dependence."
    )
    parser.add_argument("--num-systems", type=int, default=10)
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument(
        "--trigger-probability-values",
        nargs="+",
        default=["0.10", "0.20", "0.35"],
    )
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--emission-strategy", default="middle")
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
        top_k=args.top_k,
        emission_strategy=args.emission_strategy,
        output_dir=args.output_dir,
    )


def generate_reference_candidates(
    order_matrix: np.ndarray,
    network,
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


def _select_emission_position(
    candidate: ReferenceChainCandidate,
    strategy: str,
) -> int | None:
    if strategy == "middle":
        positions = emission_positions_for_reference_chain(
            candidate.chain_event_ids,
            strategy="early_middle_late",
            count=3,
        )
        if positions.size == 0:
            return None
        return int(positions[positions.size // 2])
    positions = emission_positions_for_reference_chain(
        candidate.chain_event_ids,
        strategy=strategy,
        count=1,
    )
    return int(positions[0]) if positions.size else None


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run pairwise echo-order reference-dependence diagnostics."""

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
            closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            candidates = generate_reference_candidates(
                closure,
                network,
                random_candidate_count=5,
                seed=seed + 8000,
            )
            reports = [
                evaluate_reference_chain_candidate(network, closure, candidate)
                for candidate in candidates
            ]
            ranked_rows = rank_reference_chain_candidates(reports)
            by_name = {candidate.name: candidate for candidate in candidates}
            computed: list[dict[str, object]] = []
            for rank_row in ranked_rows[: config.top_k]:
                candidate = by_name[str(rank_row["name"])]
                emission_position = _select_emission_position(
                    candidate,
                    config.emission_strategy,
                )
                if emission_position is None:
                    continue
                _, delays, reachable = echo_reachability_from_emission(
                    closure,
                    candidate.chain_event_ids,
                    emission_position,
                )
                computed.append(
                    {
                        "candidate": candidate,
                        "rank_row": rank_row,
                        "delay_ranks": delays,
                        "reachable": reachable,
                        "emission_position": float(emission_position),
                    }
                )
            for i, left in enumerate(computed):
                for j, right in enumerate(computed[i + 1 :], start=i + 1):
                    reachable_left = np.asarray(left["reachable"], dtype=bool)
                    reachable_right = np.asarray(right["reachable"], dtype=bool)
                    overlap_denominator = int(
                        np.count_nonzero(reachable_left | reachable_right)
                    )
                    overlap = (
                        float(
                            np.count_nonzero(reachable_left & reachable_right)
                            / overlap_denominator
                        )
                        if overlap_denominator
                        else float("nan")
                    )
                    comparison = compare_echo_orders(
                        left["delay_ranks"],
                        right["delay_ranks"],
                        reachable_left,
                        reachable_right,
                    )
                    left_candidate = left["candidate"]
                    right_candidate = right["candidate"]
                    left_rank = left["rank_row"]
                    right_rank = right["rank_row"]
                    rows.append(
                        {
                            "num_systems": float(config.num_systems),
                            "max_events": float(config.max_events),
                            "trigger_probability": float(trigger_probability),
                            "max_triggers_per_event": float(
                                config.max_triggers_per_event
                            ),
                            "repetition": float(repetition),
                            "left_rank": float(i + 1),
                            "right_rank": float(j + 1),
                            "left_source": left_candidate.source,
                            "right_source": right_candidate.source,
                            "left_score": float(left_rank["score"]),
                            "right_score": float(right_rank["score"]),
                            "utility_score_gap": abs(
                                float(left_rank["score"]) - float(right_rank["score"])
                            ),
                            "left_emission_position": float(
                                left["emission_position"]
                            ),
                            "right_emission_position": float(
                                right["emission_position"]
                            ),
                            "reachable_set_overlap": overlap,
                            "common_reachable_count": comparison[
                                "common_reachable_count"
                            ],
                            "common_reachable_fraction": comparison[
                                "common_reachable_fraction"
                            ],
                            "echo_order_inversion_rate": comparison[
                                "order_inversion_rate"
                            ],
                            "echo_order_agreement_rate": comparison[
                                "order_agreement_rate"
                            ],
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write echo-order reference-dependence rows."""

    output_path = output_dir / "data" / "state_change_echo_reference_dependence.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else ["empty"]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, ...]:
    """Save echo-order reference-dependence figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "state_change_echo_order_disagreement.png",
            "echo_order_inversion_rate",
            "Echo-order inversion",
        ),
        (
            "state_change_echo_reachable_overlap.png",
            "reachable_set_overlap",
            "Reachable set overlap",
        ),
    )
    paths: list[Path] = []
    probabilities = sorted({float(row["trigger_probability"]) for row in rows})
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        means = []
        for probability in probabilities:
            values = [
                float(row[key])
                for row in rows
                if float(row["trigger_probability"]) == probability
                and np.isfinite(float(row[key]))
            ]
            means.append(float(np.mean(values)) if values else float("nan"))
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(probabilities, means, marker="o")
        ax.set_xlabel("External trigger probability")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)

    gap_path = figure_dir / "state_change_utility_gap_vs_echo_disagreement.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    xs = [float(row["utility_score_gap"]) for row in rows]
    ys = [
        float(row["echo_order_inversion_rate"])
        if np.isfinite(float(row["echo_order_inversion_rate"]))
        else np.nan
        for row in rows
    ]
    ax.scatter(xs, ys, alpha=0.7)
    ax.set_xlabel("Reference-chain utility-score gap")
    ax.set_ylabel("Echo-order inversion")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(gap_path, dpi=200)
    plt.close(fig)
    paths.append(gap_path)
    return tuple(paths)


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote state-change echo reference dependence: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
