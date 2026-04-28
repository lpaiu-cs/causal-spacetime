"""Interval-profile diagnostics for observer-like chain candidates."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
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
    causal_interval_indices,
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for chain interval-profile diagnostics."""

    num_systems: int = 8
    max_events: int = 400
    trigger_probability: float = 0.25
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    candidate_sources: tuple[str, ...] = (
        "best",
        "longest_local",
        "longest_order",
        "random",
    )
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Observer-chain interval-profile diagnostic."
    )
    parser.add_argument("--num-systems", type=int, default=8)
    parser.add_argument("--max-events", type=int, default=400)
    parser.add_argument("--trigger-probability", type=float, default=0.25)
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--candidate-sources",
        nargs="+",
        default=["best", "longest_local", "longest_order", "random"],
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability=args.trigger_probability,
        max_triggers_per_event=args.max_triggers_per_event,
        repetitions=args.repetitions,
        seed=args.seed,
        candidate_sources=tuple(args.candidate_sources),
        output_dir=args.output_dir,
    )


def _interval_sizes(
    order_matrix: np.ndarray,
    chain_event_ids: np.ndarray,
    window: int,
) -> np.ndarray:
    sizes: list[float] = []
    if chain_event_ids.size <= window:
        return np.empty(0, dtype=float)
    for left, right in zip(
        chain_event_ids[:-window],
        chain_event_ids[window:],
        strict=True,
    ):
        if order_matrix[int(left), int(right)]:
            sizes.append(
                float(causal_interval_indices(order_matrix, int(left), int(right)).size)
            )
    return np.asarray(sizes, dtype=float)


def _summarize_interval_sizes(prefix: str, values: np.ndarray) -> dict[str, float]:
    mean_value = float(np.mean(values)) if values.size else 0.0
    std_value = float(np.std(values)) if values.size else 0.0
    cv_value = float(std_value / mean_value) if mean_value > 0.0 else float("nan")
    return {
        f"{prefix}_mean": mean_value,
        f"{prefix}_std": std_value,
        f"{prefix}_cv": cv_value,
        f"{prefix}_max": float(np.max(values)) if values.size else 0.0,
        f"{prefix}_count": float(values.size),
    }


def _select_candidates(
    network,
    order_matrix: np.ndarray,
    seed: int,
) -> dict[str, ObserverChainCandidate]:
    local_candidates = local_system_chain_candidates(network)
    if not local_candidates:
        local_candidates = local_system_chain_candidates(network, min_length=1)
    longest_local = max(
        local_candidates,
        key=lambda candidate: candidate.chain_event_ids.size,
    )
    longest_order = longest_chain_candidate_from_order(order_matrix)
    random_candidate = random_chain_candidate_from_order(order_matrix, seed=seed)
    all_candidates = [
        *local_candidates,
        greedy_chain_candidate_from_order(order_matrix),
        longest_order,
        random_candidate,
    ]
    reports = [
        evaluate_chain_candidate(network, order_matrix, candidate)
        for candidate in all_candidates
    ]
    ranked = rank_chain_candidates(reports)
    best_name = str(ranked[0]["name"])
    best_candidate = next(
        candidate for candidate in all_candidates if candidate.name == best_name
    )
    return {
        "best": best_candidate,
        "longest_local": longest_local,
        "longest_order": longest_order,
        "random": random_candidate,
    }


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], dict[str, np.ndarray]]:
    """Run interval-profile diagnostics."""

    rows: list[dict[str, float | str]] = []
    example_profile: dict[str, np.ndarray] = {}
    for repetition in range(config.repetitions):
        run_seed = config.seed + repetition
        network = generate_state_change_network(
            config.num_systems,
            config.max_events,
            trigger_probability=config.trigger_probability,
            max_triggers_per_event=config.max_triggers_per_event,
            seed=run_seed,
        )
        closure = transitive_closure_dag(immediate_trigger_adjacency(network))
        candidates = _select_candidates(network, closure, run_seed + 9000)
        for source in config.candidate_sources:
            if source not in candidates:
                continue
            candidate = candidates[source]
            adjacent = _interval_sizes(closure, candidate.chain_event_ids, window=1)
            window_two = _interval_sizes(closure, candidate.chain_event_ids, window=2)
            if repetition == 0 and source == config.candidate_sources[0]:
                example_profile = {
                    "adjacent": adjacent,
                    "window_two": window_two,
                }
            row: dict[str, float | str] = {
                "num_systems": float(config.num_systems),
                "max_events": float(config.max_events),
                "trigger_probability": float(config.trigger_probability),
                "max_triggers_per_event": float(config.max_triggers_per_event),
                "repetition": float(repetition),
                "candidate_source": source,
                "candidate_name": candidate.name,
                "candidate_original_source": candidate.source,
                "chain_length": float(candidate.chain_event_ids.size),
            }
            row.update(_summarize_interval_sizes("adjacent_interval", adjacent))
            row.update(_summarize_interval_sizes("window2_interval", window_two))
            rows.append(row)
    return rows, example_profile


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write interval-profile diagnostics."""

    output_path = output_dir / "data" / "observer_chain_interval_profile.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figures(
    rows: list[dict[str, float | str]],
    example_profile: dict[str, np.ndarray],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save interval-profile figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    cv_path = figure_dir / "observer_chain_interval_cv_by_source.png"
    sources = sorted({str(row["candidate_source"]) for row in rows})
    means: list[float] = []
    for source in sources:
        values = np.asarray(
            [
                float(row["adjacent_interval_cv"])
                for row in rows
                if row["candidate_source"] == source
            ],
            dtype=float,
        )
        finite_values = values[np.isfinite(values)]
        means.append(float(np.mean(finite_values)) if finite_values.size else 0.0)
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.bar(sources, means)
    ax.set_ylabel("Mean adjacent interval CV")
    ax.set_xlabel("Candidate source")
    ax.grid(True, axis="y", alpha=0.3)
    fig.autofmt_xdate(rotation=20)
    fig.tight_layout()
    fig.savefig(cv_path, dpi=200)
    plt.close(fig)

    profile_path = figure_dir / "observer_chain_interval_profile_example.png"
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    adjacent = example_profile.get("adjacent", np.empty(0, dtype=float))
    window_two = example_profile.get("window_two", np.empty(0, dtype=float))
    if adjacent.size:
        ax.plot(np.arange(adjacent.size), adjacent, marker="o", label="adjacent")
    if window_two.size:
        ax.plot(np.arange(window_two.size), window_two, marker="s", label="window 2")
    ax.set_xlabel("Chain gap index")
    ax.set_ylabel("Causal interval cardinality")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(profile_path, dpi=200)
    plt.close(fig)
    return cv_path, profile_path


def main() -> None:
    config = parse_args()
    rows, example_profile = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, example_profile, config.output_dir)
    print(f"Wrote observer-chain interval profile data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
