"""Local system chain diagnostic for state-change trigger networks."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_validation import local_chain_lengths

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for local reference-chain diagnostics."""

    num_systems: int = 8
    max_events: int = 300
    trigger_probability: float = 0.25
    max_triggers_per_event: int = 2
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="State-change reference-chain diagnostic."
    )
    parser.add_argument("--num-systems", type=int, default=8)
    parser.add_argument("--max-events", type=int, default=300)
    parser.add_argument("--trigger-probability", type=float, default=0.25)
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability=args.trigger_probability,
        max_triggers_per_event=args.max_triggers_per_event,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def _chain_is_totally_ordered(chain: list[int], closure: np.ndarray) -> bool:
    for left_index, left in enumerate(chain):
        for right in chain[left_index + 1 :]:
            if not closure[int(left), int(right)]:
                return False
    return True


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run local reference-chain diagnostics."""

    network = generate_state_change_network(
        config.num_systems,
        config.max_events,
        trigger_probability=config.trigger_probability,
        max_triggers_per_event=config.max_triggers_per_event,
        seed=config.seed,
    )
    adjacency = immediate_trigger_adjacency(network)
    closure = transitive_closure_dag(adjacency)
    chain_lengths = local_chain_lengths(network)
    observer_system_id = max(
        chain_lengths,
        key=lambda system_id: chain_lengths[system_id],
    )
    observer_chain = network.system_event_ids[observer_system_id]
    chain_array = np.asarray(observer_chain, dtype=int)
    comparable_to_chain = np.zeros(len(network.events), dtype=bool)
    both_bracketed = np.zeros(len(network.events), dtype=bool)
    for event_id in range(len(network.events)):
        has_predecessor_tick = bool(np.any(closure[chain_array, event_id]))
        has_successor_tick = bool(np.any(closure[event_id, chain_array]))
        on_chain = event_id in observer_chain
        comparable_to_chain[event_id] = (
            has_predecessor_tick or has_successor_tick or on_chain
        )
        both_bracketed[event_id] = has_predecessor_tick and has_successor_tick
    return [
        {
            "observer_system_id": float(observer_system_id),
            "observer_chain_length": float(len(observer_chain)),
            "total_event_count": float(len(network.events)),
            "chain_is_totally_ordered": float(
                _chain_is_totally_ordered(observer_chain, closure)
            ),
            "fraction_of_events_comparable_to_chain": float(
                np.mean(comparable_to_chain)
            ),
            "events_with_both_predecessor_and_successor_ticks": float(
                np.count_nonzero(both_bracketed)
            ),
        }
    ]


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write local reference-chain diagnostics."""

    output_path = output_dir / "data" / "state_change_observer_chain_diagnostic.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def save_figure(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Save local reference-chain coverage figure."""

    output_path = output_dir / "figures" / "state_change_observer_chain_coverage.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    row = rows[0]
    labels = ["comparable", "both brackets"]
    values = [
        row["fraction_of_events_comparable_to_chain"],
        row["events_with_both_predecessor_and_successor_ticks"]
        / row["total_event_count"],
    ]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.bar(labels, values)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Fraction of events")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_path = save_figure(rows, config.output_dir)
    print(f"Wrote state-change observer chain diagnostic: {data_path}")
    print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
