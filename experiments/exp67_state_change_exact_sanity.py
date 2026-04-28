"""Exact sanity checks for finite state-change causal trigger networks."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    StateChangeNetwork,
    TriggerEdge,
)
from causal_spacetime_lab.state_change_order import (
    causal_interval_indices,
    immediate_trigger_adjacency,
    is_irreflexive,
    is_transitive,
    topological_order_from_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_validation import local_chain_lengths

DEFAULT_OUTPUT = Path("outputs/data/state_change_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="State-change exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "state_change_exact_sanity.csv"


def build_hand_coded_network() -> StateChangeNetwork:
    """Build a small deterministic trigger network."""

    events = [
        StateChangeEvent(0, 0, 0, 0, 1),
        StateChangeEvent(1, 1, 0, 0, 1),
        StateChangeEvent(2, 0, 1, 1, 2),
        StateChangeEvent(3, 1, 1, 1, 2),
        StateChangeEvent(4, 1, 2, 2, 3),
    ]
    edges = [
        TriggerEdge(-1, 0, "initial_seed"),
        TriggerEdge(-1, 1, "initial_seed"),
        TriggerEdge(0, 2, "local_successor"),
        TriggerEdge(1, 2, "external_trigger"),
        TriggerEdge(1, 3, "local_successor"),
        TriggerEdge(3, 4, "local_successor"),
        TriggerEdge(2, 4, "external_trigger"),
    ]
    return StateChangeNetwork(
        events=events,
        trigger_edges=edges,
        system_event_ids={0: [0, 2], 1: [1, 3, 4]},
    )


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact sanity checks."""

    network = build_hand_coded_network()
    adjacency = immediate_trigger_adjacency(network)
    closure = transitive_closure_dag(adjacency)
    interval = causal_interval_indices(closure, 1, 4)
    chain_lengths = local_chain_lengths(network)
    topological_order = topological_order_from_adjacency(adjacency)
    return [
        _row("irreflexive", is_irreflexive(closure), 0.0),
        _row("transitive", is_transitive(closure), 0.0),
        _row(
            "expected_interval",
            np.array_equal(interval, np.asarray([2, 3])),
            ",".join(str(value) for value in interval),
        ),
        _row(
            "local_chain_lengths",
            chain_lengths == {0: 2, 1: 3},
            str(chain_lengths),
        ),
        _row(
            "topological_order_exists",
            topological_order.size == len(network.events),
            topological_order.size,
        ),
    ]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Write exact sanity rows."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    output_path = parse_args()
    rows = run_experiment()
    output_path = write_outputs(rows, output_path)
    print(f"Wrote state-change exact sanity data: {output_path}")


if __name__ == "__main__":
    main()
