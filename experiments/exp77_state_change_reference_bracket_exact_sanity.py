"""Exact sanity checks for state-change reference-chain brackets."""

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
from causal_spacetime_lab.state_change_bracket_diagnostics import (
    bracket_coverage_summary,
)
from causal_spacetime_lab.state_change_brackets import (
    assign_reference_rank_slices,
    bracket_width_rank_from_reference_brackets,
    earliest_successor_reference_position,
    latest_predecessor_reference_position,
    radar_time_rank_from_reference_brackets,
    reference_tick_brackets_from_order,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/state_change_reference_bracket_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Reference bracket exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "state_change_reference_bracket_exact_sanity.csv"


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
    return StateChangeNetwork(events, edges, {0: [0, 2], 1: [1, 3, 4]})


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact bracket sanity checks."""

    network = build_hand_coded_network()
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    reference_chain = np.asarray([1, 3, 4], dtype=int)
    predecessors, successors, accessible = reference_tick_brackets_from_order(
        closure,
        reference_chain,
    )
    time_ranks = radar_time_rank_from_reference_brackets(
        predecessors,
        successors,
        accessible,
    )
    width_ranks = bracket_width_rank_from_reference_brackets(
        predecessors,
        successors,
        accessible,
    )
    slices = assign_reference_rank_slices(time_ranks, accessible, bin_width=2)
    summary = bracket_coverage_summary(
        predecessors,
        successors,
        accessible,
        reference_chain_length=reference_chain.size,
    )
    return [
        _row(
            "latest_predecessor",
            latest_predecessor_reference_position(closure, reference_chain, 2) == 0,
            0.0,
        ),
        _row(
            "earliest_successor",
            earliest_successor_reference_position(closure, reference_chain, 2) == 2,
            2.0,
        ),
        _row(
            "predecessor_array",
            np.array_equal(predecessors, np.asarray([-1, -1, 0, -1, -1])),
            ",".join(str(value) for value in predecessors),
        ),
        _row(
            "successor_array",
            np.array_equal(successors, np.asarray([2, -1, 2, -1, -1])),
            ",".join(str(value) for value in successors),
        ),
        _row(
            "accessible_mask",
            np.array_equal(
                accessible,
                np.asarray([False, False, True, False, False]),
            ),
            "".join("1" if value else "0" for value in accessible),
        ),
        _row(
            "radar_time_rank",
            np.array_equal(time_ranks, np.asarray([-1, -1, 2, -1, -1])),
            ",".join(str(value) for value in time_ranks),
        ),
        _row(
            "bracket_width_rank",
            np.array_equal(width_ranks, np.asarray([-1, -1, 2, -1, -1])),
            ",".join(str(value) for value in width_ranks),
        ),
        _row(
            "slice_labels",
            np.array_equal(slices, np.asarray([-1, -1, 1, -1, -1])),
            ",".join(str(value) for value in slices),
        ),
        _row(
            "bracket_summary",
            summary["accessible_count"] == 1.0
            and summary["successor_only_count"] == 1.0,
            summary["accessible_fraction"],
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
    print(f"Wrote state-change reference bracket exact sanity: {output_path}")


if __name__ == "__main__":
    main()
