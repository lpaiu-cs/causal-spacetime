"""Exact sanity checks for same-emission echo-order utilities."""

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
from causal_spacetime_lab.state_change_echo import (
    echo_delay_rank_for_emission,
    echo_order_matrix_from_delay_ranks,
    echo_reachability_from_emission,
    echo_return_position_for_emission,
)
from causal_spacetime_lab.state_change_echo_diagnostics import (
    echo_delay_histogram,
    echo_order_resolution_summary,
    echo_reachability_summary,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/state_change_echo_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Same-emission echo exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "state_change_echo_exact_sanity.csv"


def build_hand_coded_network() -> StateChangeNetwork:
    """Build a deterministic trigger network with one echo-reachable target."""

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
    """Run exact same-emission echo checks."""

    network = build_hand_coded_network()
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    reference_chain = np.asarray([1, 3, 4], dtype=int)
    emission_position = 0
    returns, delays, reachable = echo_reachability_from_emission(
        closure,
        reference_chain,
        emission_position,
    )
    order_matrix = echo_order_matrix_from_delay_ranks(delays, reachable)
    summary = echo_reachability_summary(
        returns,
        delays,
        reachable,
        emission_position=emission_position,
        reference_chain_length=reference_chain.size,
    )
    resolution = echo_order_resolution_summary(delays, reachable)
    histogram = echo_delay_histogram(delays, reachable)
    return [
        _row(
            "return_position",
            echo_return_position_for_emission(
                closure,
                reference_chain,
                target_index=2,
                emission_position=emission_position,
            )
            == 2,
            2.0,
        ),
        _row(
            "delay_rank",
            echo_delay_rank_for_emission(
                closure,
                reference_chain,
                target_index=2,
                emission_position=emission_position,
            )
            == 2,
            2.0,
        ),
        _row(
            "return_array",
            np.array_equal(returns, np.asarray([-1, -1, 2, -1, -1])),
            ",".join(str(value) for value in returns),
        ),
        _row(
            "delay_array",
            np.array_equal(delays, np.asarray([-1, -1, 2, -1, -1])),
            ",".join(str(value) for value in delays),
        ),
        _row(
            "reachable_mask",
            np.array_equal(
                reachable,
                np.asarray([False, False, True, False, False]),
            ),
            "".join("1" if value else "0" for value in reachable),
        ),
        _row("echo_order_matrix", not np.any(order_matrix), float(order_matrix.sum())),
        _row(
            "reachability_summary",
            summary["reachable_count"] == 1.0
            and summary["distinct_delay_rank_count"] == 1.0,
            summary["reachable_fraction"],
        ),
        _row(
            "delay_histogram",
            histogram == {2: 1},
            ",".join(f"{key}:{value}" for key, value in histogram.items()),
        ),
        _row(
            "resolution_summary",
            resolution["reachable_count"] == 1.0
            and resolution["comparable_pair_count"] == 0.0,
            resolution["comparable_pair_count"],
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
    print(f"Wrote state-change echo exact sanity: {output_path}")


if __name__ == "__main__":
    main()
