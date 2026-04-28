"""Exact sanity checks for observer-like chain selection utilities."""

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
from causal_spacetime_lab.state_change_observer_quality import (
    evaluate_chain_candidate,
    rank_chain_candidates,
)
from causal_spacetime_lab.state_change_observers import (
    ObserverChainCandidate,
    chain_bracketing_mask,
    chain_comparability_mask,
    earliest_successor_chain_position,
    is_chain,
    latest_predecessor_chain_position,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/observer_chain_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Observer-chain exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "observer_chain_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


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


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact observer-like chain sanity checks."""

    network = build_hand_coded_network()
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    candidate = ObserverChainCandidate(
        name="manual_system_1",
        chain_event_ids=np.asarray([1, 3, 4]),
        source="manual",
    )
    comparable = chain_comparability_mask(closure, candidate.chain_event_ids)
    bracketed = chain_bracketing_mask(closure, candidate.chain_event_ids)
    report = evaluate_chain_candidate(network, closure, candidate)
    ranked = rank_chain_candidates([report])
    return [
        _row("is_chain", is_chain(closure, candidate.chain_event_ids), 0.0),
        _row(
            "comparability_mask",
            np.array_equal(comparable, np.asarray([True, True, True, True, True])),
            "".join("1" if value else "0" for value in comparable),
        ),
        _row(
            "bracketing_mask",
            np.array_equal(bracketed, np.asarray([False, False, True, False, False])),
            "".join("1" if value else "0" for value in bracketed),
        ),
        _row(
            "latest_predecessor",
            latest_predecessor_chain_position(closure, candidate.chain_event_ids, 2)
            == 0,
            0.0,
        ),
        _row(
            "earliest_successor",
            earliest_successor_chain_position(closure, candidate.chain_event_ids, 2)
            == 2,
            2.0,
        ),
        _row("quality_valid_chain", report.is_valid_chain, 0.0),
        _row(
            "quality_bracket_fraction",
            np.isclose(report.bracketed_fraction, 0.5),
            report.bracketed_fraction,
        ),
        _row("ranking_score_present", "score" in ranked[0], ranked[0]["score"]),
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
    print(f"Wrote observer-chain exact sanity data: {output_path}")


if __name__ == "__main__":
    main()
