"""Exact sanity checks for controlled echo-response motifs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
    motif_recovery_report,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_selection import (
    select_echo_reachable_targets,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_motif_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo motif exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_motif_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact controlled echo-response motif checks."""

    network, reference = generate_reference_backbone_network(8)
    network, motif_a = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(
            emission_position=2,
            planted_delay_rank=3,
            outward_steps=1,
            return_steps=1,
        ),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    report_a = motif_recovery_report(closure, reference, motif_a)
    reachable = select_echo_reachable_targets(
        closure,
        reference,
        emission_position=2,
        target_indices=np.asarray([motif_a.target_event_id], dtype=int),
    )

    network, motif_b = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(
            emission_position=1,
            planted_delay_rank=2,
            outward_steps=0,
            return_steps=0,
        ),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    report_b = motif_recovery_report(closure, reference, motif_b)
    order_report = motif_order_recovery_rate([report_a, report_b])
    return [
        _row(
            "target_reachable",
            np.array_equal(reachable, [motif_a.target_event_id]),
            1.0,
        ),
        _row(
            "target_returns_to_r5",
            motif_a.reference_return_event_id == reference[5],
            1.0,
        ),
        _row("recovered_delay_a", report_a["recovered_delay_rank"] == 3.0, 3.0),
        _row(
            "exact_recovery_a",
            report_a["exact_recovery"] == 1.0,
            report_a["exact_recovery"],
        ),
        _row("recovered_delay_b", report_b["recovered_delay_rank"] == 2.0, 2.0),
        _row(
            "motif_order_recovery",
            order_report["order_inversion_rate"] == 0.0,
            order_report["order_inversion_rate"],
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
    print(f"Wrote echo motif exact sanity: {output_path}")


if __name__ == "__main__":
    main()
