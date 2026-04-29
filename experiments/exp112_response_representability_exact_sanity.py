"""Exact sanity checks for response-order representability utilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_spectrum_semantics import (
    gated_echo_delay_from_spectrum,
)
from causal_spacetime_lab.state_change_gated_echo import (
    gated_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_response_representability import (
    has_response_order_cycle,
    response_order_topological_ranks,
    scalar_rank_representation_error,
    scalar_representability_report,
)

DEFAULT_OUTPUT = Path("outputs/data/response_representability_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Response representability sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "response_representability_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact scalar representability checks."""

    acyclic = np.asarray([[0, -1, -1], [1, 0, -1], [1, 1, 0]], dtype=int)
    cyclic = np.asarray([[0, -1, 1], [1, 0, -1], [-1, 1, 0]], dtype=int)
    ranks = response_order_topological_ranks(acyclic)
    good_error = scalar_rank_representation_error(acyclic, ranks)
    bad_error = scalar_rank_representation_error(acyclic, np.asarray([2, 1, 0]))
    clean_report = scalar_representability_report(acyclic)
    cyclic_report = scalar_representability_report(cyclic)
    network, reference = generate_reference_backbone_network(10)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 4, return_steps=0),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    gated = gated_response_signature_from_motifs(
        closure,
        reference,
        [motif],
        gate_delay_rank=3,
    )
    return [
        _row("acyclic_representable", not has_response_order_cycle(acyclic), 1.0),
        _row("cyclic_not_representable", has_response_order_cycle(cyclic), 1.0),
        _row(
            "topological_ranks_satisfy",
            good_error["violation_count"] == 0.0,
            str(ranks),
        ),
        _row(
            "representation_error_detects_violation",
            bad_error["violation_count"] > 0.0,
            bad_error["violation_count"],
        ),
        _row(
            "scalar_report_clean",
            clean_report["scalar_representable"] == 1.0,
            clean_report["rank_span"],
        ),
        _row(
            "scalar_report_cyclic",
            cyclic_report["scalar_representable"] == 0.0,
            cyclic_report["directed_3cycle_count"],
        ),
        _row(
            "gated_delay_examples",
            gated_echo_delay_from_spectrum([2, 4, 5], 3) == 4
            and gated_echo_delay_from_spectrum([2], 3) is None,
            1.0,
        ),
        _row(
            "gated_response_signature",
            gated.reachable_mask.size == 1
            and bool(gated.reachable_mask[0])
            and int(gated.delay_ranks[0]) == 4,
            str(gated.delay_ranks),
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
    output_path = write_outputs(run_experiment(), parse_args())
    print(f"Wrote response representability exact sanity: {output_path}")


if __name__ == "__main__":
    main()
