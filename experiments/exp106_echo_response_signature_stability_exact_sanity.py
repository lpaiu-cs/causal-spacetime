"""Exact sanity checks for response-signature stability utilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_motifs import insert_multiple_echo_motifs
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_response_layers import (
    EchoResponseLayerSpec,
    build_layered_echo_motif_specs,
    planted_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_response_signature import (
    echo_response_signature_from_motif_rows,
    response_order_sign_matrix,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    consensus_response_order_matrix,
    response_order_cycle_count,
    stable_response_order_core,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_response_signature_stability_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Echo response signature stability sanity."
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return (
        args.output_dir
        / "data"
        / "echo_response_signature_stability_exact_sanity.csv"
    )


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact response-signature stability checks."""

    signs = response_order_sign_matrix([2, 5, -1], [True, True, False])
    rows = [
        {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
        {"target_event_id": 2.0, "recovered_delay_rank": 5.0},
    ]
    signature = echo_response_signature_from_motif_rows(rows)
    identical = compare_response_signatures(signature, signature)
    variant_a = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 5.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 8.0},
        ],
        label="a",
    )
    variant_b = echo_response_signature_from_motif_rows(
        [
            {"target_event_id": 1.0, "recovered_delay_rank": 2.0},
            {"target_event_id": 2.0, "recovered_delay_rank": 6.0},
            {"target_event_id": 3.0, "recovered_delay_rank": 8.0},
        ],
        label="b",
    )
    core = stable_response_order_core([variant_a, variant_b])
    consensus = consensus_response_order_matrix([variant_a, variant_b])
    cyclic = np.asarray([[0, -1, 1], [1, 0, -1], [-1, 1, 0]], dtype=int)
    acyclic = np.asarray([[0, -1, -1], [1, 0, -1], [1, 1, 0]], dtype=int)
    network, reference = generate_reference_backbone_network(12)
    specs = build_layered_echo_motif_specs(
        reference,
        2,
        [
            EchoResponseLayerSpec(3, 2),
            EchoResponseLayerSpec(5, 2),
        ],
        seed=0,
    )
    network, motifs = insert_multiple_echo_motifs(network, reference, specs)
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    planted = planted_response_signature_from_motifs(motifs)
    return [
        _row("response_order_sign_matrix", signs[0, 1] == -1, str(signs)),
        _row(
            "signature_from_rows",
            np.array_equal(signature.delay_ranks, np.asarray([2, 5])),
            str(signature.delay_ranks),
        ),
        _row(
            "compare_identical",
            identical["pair_agreement_fraction"] == 1.0,
            identical["pair_agreement_fraction"],
        ),
        _row(
            "stable_core",
            float(core["stable_pair_fraction"]) == 1.0,
            float(core["stable_pair_fraction"]),
        ),
        _row(
            "consensus_matrix",
            np.array_equal(consensus, variant_a.order_sign_matrix),
            str(consensus),
        ),
        _row("cycle_count_cyclic", response_order_cycle_count(cyclic) == 1, 1.0),
        _row("cycle_count_acyclic", response_order_cycle_count(acyclic) == 0, 0.0),
        _row("layered_spec_count", len(specs) == 4, float(len(specs))),
        _row(
            "planted_signature_clean",
            planted.target_event_ids.size == 4
            and np.array_equal(np.sort(planted.delay_ranks), np.asarray([3, 3, 5, 5])),
            str(planted.delay_ranks),
        ),
        _row("clean_closure_exists", closure.shape[0] == len(network.events), 1.0),
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
    print(f"Wrote echo response signature stability exact sanity: {output_path}")


if __name__ == "__main__":
    main()
