"""Exact sanity checks for return-spectrum stability utilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_coarse_graining import (
    expected_coarse_delay_rank_for_motif,
    protected_indices_for_reference_and_motifs,
    remap_echo_motif_record_for_event_thinning,
    remap_reference_chain,
    restrict_transitive_order_to_retained_events,
    return_spectrum_stability_report,
    sample_retained_indices,
    spectrum_jaccard,
    subsample_reference_chain_positions,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_edge_coarse_graining import (
    protected_edge_keys_for_motifs,
    protected_source_target_pairs_for_reference_chain,
    thin_immediate_trigger_edges,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_return_spectrum_stability_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo spectrum stability sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_return_spectrum_stability_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact return-spectrum stability sanity checks."""

    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(network, reference, EchoMotifSpec(2, 3))
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    protected = protected_indices_for_reference_and_motifs(reference, [motif])
    retained = sample_retained_indices(
        len(network.events),
        keep_probability=0.0,
        protected_indices=protected,
        seed=0,
    )
    thinning = restrict_transitive_order_to_retained_events(closure, retained)
    coarse_reference = remap_reference_chain(reference, thinning.old_to_new)
    coarse_motif = remap_echo_motif_record_for_event_thinning(
        motif,
        thinning.old_to_new,
    )
    subsampling = subsample_reference_chain_positions(
        reference,
        stride=3,
        protected_positions=np.asarray([motif.emission_position]),
    )
    expected_delay = expected_coarse_delay_rank_for_motif(motif, subsampling)
    protected_pairs = protected_source_target_pairs_for_reference_chain(reference)
    protected_pairs |= protected_edge_keys_for_motifs([motif])
    thinned_network, _ = thin_immediate_trigger_edges(
        network,
        removal_probability=1.0,
        protected_source_target_pairs=protected_pairs,
        seed=0,
    )
    topological_order = topological_order_from_adjacency(
        immediate_trigger_adjacency(thinned_network)
    )
    relation_preserved = thinning.restricted_order_matrix[
        thinning.old_to_new[int(reference[2])],
        thinning.old_to_new[int(motif.target_event_id)],
    ]
    stability = return_spectrum_stability_report(
        np.asarray([3, 4, 5]),
        np.asarray([3, 4, 5]),
    )
    return [
        _row("protected_includes_reference", set(reference).issubset(protected), 1.0),
        _row("protected_includes_target", motif.target_event_id in protected, 1.0),
        _row(
            "sample_keeps_protected",
            set(protected).issubset(set(retained)),
            str(retained),
        ),
        _row("restricted_relation_preserved", bool(relation_preserved), 1.0),
        _row(
            "remap_reference_chain_size",
            coarse_reference.size == reference.size,
            1.0,
        ),
        _row("remap_motif_retained", coarse_motif is not None, 1.0),
        _row(
            "spectrum_jaccard_identical",
            spectrum_jaccard([1, 2], [2, 1]) == 1.0,
            1.0,
        ),
        _row("spectrum_jaccard_empty", spectrum_jaccard([], []) == 1.0, 1.0),
        _row(
            "stability_report_fields",
            stability["exact_earliest_preserved"] == 1.0,
            str(sorted(stability)),
        ),
        _row(
            "subsampling_keeps_protected",
            motif.emission_position in subsampling.retained_reference_positions,
            str(subsampling.retained_reference_positions),
        ),
        _row("expected_coarse_delay", expected_delay == 2, float(expected_delay)),
        _row(
            "edge_thinning_acyclic",
            topological_order.size == len(network.events),
            1.0,
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
    print(f"Wrote echo return-spectrum stability exact sanity: {output_path}")


if __name__ == "__main__":
    main()
