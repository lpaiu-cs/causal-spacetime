"""Exact sanity checks for echo coarse-graining utilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_coarse_graining import (
    expected_coarse_delay_rank_for_motif,
    protected_indices_for_reference_and_motifs,
    remap_echo_motif_record_for_event_thinning,
    remap_reference_chain,
    restrict_transitive_order_to_retained_events,
    return_spectrum_stability_report,
    subsample_reference_chain_positions,
)
from causal_spacetime_lab.state_change_echo_interference import return_delay_spectrum
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_edge_coarse_graining import (
    compare_recovery_before_after_edge_thinning,
    thin_immediate_trigger_edges,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_coarse_graining_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo coarse-graining exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_coarse_graining_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact coarse-graining sanity checks."""

    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    baseline_spectrum = return_delay_spectrum(
        closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    protected = protected_indices_for_reference_and_motifs(reference, [motif])
    thinning = restrict_transitive_order_to_retained_events(closure, protected)
    coarse_reference = remap_reference_chain(reference, thinning.old_to_new)
    coarse_motif = remap_echo_motif_record_for_event_thinning(
        motif,
        thinning.old_to_new,
    )
    assert coarse_motif is not None
    coarse_spectrum = return_delay_spectrum(
        thinning.restricted_order_matrix,
        coarse_reference,
        coarse_motif.target_event_id,
        coarse_motif.emission_position,
    )
    stability = return_spectrum_stability_report(baseline_spectrum, coarse_spectrum)
    subsampling = subsample_reference_chain_positions(
        reference,
        stride=2,
        protected_positions=np.asarray([motif.emission_position]),
    )
    expected_delay = expected_coarse_delay_rank_for_motif(motif, subsampling)
    thinned_network, removed = thin_immediate_trigger_edges(
        network,
        removal_probability=1.0,
        seed=0,
    )
    thinned_closure = transitive_closure_dag(
        immediate_trigger_adjacency(thinned_network)
    )
    after_delay = echo_delay_rank_for_emission(
        thinned_closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    edge_rows = compare_recovery_before_after_edge_thinning(
        closure,
        thinned_closure,
        reference,
        [motif],
    )
    return [
        _row(
            "event_thinning_preserves_spectrum",
            np.array_equal(baseline_spectrum, coarse_spectrum),
            str(coarse_spectrum),
        ),
        _row(
            "spectrum_jaccard_one",
            stability["spectrum_jaccard"] == 1.0,
            stability["spectrum_jaccard"],
        ),
        _row("expected_coarse_delay", expected_delay == 2, float(expected_delay)),
        _row("edge_thinning_removed_edges", removed > 0, float(removed)),
        _row("edge_thinning_missing_recovery", after_delay is None, "missing"),
        _row(
            "edge_compare_became_missing",
            edge_rows[0]["effect"] == "became_missing",
            str(edge_rows[0]["effect"]),
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
    print(f"Wrote echo coarse-graining exact sanity: {output_path}")


if __name__ == "__main__":
    main()
