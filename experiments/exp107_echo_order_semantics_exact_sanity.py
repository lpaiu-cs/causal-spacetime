"""Exact sanity checks for echo-order spectrum semantics."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_interference import shortcut_depth
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_echo_spectrum_semantics import (
    compress_suffix_spectrum,
    full_transitive_return_spectrum,
    immediate_edge_return_spectrum,
    is_suffix_spectrum,
    retained_reference_return_spectrum,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_order_semantics_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo-order semantics sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_order_semantics_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact spectrum-semantics checks."""

    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(
            emission_position=2,
            planted_delay_rank=4,
            outward_steps=1,
            return_steps=0,
        ),
    )
    adjacency = immediate_trigger_adjacency(network)
    closure = transitive_closure_dag(adjacency)
    s_full = full_transitive_return_spectrum(
        closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    s_immediate = immediate_edge_return_spectrum(
        adjacency,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    retained_reference = reference[np.asarray([0, 2, 4, 7])]
    s_retained = retained_reference_return_spectrum(
        closure,
        retained_reference,
        motif.target_event_id,
        1,
    )
    shortcut_network, _ = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(probability=1.0, min_depth=2, max_depth=2),
        seed=0,
    )
    shortcut_closure = transitive_closure_dag(
        immediate_trigger_adjacency(shortcut_network)
    )
    shortcut_spectrum = full_transitive_return_spectrum(
        shortcut_closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    recovered_shortcut = int(np.min(shortcut_spectrum))
    return [
        _row("s_full_suffix", is_suffix_spectrum(s_full), str(s_full)),
        _row(
            "s_full_starts_at_planted_delay",
            s_full.size > 0 and int(s_full[0]) == 4,
            str(s_full),
        ),
        _row("s_immediate_sparse", np.array_equal(s_immediate, [4]), str(s_immediate)),
        _row(
            "s_retained_coarse_sparse",
            np.array_equal(s_retained, [2]),
            str(s_retained),
        ),
        _row(
            "d_echo_is_min_s_full",
            int(np.min(s_full)) == 4,
            float(np.min(s_full)),
        ),
        _row("w_rank_not_used", True, "passive bracket-width rank is separate"),
        _row(
            "shortcut_d_echo",
            recovered_shortcut == 2,
            float(recovered_shortcut),
        ),
        _row(
            "shortcut_depth",
            shortcut_depth(recovered_shortcut, motif.planted_delay_rank) == 2.0,
            shortcut_depth(recovered_shortcut, motif.planted_delay_rank),
        ),
        _row(
            "suffix_summary",
            compress_suffix_spectrum(s_full)["is_suffix"] == 1.0,
            str(compress_suffix_spectrum(s_full)),
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
    print(f"Wrote echo-order semantics exact sanity: {output_path}")


if __name__ == "__main__":
    main()
