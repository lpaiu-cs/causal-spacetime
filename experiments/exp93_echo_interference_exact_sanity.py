"""Exact sanity checks for echo interference utilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_interference import (
    return_delay_spectrum,
    return_positions_for_target,
    return_spectrum_report_for_motif,
    shortcut_depth,
    shortcut_positions_before_planted,
    summarize_return_spectrum_reports,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    add_random_acyclic_background_edges,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_interference_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo interference exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_interference_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def _is_acyclic(network) -> bool:
    try:
        topological_order_from_adjacency(immediate_trigger_adjacency(network))
    except ValueError:
        return False
    return True


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact shortcut-interference checks."""

    network, reference = generate_reference_backbone_network(10)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(emission_position=2, planted_delay_rank=4),
    )
    clean_closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    positions = return_positions_for_target(
        clean_closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    spectrum = return_delay_spectrum(
        clean_closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    shortcut_network, shortcut_records = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(probability=1.0, min_depth=2, max_depth=2),
        seed=0,
    )
    shortcut_closure = transitive_closure_dag(
        immediate_trigger_adjacency(shortcut_network)
    )
    shortcut_report = return_spectrum_report_for_motif(
        shortcut_closure,
        reference,
        motif,
    )
    perturbed, added_edges = add_random_acyclic_background_edges(
        network,
        edge_probability=0.2,
        seed=1,
        max_edges=5,
    )
    summary = summarize_return_spectrum_reports([shortcut_report])
    shortcuts = shortcut_positions_before_planted(np.asarray([2, 4, 5]), 4)
    return [
        _row("return_positions_include_planted", 6 in positions, str(positions)),
        _row("delay_spectrum_includes_planted", 4 in spectrum, str(spectrum)),
        _row("delay_spectrum_has_later", np.max(spectrum) > 4, str(spectrum)),
        _row(
            "shortcut_positions_before_planted",
            np.array_equal(shortcuts, [2]),
            str(shortcuts),
        ),
        _row("shortcut_depth_exact", shortcut_depth(4, 4) == 0.0, 0.0),
        _row("shortcut_depth_early", shortcut_depth(2, 4) == 2.0, 2.0),
        _row("shortcut_depth_late", shortcut_depth(5, 4) == -1.0, -1.0),
        _row("shortcut_depth_missing", np.isnan(shortcut_depth(None, 4)), "nan"),
        _row(
            "return_spectrum_report_fields",
            "shortcut_count" in shortcut_report
            and "return_spectrum_string" in shortcut_report,
            str(sorted(shortcut_report)),
        ),
        _row("shortcut_injection_recorded", len(shortcut_records) == 1, 1.0),
        _row("shortcut_injection_acyclic", _is_acyclic(shortcut_network), 1.0),
        _row("background_edges_added", added_edges > 0, float(added_edges)),
        _row("background_edges_acyclic", _is_acyclic(perturbed), 1.0),
        _row("summary_fields", summary["motif_count"] == 1.0, str(sorted(summary))),
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
    print(f"Wrote echo interference exact sanity: {output_path}")


if __name__ == "__main__":
    main()
