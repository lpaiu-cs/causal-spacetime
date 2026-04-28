"""Exact sanity checks for echo shortcut return classification."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from causal_spacetime_lab.state_change_echo_interference import (
    return_delay_spectrum,
    return_spectrum_report_for_motif,
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
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT = Path("outputs/data/echo_shortcut_exact_sanity.csv")


def parse_args() -> Path:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo shortcut exact sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return args.output_dir / "data" / "echo_shortcut_exact_sanity.csv"


def _row(name: str, passed: bool, value: float | str) -> dict[str, float | str]:
    return {"check": name, "passed": float(passed), "value": value}


def run_experiment() -> list[dict[str, float | str]]:
    """Run exact shortcut classification checks."""

    network, reference = generate_reference_backbone_network(10)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(emission_position=2, planted_delay_rank=4),
    )
    clean_closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    clean_spectrum = return_delay_spectrum(
        clean_closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )
    network, shortcut_records = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(
            probability=1.0,
            min_depth=2,
            max_depth=2,
        ),
        seed=0,
    )
    shortcut_closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    report = return_spectrum_report_for_motif(shortcut_closure, reference, motif)
    return [
        _row("clean_spectrum_has_planted", 4 in clean_spectrum, str(clean_spectrum)),
        _row(
            "clean_spectrum_has_later",
            np.max(clean_spectrum) > 4,
            str(clean_spectrum),
        ),
        _row(
            "shortcut_recorded",
            len(shortcut_records) == 1,
            float(len(shortcut_records)),
        ),
        _row("recovered_delay_shortcut", report["recovered_delay_rank"] == 2.0, 2.0),
        _row(
            "shortcut_count",
            float(report["shortcut_count"]) >= 1.0,
            report["shortcut_count"],
        ),
        _row(
            "shortcut_depth",
            report["shortcut_depth"] == 2.0,
            report["shortcut_depth"],
        ),
        _row(
            "early_shortcut",
            report["early_shortcut"] == 1.0,
            report["early_shortcut"],
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
    print(f"Wrote echo shortcut exact sanity: {output_path}")


if __name__ == "__main__":
    main()
