"""Finite state-change causal trigger network diagnostics."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    is_irreflexive,
    is_transitive,
    local_finiteness_report,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_validation import (
    branching_statistics,
    state_change_network_summary,
    trigger_graph_summary,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for finite state-change toy-model diagnostics."""

    num_systems_values: tuple[int, ...] = (5, 10)
    max_events_values: tuple[int, ...] = (100, 300, 600)
    trigger_probability_values: tuple[float, ...] = (0.15, 0.30)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="State-change toy model.")
    parser.add_argument("--num-systems-values", nargs="+", default=["5", "10"])
    parser.add_argument("--max-events-values", nargs="+", default=["100", "300", "600"])
    parser.add_argument(
        "--trigger-probability-values",
        nargs="+",
        default=["0.15", "0.30"],
    )
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems_values=_parse_int_values(args.num_systems_values),
        max_events_values=_parse_int_values(args.max_events_values),
        trigger_probability_values=_parse_float_values(
            args.trigger_probability_values
        ),
        max_triggers_per_event=args.max_triggers_per_event,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float]]:
    """Run finite state-change toy-model diagnostics."""

    rows: list[dict[str, float]] = []
    for num_systems in config.num_systems_values:
        for max_events in config.max_events_values:
            for trigger_probability in config.trigger_probability_values:
                for repetition in range(config.repetitions):
                    network = generate_state_change_network(
                        num_systems,
                        max_events,
                        trigger_probability=trigger_probability,
                        max_triggers_per_event=config.max_triggers_per_event,
                        seed=(
                            config.seed
                            + 100000 * num_systems
                            + 100 * max_events
                            + repetition
                        ),
                    )
                    adjacency = immediate_trigger_adjacency(network)
                    closure = transitive_closure_dag(adjacency)
                    row = {
                        "num_systems": float(num_systems),
                        "max_events": float(max_events),
                        "trigger_probability": float(trigger_probability),
                        "max_triggers_per_event": float(
                            config.max_triggers_per_event
                        ),
                        "repetition": float(repetition),
                        "is_irreflexive": float(is_irreflexive(closure)),
                        "is_transitive": float(is_transitive(closure)),
                    }
                    row.update(state_change_network_summary(network))
                    row.update(trigger_graph_summary(adjacency, closure))
                    row.update(local_finiteness_report(closure))
                    row.update(branching_statistics(network))
                    rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, float]], output_dir: Path) -> Path:
    """Write state-change toy-model summary rows."""

    output_path = output_dir / "data" / "state_change_toy_model_summary.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_max_events(
    rows: list[dict[str, float]],
    key: str,
    trigger_probability: float,
) -> tuple[list[float], list[float]]:
    max_events_values = sorted(
        {
            row["max_events"]
            for row in rows
            if row["trigger_probability"] == trigger_probability
        }
    )
    means = [
        float(
            np.mean(
                [
                    row[key]
                    for row in rows
                    if row["trigger_probability"] == trigger_probability
                    and row["max_events"] == max_events
                ]
            )
        )
        for max_events in max_events_values
    ]
    return max_events_values, means


def save_figures(
    rows: list[dict[str, float]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save state-change toy-model figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "state_change_relation_density.png",
            "relation_density",
            "Relation density",
        ),
        (
            "state_change_max_interval_size.png",
            "max_interval_size",
            "Max interval size",
        ),
        (
            "state_change_events_per_system.png",
            "mean_events_per_system",
            "Mean events per system",
        ),
    )
    paths: list[Path] = []
    trigger_values = sorted({row["trigger_probability"] for row in rows})
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        for trigger_probability in trigger_values:
            xs, ys = _mean_by_max_events(rows, key, trigger_probability)
            ax.plot(xs, ys, marker="o", label=f"p={trigger_probability:.2f}")
        ax.set_xlabel("Maximum stored events")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)
    return tuple(paths)  # type: ignore[return-value]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote state-change toy-model summary: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
