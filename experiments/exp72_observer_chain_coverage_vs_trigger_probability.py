"""Observer-like chain coverage versus external trigger probability."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_observer_quality import (
    evaluate_chain_candidate,
    rank_chain_candidates,
)
from causal_spacetime_lab.state_change_observers import (
    greedy_chain_candidate_from_order,
    local_system_chain_candidates,
    longest_chain_candidate_from_order,
    random_chain_candidate_from_order,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for chain coverage versus trigger probability."""

    num_systems: int = 10
    max_events: int = 500
    trigger_probability_values: tuple[float, ...] = (0.05, 0.10, 0.20, 0.35, 0.50)
    max_triggers_per_event: int = 2
    repetitions: int = 5
    seed: int = 0
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_float_values(values: list[str]) -> tuple[float, ...]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Observer-chain coverage versus trigger probability."
    )
    parser.add_argument("--num-systems", type=int, default=10)
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument(
        "--trigger-probability-values",
        nargs="+",
        default=["0.05", "0.10", "0.20", "0.35", "0.50"],
    )
    parser.add_argument("--max-triggers-per-event", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability_values=_parse_float_values(
            args.trigger_probability_values
        ),
        max_triggers_per_event=args.max_triggers_per_event,
        repetitions=args.repetitions,
        seed=args.seed,
        output_dir=args.output_dir,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run coverage diagnostics over trigger probabilities."""

    rows: list[dict[str, float | str]] = []
    for trigger_probability in config.trigger_probability_values:
        for repetition in range(config.repetitions):
            run_seed = config.seed + int(10000 * trigger_probability) + repetition
            network = generate_state_change_network(
                config.num_systems,
                config.max_events,
                trigger_probability=trigger_probability,
                max_triggers_per_event=config.max_triggers_per_event,
                seed=run_seed,
            )
            closure = transitive_closure_dag(immediate_trigger_adjacency(network))
            local_candidates = local_system_chain_candidates(network)
            candidates = [
                *local_candidates,
                greedy_chain_candidate_from_order(closure),
                longest_chain_candidate_from_order(closure),
                random_chain_candidate_from_order(closure, seed=run_seed + 17),
            ]
            reports = [
                evaluate_chain_candidate(network, closure, candidate)
                for candidate in candidates
            ]
            ranked = rank_chain_candidates(reports)
            local_reports = [
                report for report in reports if report.source == "local_system"
            ]
            row = {
                "num_systems": float(config.num_systems),
                "max_events": float(config.max_events),
                "trigger_probability": float(trigger_probability),
                "max_triggers_per_event": float(config.max_triggers_per_event),
                "repetition": float(repetition),
                "best_score": float(ranked[0]["score"]),
                "best_source": str(ranked[0]["source"]),
                "best_comparable_fraction": float(
                    ranked[0]["comparable_fraction"]
                ),
                "best_bracketed_fraction": float(ranked[0]["bracketed_fraction"]),
                "mean_local_comparable_fraction": float(
                    np.mean(
                        [report.comparable_fraction for report in local_reports]
                    )
                )
                if local_reports
                else float("nan"),
                "mean_local_bracketed_fraction": float(
                    np.mean([report.bracketed_fraction for report in local_reports])
                )
                if local_reports
                else float("nan"),
            }
            rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write coverage diagnostics."""

    output_path = (
        output_dir / "data" / "observer_chain_coverage_vs_trigger_probability.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_probability(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[float], list[float]]:
    probabilities = sorted({float(row["trigger_probability"]) for row in rows})
    means = [
        float(
            np.mean(
                [
                    float(row[key])
                    for row in rows
                    if float(row["trigger_probability"]) == probability
                ]
            )
        )
        for probability in probabilities
    ]
    return probabilities, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save coverage and bracketing figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specs = (
        (
            "observer_chain_coverage_vs_trigger_probability.png",
            "best_comparable_fraction",
            "Best comparable fraction",
        ),
        (
            "observer_chain_bracketing_vs_trigger_probability.png",
            "best_bracketed_fraction",
            "Best bracketed fraction",
        ),
    )
    paths: list[Path] = []
    for filename, key, ylabel in specs:
        path = figure_dir / filename
        xs, ys = _mean_by_probability(rows, key)
        fig, ax = plt.subplots(figsize=(7.0, 4.7))
        ax.plot(xs, ys, marker="o")
        ax.set_xlabel("External trigger probability")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, alpha=0.3)
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
    print(f"Wrote observer-chain coverage data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
