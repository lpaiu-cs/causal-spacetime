"""Reference-protocol dependence of echo shortcut classifications."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.ordinal import order_inversion_rate
from causal_spacetime_lab.state_change_echo_interference import (
    return_spectrum_report_for_motif,
    summarize_return_spectrum_reports,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_multiple_echo_motifs,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_background_state_change_network_with_reference,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_reference_chains import (
    ReferenceChainCandidate,
    greedy_reference_chain_candidate_from_order,
    local_system_reference_chain_candidates,
    longest_reference_chain_candidate_from_order,
    random_reference_chain_candidate_from_order,
)
from causal_spacetime_lab.state_change_reference_quality import (
    evaluate_reference_chain_candidate,
    rank_reference_chain_candidates,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for reference-dependent shortcut diagnostics."""

    num_systems: int = 8
    max_events: int = 300
    trigger_probability: float = 0.20
    motif_count: int = 30
    shortcut_probability: float = 0.3
    repetitions: int = 5
    seed: int = 0
    top_k: int = 5
    delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Echo shortcut reference dependence.")
    parser.add_argument("--num-systems", type=int, default=8)
    parser.add_argument("--max-events", type=int, default=300)
    parser.add_argument("--trigger-probability", type=float, default=0.20)
    parser.add_argument("--motif-count", type=int, default=30)
    parser.add_argument("--shortcut-probability", type=float, default=0.3)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--delay-ranks", nargs="+", default=["3", "5", "8", "13"])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability=args.trigger_probability,
        motif_count=args.motif_count,
        shortcut_probability=args.shortcut_probability,
        repetitions=args.repetitions,
        seed=args.seed,
        top_k=args.top_k,
        delay_ranks=_parse_int_values(args.delay_ranks),
        output_dir=args.output_dir,
    )


def _random_specs(
    reference_length: int,
    config: ExperimentConfig,
    seed: int,
) -> list[EchoMotifSpec]:
    rng = np.random.default_rng(seed)
    valid_delays = np.asarray(
        [delay for delay in config.delay_ranks if delay < reference_length],
        dtype=int,
    )
    specs: list[EchoMotifSpec] = []
    for _ in range(config.motif_count):
        delay = int(rng.choice(valid_delays))
        emission = int(rng.integers(0, reference_length - delay))
        specs.append(EchoMotifSpec(emission, delay))
    return specs


def _candidate_references(
    closure: np.ndarray,
    network,
    primary_reference: np.ndarray,
    top_k: int,
    seed: int,
) -> list[tuple[ReferenceChainCandidate, dict[str, float | str]]]:
    candidates = [
        ReferenceChainCandidate("primary_reference", primary_reference, "manual"),
        *local_system_reference_chain_candidates(network, min_length=2),
        greedy_reference_chain_candidate_from_order(closure),
        longest_reference_chain_candidate_from_order(closure),
    ]
    for index in range(3):
        candidates.append(
            random_reference_chain_candidate_from_order(
                closure,
                seed=seed + index,
                name=f"random_reference_chain_{index}",
            )
        )
    reports = [
        evaluate_reference_chain_candidate(network, closure, candidate)
        for candidate in candidates
    ]
    ranked = rank_reference_chain_candidates(reports)
    by_name = {candidate.name: candidate for candidate in candidates}
    selected: list[tuple[ReferenceChainCandidate, dict[str, float | str]]] = []
    seen: set[str] = set()
    for row in ranked:
        name = str(row["name"])
        if name in seen:
            continue
        seen.add(name)
        selected.append((by_name[name], row))
        if len(selected) >= top_k:
            break
    if "primary_reference" not in seen:
        primary_row = next(row for row in ranked if row["name"] == "primary_reference")
        selected.insert(0, (by_name["primary_reference"], primary_row))
    return selected


def _reports_for_reference(
    closure: np.ndarray,
    reference: np.ndarray,
    motifs,
) -> list[dict[str, float | str]]:
    reports: list[dict[str, float | str]] = []
    if reference.size < 2:
        return reports
    for motif in motifs:
        emission = min(motif.emission_position, reference.size - 2)
        reports.append(
            return_spectrum_report_for_motif(
                closure,
                reference,
                replace(motif, emission_position=emission),
            )
        )
    return reports


def _recovered_array(rows: list[dict[str, float | str]]) -> np.ndarray:
    return np.asarray([float(row["recovered_delay_rank"]) for row in rows], dtype=float)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run reference-dependent shortcut classification diagnostics."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        seed = config.seed + repetition
        network, primary_reference = (
            generate_background_state_change_network_with_reference(
                config.num_systems,
                config.max_events,
                config.trigger_probability,
                reference_source="highest_utility",
                seed=seed,
            )
        )
        if primary_reference.size < 3:
            continue
        specs = _random_specs(primary_reference.size, config, seed + 1000)
        network, motifs = insert_multiple_echo_motifs(
            network,
            primary_reference,
            specs,
        )
        network, shortcut_records = inject_shortcut_returns(
            network,
            primary_reference,
            motifs,
            ShortcutInjectionSpec(probability=config.shortcut_probability),
            seed=seed + 2000,
        )
        closure = transitive_closure_dag(immediate_trigger_adjacency(network))
        primary_reports = _reports_for_reference(closure, primary_reference, motifs)
        primary_summary = summarize_return_spectrum_reports(primary_reports)
        primary_recovered = _recovered_array(primary_reports)
        primary_visible = np.isfinite(primary_recovered)
        for candidate, rank_row in _candidate_references(
            closure,
            network,
            primary_reference,
            config.top_k,
            seed + 3000,
        ):
            reports = _reports_for_reference(closure, candidate.chain_event_ids, motifs)
            summary = summarize_return_spectrum_reports(reports)
            recovered = _recovered_array(reports)
            visible = np.isfinite(recovered)
            common = visible & primary_visible
            inversion = (
                order_inversion_rate(
                    primary_recovered[common],
                    recovered[common],
                    ignore_ties=True,
                )
                if np.count_nonzero(common) >= 2
                else float("nan")
            )
            rows.append(
                {
                    "num_systems": float(config.num_systems),
                    "max_events": float(config.max_events),
                    "trigger_probability": float(config.trigger_probability),
                    "repetition": float(repetition),
                    "candidate_name": candidate.name,
                    "candidate_source": candidate.source,
                    "utility_score": float(rank_row["score"]),
                    "utility_rank": float(rank_row["rank"]),
                    "chain_length": float(candidate.chain_event_ids.size),
                    "motif_count": float(len(motifs)),
                    "shortcut_injection_count": float(len(shortcut_records)),
                    "primary_exact_recovery_fraction": primary_summary[
                        "exact_recovery_fraction"
                    ],
                    "primary_shortcut_fraction": primary_summary[
                        "shortcut_fraction"
                    ],
                    "alternative_visible_fraction": float(np.mean(visible))
                    if visible.size
                    else 0.0,
                    "alternative_shortcut_fraction": summary["shortcut_fraction"],
                    "alternative_exact_recovery_fraction": summary[
                        "exact_recovery_fraction"
                    ],
                    "common_visible_count": float(np.count_nonzero(common)),
                    "order_inversion_against_primary": inversion,
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write reference-dependent shortcut rows."""

    output_path = output_dir / "data" / "echo_shortcut_reference_dependence.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _mean_by_source(
    rows: list[dict[str, float | str]],
    key: str,
) -> tuple[list[str], list[float]]:
    sources = sorted({str(row["candidate_source"]) for row in rows})
    means: list[float] = []
    for source in sources:
        values = [
            float(row[key])
            for row in rows
            if row["candidate_source"] == source and np.isfinite(float(row[key]))
        ]
        means.append(float(np.mean(values)) if values else float("nan"))
    return sources, means


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Save reference-dependent shortcut figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    visibility_path = figure_dir / "echo_shortcut_visibility_by_reference.png"
    disagreement_path = figure_dir / "echo_shortcut_order_disagreement_by_reference.png"
    for path, key, ylabel in (
        (visibility_path, "alternative_visible_fraction", "Visible motif fraction"),
        (
            disagreement_path,
            "order_inversion_against_primary",
            "Order inversion against primary",
        ),
    ):
        sources, means = _mean_by_source(rows, key)
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.bar(sources, means)
        ax.set_xlabel("Reference-chain source")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, axis="y", alpha=0.3)
        fig.autofmt_xdate(rotation=20)
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return visibility_path, disagreement_path


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo shortcut reference dependence: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
