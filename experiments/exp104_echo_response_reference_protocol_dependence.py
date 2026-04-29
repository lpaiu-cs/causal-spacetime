"""Reference-protocol dependence of echo-response order signatures."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_motifs import insert_multiple_echo_motifs
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_background_state_change_network_with_reference as generate_background,
)
from causal_spacetime_lab.state_change_echo_selection import (
    emission_positions_for_reference_chain,
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
from causal_spacetime_lab.state_change_response_layers import (
    EchoResponseLayerSpec,
    build_layered_echo_motif_specs,
)
from causal_spacetime_lab.state_change_response_signature import (
    EchoResponseSignature,
    echo_response_signature_from_motifs,
    response_order_sign_matrix,
    signature_reachable_fraction,
    signature_tie_fraction,
)
from causal_spacetime_lab.state_change_response_signature_comparison import (
    compare_response_signatures,
    stable_response_order_core,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for reference-protocol dependence diagnostics."""

    num_systems: int = 8
    max_events: int = 300
    trigger_probability: float = 0.20
    emission_position_fraction: float = 0.2
    layer_delay_ranks: tuple[int, ...] = (3, 5, 8, 13)
    targets_per_layer: int = 10
    repetitions: int = 5
    seed: int = 0
    top_k: int = 5
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str]) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Echo response reference-protocol dependence."
    )
    parser.add_argument("--num-systems", type=int, default=8)
    parser.add_argument("--max-events", type=int, default=300)
    parser.add_argument("--trigger-probability", type=float, default=0.20)
    parser.add_argument("--emission-position-fraction", type=float, default=0.2)
    parser.add_argument(
        "--layer-delay-ranks",
        nargs="+",
        default=["3", "5", "8", "13"],
    )
    parser.add_argument("--targets-per-layer", type=int, default=10)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        num_systems=args.num_systems,
        max_events=args.max_events,
        trigger_probability=args.trigger_probability,
        emission_position_fraction=args.emission_position_fraction,
        layer_delay_ranks=_parse_int_values(args.layer_delay_ranks),
        targets_per_layer=args.targets_per_layer,
        repetitions=args.repetitions,
        seed=args.seed,
        top_k=args.top_k,
        output_dir=args.output_dir,
    )


def _primary_emission_position(reference: np.ndarray, config: ExperimentConfig) -> int:
    max_delay = max(config.layer_delay_ranks)
    upper = reference.size - max_delay - 1
    if upper < 1:
        raise ValueError("primary reference chain is too short for layer delays")
    proposed = int(round(config.emission_position_fraction * (reference.size - 1)))
    return int(min(max(1, proposed), upper))


def _layer_specs(config: ExperimentConfig) -> list[EchoResponseLayerSpec]:
    return [
        EchoResponseLayerSpec(delay, config.targets_per_layer)
        for delay in config.layer_delay_ranks
    ]


def _candidate_map(
    network,
    closure: np.ndarray,
    seed: int,
) -> dict[str, ReferenceChainCandidate]:
    candidates = local_system_reference_chain_candidates(network, min_length=3)
    candidates.extend(
        [
            greedy_reference_chain_candidate_from_order(closure),
            longest_reference_chain_candidate_from_order(closure),
            random_reference_chain_candidate_from_order(closure, seed=seed),
        ]
    )
    return {candidate.name: candidate for candidate in candidates}


def _signature_for_reference(
    closure: np.ndarray,
    reference: np.ndarray,
    target_ids: np.ndarray,
    emission_position: int,
    label: str,
) -> EchoResponseSignature:
    delays: list[int] = []
    reachable: list[bool] = []
    for target_id in target_ids:
        recovered = echo_delay_rank_for_emission(
            closure,
            reference,
            int(target_id),
            emission_position,
        )
        delays.append(int(recovered) if recovered is not None else -1)
        reachable.append(recovered is not None)
    delay_array = np.asarray(delays, dtype=int)
    reachable_array = np.asarray(reachable, dtype=bool)
    return EchoResponseSignature(
        target_event_ids=np.asarray(target_ids, dtype=int),
        delay_ranks=delay_array,
        reachable_mask=reachable_array,
        order_sign_matrix=response_order_sign_matrix(delay_array, reachable_array),
        label=label,
    )


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run reference-protocol dependence diagnostics."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        seed = config.seed + repetition
        network, primary_reference = generate_background(
            config.num_systems,
            config.max_events,
            config.trigger_probability,
            reference_source="highest_utility",
            seed=seed,
        )
        if primary_reference.size <= max(config.layer_delay_ranks) + 2:
            continue
        emission = _primary_emission_position(primary_reference, config)
        specs = build_layered_echo_motif_specs(
            primary_reference,
            emission,
            _layer_specs(config),
            seed=seed + 1000,
        )
        network, motifs = insert_multiple_echo_motifs(
            network,
            primary_reference,
            specs,
        )
        closure = transitive_closure_dag(immediate_trigger_adjacency(network))
        primary_signature = echo_response_signature_from_motifs(
            closure,
            primary_reference,
            motifs,
            label="primary",
        )
        target_ids = primary_signature.target_event_ids
        candidates = _candidate_map(network, closure, seed + 2000)
        reports = [
            evaluate_reference_chain_candidate(network, closure, candidate)
            for candidate in candidates.values()
        ]
        ranked = rank_reference_chain_candidates(reports)
        selected_names = [str(row["name"]) for row in ranked[: config.top_k]]
        for rank, name in enumerate(selected_names, start=1):
            candidate = candidates[name]
            emissions = emission_positions_for_reference_chain(
                candidate.chain_event_ids,
                strategy="interior_quantiles",
                count=1,
            )
            if emissions.size == 0:
                continue
            alt_signature = _signature_for_reference(
                closure,
                candidate.chain_event_ids,
                target_ids,
                int(emissions[0]),
                label=name,
            )
            comparison = compare_response_signatures(primary_signature, alt_signature)
            core = stable_response_order_core([primary_signature, alt_signature])
            rows.append(
                {
                    "repetition": float(repetition),
                    "reference_rank": float(rank),
                    "reference_name": name,
                    "reference_source": candidate.source,
                    "reference_length": float(candidate.chain_event_ids.size),
                    "reachable_fraction": signature_reachable_fraction(alt_signature),
                    "pair_agreement_fraction": comparison[
                        "pair_agreement_fraction"
                    ],
                    "pair_disagreement_fraction": comparison[
                        "pair_disagreement_fraction"
                    ],
                    "tie_fraction": signature_tie_fraction(alt_signature),
                    "stable_pair_fraction": float(core["stable_pair_fraction"]),
                }
            )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write experiment CSV."""

    path = output_dir / "data" / "echo_response_reference_protocol_dependence.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _group_means(
    rows: list[dict[str, float | str]],
    value: str,
) -> tuple[list[str], list[float]]:
    groups = sorted({str(row["reference_source"]) for row in rows})
    means: list[float] = []
    for group in groups:
        values = np.asarray(
            [
                float(row[value])
                for row in rows
                if row["reference_source"] == group
            ],
            dtype=float,
        )
        values = values[np.isfinite(values)]
        means.append(float(np.mean(values)) if values.size else float("nan"))
    return groups, means


def save_figures(rows: list[dict[str, float | str]], output_dir: Path) -> list[Path]:
    """Save summary figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for value, filename, ylabel in [
        (
            "pair_agreement_fraction",
            "echo_response_reference_pair_agreement.png",
            "pair agreement with primary",
        ),
        (
            "reachable_fraction",
            "echo_response_reference_reachability.png",
            "reachable target fraction",
        ),
    ]:
        labels, means = _group_means(rows, value)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(labels, means)
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 1.02)
        ax.tick_params(axis="x", rotation=25)
        ax.grid(True, axis="y", alpha=0.3)
        path = figure_dir / filename
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote echo response reference-protocol dependence: {data_path}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
