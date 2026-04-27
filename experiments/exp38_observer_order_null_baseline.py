"""Observer-derived distance order compared with null baselines."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    append_oriented_protocol_chains,
    common_safe_tau_range_for_oriented_protocol_1p1,
    reconstruct_oriented_chart_from_order,
)
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    procrustes_align,
    sample_quadruplet_constraints_from_distance_values,
    squared_distance_matrix,
)
from causal_spacetime_lab.representation_validation import (
    random_quadruplet_constraints,
    representation_generalization_report,
    shuffle_constraint_sides,
    split_constraints,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for observer-order null baseline diagnostics."""

    T: float = 2.0
    n_events: int = 400
    tick_values: tuple[int, ...] = (32, 64, 128)
    constraint_counts: tuple[int, ...] = (1000, 3000)
    repetitions: int = 5
    beacon_separation: float = 0.15
    seed: int = 0
    steps: int = 800
    restarts: int = 2
    learning_rate: float = 0.1
    output_dir: Path = DEFAULT_OUTPUT_DIR


def _parse_int_values(values: list[str], name: str) -> tuple[int, ...]:
    parsed: list[int] = []
    for value in values:
        parsed.extend(int(part) for part in value.split(",") if part)
    if not parsed:
        raise argparse.ArgumentTypeError(f"at least one {name} value is required")
    return tuple(parsed)


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Compare observer-derived order with null baselines."
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--N", type=int, default=400)
    parser.add_argument("--tick-values", nargs="+", default=["32", "64", "128"])
    parser.add_argument("--constraint-counts", nargs="+", default=["1000", "3000"])
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=800)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=0.1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_events=args.N,
        tick_values=_parse_int_values(args.tick_values, "tick"),
        constraint_counts=_parse_int_values(args.constraint_counts, "constraint"),
        repetitions=args.repetitions,
        beacon_separation=args.beacon_separation,
        seed=args.seed,
        steps=args.steps,
        restarts=args.restarts,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )


def _reconstruct_positions(
    support_events: np.ndarray,
    tick_count: int,
    config: ExperimentConfig,
) -> tuple[np.ndarray, np.ndarray]:
    spec = ObserverProtocolSpec(
        name="stationary_oriented",
        beta=0.0,
        origin_lab_time=0.0,
        origin_lab_position=0.0,
        beacon_separation=config.beacon_separation,
    )
    tau_range = common_safe_tau_range_for_oriented_protocol_1p1(config.T, spec)
    chain_events: list[np.ndarray] = []
    indices, _ = append_oriented_protocol_chains(
        chain_events,
        spec,
        tick_count,
        tau_range,
        support_events.shape[0],
    )
    combined = np.vstack((support_events, *chain_events))
    chart = reconstruct_oriented_chart_from_order(
        causal_matrix_1p1(combined),
        np.arange(support_events.shape[0], dtype=int),
        indices.primary,
        indices.beacon,
        indices.clock_times,
        config.beacon_separation,
        spec.name,
    )
    accessible = chart.accessible
    return support_events[accessible, 1], chart.reconstructed_coords[accessible, 1]


def _constraints_for_model(
    observer_constraints: np.ndarray,
    accessible_count: int,
    model: str,
    seed: int,
) -> np.ndarray:
    if model == "observer":
        return observer_constraints
    if model == "shuffled":
        return shuffle_constraint_sides(observer_constraints, seed=seed)
    if model == "random":
        return random_quadruplet_constraints(
            accessible_count,
            observer_constraints.shape[0],
            seed=seed,
        )
    raise ValueError(f"unknown model: {model}")


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run observer-derived order/null baseline comparison."""

    rows: list[dict[str, float | str]] = []
    for repetition in range(config.repetitions):
        support_events = sprinkle_1p1_causal_diamond(
            config.n_events,
            T=config.T,
            seed=config.seed + repetition,
        )
        for tick_count in config.tick_values:
            true_x, reconstructed_x = _reconstruct_positions(
                support_events,
                tick_count,
                config,
            )
            accessible_count = true_x.size
            if accessible_count < 4:
                continue
            distance_matrix = squared_distance_matrix(reconstructed_x[:, None])
            for constraint_count in config.constraint_counts:
                observer_constraints = (
                    sample_quadruplet_constraints_from_distance_values(
                        distance_matrix,
                        constraint_count,
                        seed=config.seed + 10_000 * tick_count + repetition,
                    )
                )
                for model_index, model in enumerate(("observer", "shuffled", "random")):
                    constraints = _constraints_for_model(
                        observer_constraints,
                        accessible_count,
                        model,
                        seed=config.seed + 100_000 * model_index + repetition,
                    )
                    train, test = split_constraints(
                        constraints,
                        train_fraction=0.7,
                        seed=config.seed + 1_000_000 * model_index + repetition,
                    )
                    embedding, diagnostics = fit_ordinal_embedding_gradient_descent(
                        accessible_count,
                        1,
                        train,
                        steps=config.steps,
                        restarts=config.restarts,
                        learning_rate=config.learning_rate,
                        seed=config.seed + 20_000 * tick_count + model_index,
                        batch_size=min(2048, max(1, train.shape[0])),
                    )
                    report = representation_generalization_report(
                        embedding,
                        train,
                        test,
                    )
                    aligned, alignment = procrustes_align(embedding, true_x[:, None])
                    rows.append(
                        {
                            "model": model,
                            "tick_count": float(tick_count),
                            "constraint_count": float(constraint_count),
                            "repetition": float(repetition),
                            "accessible_count": float(accessible_count),
                            "accessible_fraction": float(
                                accessible_count / config.n_events
                            ),
                            "final_train_loss": diagnostics["final_loss"],
                            **report,
                            "distance_order_error_against_true_x": (
                                embedding_distance_order_error(
                                    aligned,
                                    true_x[:, None],
                                    seed=config.seed + repetition,
                                )
                            ),
                            "alignment_rmse": alignment["rmse"],
                        }
                    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write observer-order baseline rows."""

    output_path = output_dir / "data" / "observer_order_null_baseline.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def _plot(
    rows: list[dict[str, float | str]],
    key: str,
    ylabel: str,
    output_path: Path,
) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    for model in sorted({str(row["model"]) for row in rows}):
        subset = [row for row in rows if row["model"] == model]
        x_values = sorted({float(row["tick_count"]) for row in subset})
        y_values = []
        for x in x_values:
            values = [
                float(row[key])
                for row in subset
                if float(row["tick_count"]) == x
            ]
            y_values.append(float(np.nanmean(values)))
        ax.plot(x_values, y_values, marker="o", label=model)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel + " vs tick count")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Save observer-order baseline figures."""

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return (
        _plot(
            rows,
            "test_violation_rate",
            "Held-out violation rate",
            figure_dir / "observer_order_vs_null_test_violation.png",
        ),
        _plot(
            rows,
            "alignment_rmse",
            "Alignment RMSE",
            figure_dir / "observer_order_vs_null_alignment_rmse.png",
        ),
        _plot(
            rows,
            "distance_order_error_against_true_x",
            "Distance-order error against true x",
            figure_dir / "observer_order_vs_null_distance_order_error.png",
        ),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    data_path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote observer-order null baseline data: {data_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
