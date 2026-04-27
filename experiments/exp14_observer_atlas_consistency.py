"""Observer-atlas consistency for affine oriented radar protocols."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from causal_spacetime_lab.atlas_validation import chart_interval_disagreement
from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    ProtocolChainIndices,
    ReconstructedChart,
    affine_lab_to_rest_1p1,
    append_oriented_protocol_chains,
    common_safe_tau_range_for_oriented_protocol_1p1,
    reconstruct_oriented_chart_from_order,
)
from causal_spacetime_lab.poincare_maps import (
    compose_affine_lorentz_maps_1p1,
    expected_relative_beta_1p1,
    fit_affine_lorentz_beta_grid_1p1,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

DEFAULT_N_VALUES = (300, 600, 1200)
DEFAULT_TICK_VALUES = (32, 64, 128)
DEFAULT_OUTPUT_DIR = Path("outputs")
DEFAULT_TRANSITIONS = (
    ("A_lab", "B_moving_pos"),
    ("A_lab", "C_moving_neg"),
    ("B_moving_pos", "C_moving_neg"),
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for observer-atlas consistency validation."""

    T: float = 2.0
    n_values: tuple[int, ...] = DEFAULT_N_VALUES
    tick_values: tuple[int, ...] = DEFAULT_TICK_VALUES
    repetitions: int = 5
    seed: int = 0
    beacon_separation: float = 0.15
    output_dir: Path = DEFAULT_OUTPUT_DIR
    num_invariant_pairs: int = 500


def default_protocols(beacon_separation: float) -> tuple[ObserverProtocolSpec, ...]:
    """Return the default three-protocol atlas."""

    return (
        ObserverProtocolSpec(
            name="A_lab",
            beta=0.0,
            origin_lab_time=0.0,
            origin_lab_position=0.0,
            beacon_separation=beacon_separation,
        ),
        ObserverProtocolSpec(
            name="B_moving_pos",
            beta=0.3,
            origin_lab_time=0.05,
            origin_lab_position=-0.05,
            beacon_separation=beacon_separation,
        ),
        ObserverProtocolSpec(
            name="C_moving_neg",
            beta=-0.25,
            origin_lab_time=-0.04,
            origin_lab_position=0.04,
            beacon_separation=beacon_separation,
        ),
    )


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
        description=(
            "Validate overlap and affine Lorentz/Poincare transition maps for "
            "a small atlas of oriented observer protocols."
        )
    )
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument(
        "--n-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_VALUES],
    )
    parser.add_argument(
        "--tick-values",
        nargs="+",
        default=[str(value) for value in DEFAULT_TICK_VALUES],
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--beacon-separation", type=float, default=0.15)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--num-invariant-pairs", type=int, default=500)
    args = parser.parse_args()
    return ExperimentConfig(
        T=args.T,
        n_values=_parse_int_values(args.n_values, "N"),
        tick_values=_parse_int_values(args.tick_values, "tick"),
        repetitions=args.repetitions,
        seed=args.seed,
        beacon_separation=args.beacon_separation,
        output_dir=args.output_dir,
        num_invariant_pairs=args.num_invariant_pairs,
    )


def _validate_config(config: ExperimentConfig) -> None:
    if config.T <= 0:
        raise ValueError("T must be positive")
    if any(value <= 0 for value in config.n_values):
        raise ValueError("all n_values must be positive")
    if any(value < 2 for value in config.tick_values):
        raise ValueError("all tick_values must be at least 2")
    if config.repetitions <= 0:
        raise ValueError("repetitions must be positive")
    if config.beacon_separation <= 0:
        raise ValueError("beacon_separation must be positive")
    if config.num_invariant_pairs < 0:
        raise ValueError("num_invariant_pairs must be non-negative")


def _build_protocol_chains(
    protocols: tuple[ObserverProtocolSpec, ...],
    tick_count: int,
    n_events: int,
    T_global: float,
) -> tuple[list[np.ndarray], dict[str, ProtocolChainIndices]]:
    chain_events: list[np.ndarray] = []
    chain_indices: dict[str, ProtocolChainIndices] = {}
    next_index = n_events
    for spec in protocols:
        tau_range = common_safe_tau_range_for_oriented_protocol_1p1(
            T_global,
            spec,
        )
        indices, next_index = append_oriented_protocol_chains(
            chain_events,
            spec,
            tick_count,
            tau_range,
            next_index,
        )
        chain_indices[spec.name] = indices
    return chain_events, chain_indices


def _reconstruct_charts(
    causal_matrix: np.ndarray,
    target_indices: np.ndarray,
    protocols: tuple[ObserverProtocolSpec, ...],
    chain_indices: dict[str, ProtocolChainIndices],
) -> dict[str, ReconstructedChart]:
    charts: dict[str, ReconstructedChart] = {}
    for spec in protocols:
        indices = chain_indices[spec.name]
        charts[spec.name] = reconstruct_oriented_chart_from_order(
            causal_matrix,
            target_indices,
            indices.primary,
            indices.beacon,
            indices.clock_times,
            spec.beacon_separation,
            spec.name,
        )
    return charts


def _chart_event_rows(
    support_events: np.ndarray,
    charts: dict[str, ReconstructedChart],
    protocols: tuple[ObserverProtocolSpec, ...],
    n_events: int,
    tick_count: int,
    repetition: int,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for spec in protocols:
        chart = charts[spec.name]
        true_coords = affine_lab_to_rest_1p1(
            support_events,
            spec.beta,
            spec.origin_lab,
        )
        for row_index, target_index in enumerate(chart.target_indices):
            accessible = bool(chart.accessible[row_index])
            reconstructed = chart.reconstructed_coords[row_index]
            true_coord = true_coords[row_index]
            rows.append(
                {
                    "N": float(n_events),
                    "tick_count": float(tick_count),
                    "repetition": float(repetition),
                    "protocol_name": spec.name,
                    "target_index": float(target_index),
                    "accessible": float(accessible),
                    "T_hat": float(reconstructed[0]),
                    "X_hat": float(reconstructed[1]),
                    "true_protocol_T": float(true_coord[0]),
                    "true_protocol_X": float(true_coord[1]),
                    "T_error": float(reconstructed[0] - true_coord[0])
                    if accessible
                    else float("nan"),
                    "X_error": float(reconstructed[1] - true_coord[1])
                    if accessible
                    else float("nan"),
                }
            )
    return rows


def _fit_transition_row(
    charts: dict[str, ReconstructedChart],
    specs: dict[str, ObserverProtocolSpec],
    source_name: str,
    target_name: str,
    n_events: int,
    tick_count: int,
    repetition: int,
    invariant_pairs: int,
    seed: int,
) -> tuple[dict[str, float | str], tuple[float, np.ndarray] | None]:
    source_chart = charts[source_name]
    target_chart = charts[target_name]
    source_spec = specs[source_name]
    target_spec = specs[target_name]
    overlap = source_chart.accessible & target_chart.accessible
    overlap &= np.all(np.isfinite(source_chart.reconstructed_coords), axis=1)
    overlap &= np.all(np.isfinite(target_chart.reconstructed_coords), axis=1)
    overlap_count = int(np.count_nonzero(overlap))
    expected_beta = expected_relative_beta_1p1(source_spec.beta, target_spec.beta)
    invariant = chart_interval_disagreement(
        source_chart.reconstructed_coords,
        target_chart.reconstructed_coords,
        overlap,
        invariant_pairs,
        seed=seed,
    )

    fit: tuple[float, np.ndarray] | None = None
    fitted_beta = float("nan")
    translation = np.asarray([np.nan, np.nan], dtype=float)
    transition_rmse = float("nan")
    if overlap_count >= 3:
        fitted_beta, translation, transition_rmse = fit_affine_lorentz_beta_grid_1p1(
            source_chart.reconstructed_coords[overlap],
            target_chart.reconstructed_coords[overlap],
        )
        fit = (fitted_beta, translation)

    row = {
        "N": float(n_events),
        "tick_count": float(tick_count),
        "repetition": float(repetition),
        "source_protocol": source_name,
        "target_protocol": target_name,
        "source_beta": float(source_spec.beta),
        "target_beta": float(target_spec.beta),
        "expected_relative_beta": expected_beta,
        "fitted_beta": fitted_beta,
        "fitted_beta_error": fitted_beta - expected_beta
        if np.isfinite(fitted_beta)
        else float("nan"),
        "fitted_translation_t": float(translation[0]),
        "fitted_translation_x": float(translation[1]),
        "transition_rmse": transition_rmse,
        "accessible_overlap_count": float(overlap_count),
        "accessible_overlap_fraction": overlap_count / n_events,
        "invariant_pair_count": invariant["pair_count"],
        "invariant_interval_rmse": invariant["interval_rmse"],
        "invariant_interval_mae": invariant["interval_mae"],
        "invariant_interval_bias": invariant["interval_bias"],
    }
    return row, fit


def _loop_row(
    fits: dict[tuple[str, str], tuple[float, np.ndarray] | None],
    n_events: int,
    tick_count: int,
    repetition: int,
) -> dict[str, float | str]:
    direct = fits.get(("A_lab", "C_moving_neg"))
    first = fits.get(("A_lab", "B_moving_pos"))
    second = fits.get(("B_moving_pos", "C_moving_neg"))
    if direct is None or first is None or second is None:
        beta_composed = float("nan")
        beta_direct = float("nan")
        beta_error = float("nan")
        translation_error = float("nan")
    else:
        beta_composed, translation_composed = compose_affine_lorentz_maps_1p1(
            first[0],
            first[1],
            second[0],
            second[1],
        )
        beta_direct = direct[0]
        beta_error = beta_composed - beta_direct
        translation_error = float(np.linalg.norm(translation_composed - direct[1]))

    return {
        "N": float(n_events),
        "tick_count": float(tick_count),
        "repetition": float(repetition),
        "loop": "A_lab->B_moving_pos->C_moving_neg",
        "beta_composed": beta_composed,
        "beta_direct": beta_direct,
        "beta_composition_error": beta_error,
        "translation_composition_error_norm": translation_error,
    }


def run_experiment(
    config: ExperimentConfig,
) -> tuple[
    list[dict[str, float | str]],
    list[dict[str, float | str]],
    list[dict[str, float | str]],
]:
    """Run observer-atlas consistency validation."""

    _validate_config(config)
    protocols = default_protocols(config.beacon_separation)
    specs = {spec.name: spec for spec in protocols}
    chart_rows: list[dict[str, float | str]] = []
    transition_rows: list[dict[str, float | str]] = []
    loop_rows: list[dict[str, float | str]] = []

    for n_index, n_events in enumerate(config.n_values):
        for repetition in range(config.repetitions):
            support_events = sprinkle_1p1_causal_diamond(
                n_events,
                T=config.T,
                seed=config.seed + 100_000 * n_index + 1_000 * repetition,
            )
            target_indices = np.arange(n_events, dtype=int)
            for tick_count in config.tick_values:
                chain_events, chain_indices = _build_protocol_chains(
                    protocols,
                    tick_count,
                    n_events,
                    config.T,
                )
                combined_events = np.vstack((support_events, *chain_events))
                causal_matrix = causal_matrix_1p1(combined_events)
                charts = _reconstruct_charts(
                    causal_matrix,
                    target_indices,
                    protocols,
                    chain_indices,
                )
                chart_rows.extend(
                    _chart_event_rows(
                        support_events,
                        charts,
                        protocols,
                        n_events,
                        tick_count,
                        repetition,
                    )
                )
                fits: dict[tuple[str, str], tuple[float, np.ndarray] | None] = {}
                for transition_index, (source_name, target_name) in enumerate(
                    DEFAULT_TRANSITIONS
                ):
                    row, fit = _fit_transition_row(
                        charts,
                        specs,
                        source_name,
                        target_name,
                        n_events,
                        tick_count,
                        repetition,
                        config.num_invariant_pairs,
                        seed=(
                            config.seed
                            + 10_000 * n_index
                            + 100 * repetition
                            + transition_index
                        ),
                    )
                    transition_rows.append(row)
                    fits[(source_name, target_name)] = fit
                loop_rows.append(_loop_row(fits, n_events, tick_count, repetition))

    return chart_rows, transition_rows, loop_rows


def _write_csv(rows: list[dict[str, float | str]], output_path: Path) -> Path:
    if not rows:
        raise RuntimeError("no rows to write")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def write_outputs(
    chart_rows: list[dict[str, float | str]],
    transition_rows: list[dict[str, float | str]],
    loop_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    """Write observer-atlas CSV outputs."""

    data_dir = output_dir / "data"
    chart_path = data_dir / "observer_atlas_chart_events.csv"
    transition_path = data_dir / "observer_atlas_transition_summary.csv"
    loop_path = data_dir / "observer_atlas_loop_summary.csv"
    return (
        _write_csv(chart_rows, chart_path),
        _write_csv(transition_rows, transition_path),
        _write_csv(loop_rows, loop_path),
    )


def _transition_label(row: dict[str, float | str]) -> str:
    return f"{row['source_protocol']}->{row['target_protocol']}"


def _mean_by_tick(
    rows: list[dict[str, float | str]],
    value_key: str,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    labels = sorted({_transition_label(row) for row in rows})
    result: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for label in labels:
        label_rows = [row for row in rows if _transition_label(row) == label]
        ticks = sorted({int(row["tick_count"]) for row in label_rows})
        values: list[float] = []
        for tick in ticks:
            tick_values = np.asarray(
                [
                    float(row[value_key])
                    for row in label_rows
                    if int(row["tick_count"]) == tick
                ],
                dtype=float,
            )
            values.append(float(np.nanmean(tick_values)))
        result[label] = (np.asarray(ticks, dtype=float), np.asarray(values))
    return result


def save_transition_beta_error_vs_ticks(
    transition_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save fitted relative beta error versus tick count."""

    output_path = (
        output_dir
        / "figures"
        / "observer_atlas_transition_beta_error_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, (ticks, values) in _mean_by_tick(
        transition_rows,
        "fitted_beta_error",
    ).items():
        ax.plot(ticks, values, marker="o", label=label)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Fitted beta minus expected beta")
    ax.set_title("Atlas transition beta error")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_transition_rmse_vs_ticks(
    transition_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save affine transition residual RMSE versus tick count."""

    output_path = output_dir / "figures" / "observer_atlas_transition_rmse_vs_ticks.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, (ticks, values) in _mean_by_tick(
        transition_rows,
        "transition_rmse",
    ).items():
        ax.plot(ticks, values, marker="o", label=label)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Affine transition RMSE")
    ax.set_title("Atlas transition residual")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_invariant_disagreement_vs_ticks(
    transition_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save invariant interval disagreement versus tick count."""

    output_path = (
        output_dir
        / "figures"
        / "observer_atlas_invariant_disagreement_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, (ticks, values) in _mean_by_tick(
        transition_rows,
        "invariant_interval_rmse",
    ).items():
        ax.plot(ticks, values, marker="o", label=label)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Invariant interval RMSE")
    ax.set_title("Cross-chart invariant interval disagreement")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_overlap_fraction_vs_ticks(
    transition_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save accessible overlap fraction versus tick count."""

    output_path = (
        output_dir / "figures" / "observer_atlas_overlap_fraction_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, (ticks, values) in _mean_by_tick(
        transition_rows,
        "accessible_overlap_fraction",
    ).items():
        ax.plot(ticks, values, marker="o", label=label)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Accessible overlap fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Overlap of reconstructed charts")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_loop_consistency_vs_ticks(
    loop_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Save loop consistency errors versus tick count."""

    output_path = (
        output_dir / "figures" / "observer_atlas_loop_consistency_vs_ticks.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ticks = sorted({int(row["tick_count"]) for row in loop_rows})
    beta_errors: list[float] = []
    translation_errors: list[float] = []
    for tick in ticks:
        rows = [row for row in loop_rows if int(row["tick_count"]) == tick]
        beta_errors.append(
            float(
                np.nanmean(
                    [abs(float(row["beta_composition_error"])) for row in rows]
                )
            )
        )
        translation_errors.append(
            float(
                np.nanmean(
                    [
                        float(row["translation_composition_error_norm"])
                        for row in rows
                    ]
                )
            )
        )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(ticks, beta_errors, marker="o", label="|beta composition error|")
    ax.plot(
        ticks,
        translation_errors,
        marker="s",
        label="translation composition error norm",
    )
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Observer tick count")
    ax.set_ylabel("Loop consistency error")
    ax.set_title("Atlas loop consistency")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path


def save_figures(
    transition_rows: list[dict[str, float | str]],
    loop_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path, Path, Path, Path]:
    """Save observer-atlas figures."""

    return (
        save_transition_beta_error_vs_ticks(transition_rows, output_dir),
        save_transition_rmse_vs_ticks(transition_rows, output_dir),
        save_invariant_disagreement_vs_ticks(transition_rows, output_dir),
        save_overlap_fraction_vs_ticks(transition_rows, output_dir),
        save_loop_consistency_vs_ticks(loop_rows, output_dir),
    )


def main() -> None:
    config = parse_args()
    chart_rows, transition_rows, loop_rows = run_experiment(config)
    chart_path, transition_path, loop_path = write_outputs(
        chart_rows,
        transition_rows,
        loop_rows,
        config.output_dir,
    )
    figure_paths = save_figures(transition_rows, loop_rows, config.output_dir)
    print(f"Wrote chart-event results: {chart_path}")
    print(f"Wrote transition summary: {transition_path}")
    print(f"Wrote loop summary: {loop_path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")
    print(
        "Atlas transition maps use supplied observer origins, clocks, and "
        "beacon separations; hidden coordinates are used only for validation."
    )


if __name__ == "__main__":
    main()
