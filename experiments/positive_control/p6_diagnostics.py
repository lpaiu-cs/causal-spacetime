"""P6b same-data comparison of the frozen instrument and cheap diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from pc_common import (
    git_describe,
    parse_seed_spec,
    primary_scene_config,
    write_rows_csv,
)

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.positive_control.dynamics import transitive_percolation
from causal_spacetime_lab.positive_control.epsilon_sweep import epsilon_diluted_order
from causal_spacetime_lab.positive_control.rewire import geometric_post_closure_density
from causal_spacetime_lab.positive_control.scene import build_positive_control_scene
from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    chain_observables,
    mcmc_2d_order_fast,
    perm_to_causal_matrix,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"
FROZEN = ROOT / "docs" / "prereg" / "frozen"
CONSTANTS = FROZEN / "p6b_test_constants.json"
STRING_FIELDS = {"source", "condition", "reference_group", "status", "code_version"}
RANK_DIAGNOSTICS = (
    "instrument_margin",
    "mm_distance",
    "abundance_distance",
    "height_distance",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _normalize_p5_row(row: dict[str, str]) -> dict[str, str]:
    """Repair legacy short crystal rows in the frozen concatenated P5 CSV."""
    if row.get("seed") is not None:
        return row
    if row.get("beta") != "32.0" or not str(row.get("n_targets", "")).startswith(
        "structural_block:"
    ):
        raise ValueError("unrecognized malformed P5 frozen row")
    return {
        **row,
        "sample": row["heldout"],
        "seed": row["min_chain_len"],
        "status": row["n_targets"],
        "heldout": "",
        "min_chain_len": "",
        "n_targets": "",
    }


def _float(row: dict, key: str) -> float | None:
    value = row.get(key)
    return float(value) if value not in (None, "") else None


def _digest(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()[:16]


def _instrument_margin(source: str, row: dict) -> float:
    if row["status"] != "ok":
        return -5.0
    heldout = _float(row, "heldout")
    if heldout is None:
        heldout = _float(row, "heldout_violation")
    if heldout is None:
        raise ValueError(f"{source} row has no heldout value")
    truth = _float(row, "truth")
    if truth is None:
        truth = _float(row, "truth_order_error")
    null_gap = _float(row, "null_gap")
    if source == "P1":
        assert truth is not None
        return min((0.15 - truth) / 0.15, (0.05 - heldout) / 0.05)
    margins = [(0.10 - heldout) / 0.10]
    if null_gap is not None:
        margins.append((null_gap - 0.10) / 0.10)
    if truth is not None:
        margins.append((0.40 - truth) / 0.40)
    return min(margins)


def _base_row(
    source: str,
    condition: str,
    reference_group: str,
    seed: int,
    causal: np.ndarray,
) -> dict:
    return {
        "source": source,
        "condition": condition,
        "reference_group": reference_group,
        "seed": float(seed),
        **chain_observables(causal),
        "code_version": git_describe(),
    }


def collect_references() -> None:
    constants = json.loads(CONSTANTS.read_text(encoding="utf-8"))
    rows: list[dict] = []
    for seed in parse_seed_spec(constants["reference_seeds_p1_scene"]):
        scene = build_positive_control_scene(primary_scene_config(seed))
        rows.append(_base_row("reference", "p1_scene", "p1_scene", seed, scene.causal))
    for seed in parse_seed_spec(constants["reference_seeds_p3_sprinkled"]):
        events = sprinkle_1p1_causal_diamond(1500, T=2.0, seed=seed)
        causal = causal_matrix_1p1(events)
        rows.append(_base_row("reference", "p3_sprinkled", "p3_n1500", seed, causal))
    for seed in parse_seed_spec(constants["reference_seeds_p5_uniform"]):
        pi = np.random.default_rng(seed).permutation(600)
        rows.append(
            _base_row(
                "reference",
                "p5_uniform",
                "orders_n600",
                seed,
                perm_to_causal_matrix(pi),
            )
        )
    write_rows_csv(OUT / "p6b_raw_references.csv", rows)


def collect_p1() -> None:
    frozen = _read_csv(FROZEN / "p1_stage_b_epsilon_sweep.csv")
    lookup = {(int(float(r["seed"])), float(r["epsilon"])): r for r in frozen}
    epsilons = sorted({float(row["epsilon"]) for row in frozen})
    rows: list[dict] = []
    for seed in sorted({int(float(row["seed"])) for row in frozen}):
        scene = build_positive_control_scene(primary_scene_config(seed))
        target_density, _, _ = geometric_post_closure_density(scene)
        for epsilon in epsilons:
            if epsilon == 0.0:
                causal = scene.causal
                density = target_density
            else:
                causal, density, _ = epsilon_diluted_order(
                    scene, epsilon, seed + 41, target_density
                )
            old = lookup[(seed, epsilon)]
            if old["status"] == "ok" and _digest(causal) != old["causal_digest"]:
                raise RuntimeError(
                    f"P1 causal digest mismatch at seed={seed}, eps={epsilon}"
                )
            if epsilon <= 0.15:
                label: float | str = 1.0
            elif epsilon >= 0.90:
                label = 0.0
            else:
                label = ""
            row = _base_row("P1", f"epsilon={epsilon:g}", "p1_scene", seed, causal)
            row.update(
                {
                    "epsilon": epsilon,
                    "label": label,
                    "status": old["status"],
                    "heldout": old.get("heldout_violation", ""),
                    "truth": old.get("truth_order_error", ""),
                    "instrument_margin": _instrument_margin("P1", old),
                    "achieved_density": density,
                }
            )
            rows.append(row)
        print(f"P1 seed={seed}: {len(epsilons)} orders", flush=True)
    write_rows_csv(OUT / "p6b_raw_p1.csv", rows)


def collect_p3() -> None:
    frozen = _read_csv(FROZEN / "p3_stage_b_dynamics.csv")
    lookup = {(int(float(r["seed"])), float(r["p"])): r for r in frozen}
    rows: list[dict] = []
    for p in sorted({float(row["p"]) for row in frozen}):
        for seed in sorted({int(float(row["seed"])) for row in frozen}):
            causal, _ = transitive_percolation(1500, p, seed)
            old = lookup[(seed, p)]
            row = _base_row("P3", f"p={p:g}", "p3_n1500", seed, causal)
            row.update(
                {
                    "p": p,
                    "label": 0.0,
                    "status": old["status"],
                    "heldout": old.get("heldout", ""),
                    "null_gap": old.get("null_gap", ""),
                    "instrument_margin": _instrument_margin("P3", old),
                }
            )
            rows.append(row)
        print(f"P3 p={p:g}: 20 orders", flush=True)
    write_rows_csv(OUT / "p6b_raw_p3.csv", rows)


def collect_p6() -> None:
    frozen = _read_csv(FROZEN / "p6_stage_b_layered.csv")
    rows: list[dict] = []
    for old in frozen:
        layer_count = int(float(old["layer_count"]))
        moves = int(float(old["moves"]))
        row = {
            "source": "P6",
            "condition": f"k={layer_count},moves={moves}",
            "reference_group": "orders_n600",
            "seed": old["seed"],
            "layer_count": old["layer_count"],
            "moves": old["moves"],
            "label": 0.0,
            "status": old["status"],
            "heldout": old.get("heldout", ""),
            "null_gap": old.get("null_gap", ""),
            "truth": old.get("truth", ""),
            "instrument_margin": _instrument_margin("P6", old),
            "R": old["R"],
            "n0": old["n0"],
            "n1": old["n1"],
            "n2": old["n2"],
            "mm_dim": old["mm_dim"],
            "height": old["height"],
            "code_version": git_describe(),
        }
        rows.append(row)
    write_rows_csv(OUT / "p6b_raw_p6.csv", rows)


def collect_p5(beta: float, seed: int, accelerated: bool = False) -> None:
    constants = json.loads((FROZEN / "p5_test_constants.json").read_text())
    frozen = [
        _normalize_p5_row(row) for row in _read_csv(FROZEN / "p5_stage_b_all.csv")
    ]
    expected = [
        row for row in frozen if float(row["beta"]) == beta and int(row["seed"]) == seed
    ]
    expected.sort(key=lambda row: int(row["sample"]))
    if not expected:
        raise SystemExit(f"no frozen P5 cells for beta={beta:g}, seed={seed}")
    crystal = beta == float(constants["crystal_beta"])
    steps = constants["steps_crystal" if crystal else "steps_continuum"]
    sample_every = constants[
        "sample_every_crystal" if crystal else "sample_every_continuum"
    ]
    if crystal:
        pi0 = bipartite_perm(600)
    else:
        pi0 = np.random.default_rng(seed * 31 + int(round(beta * 10))).permutation(600)
    sampler = mcmc_2d_order_fast
    if accelerated:
        from causal_spacetime_lab.positive_control.accelerated_two_orders import (
            mcmc_2d_order_replay_accelerated,
        )

        sampler = mcmc_2d_order_replay_accelerated
    samples, acceptance, _ = sampler(
        pi0,
        beta=beta,
        eps=0.02,
        steps=steps,
        seed=seed * 1000 + int(round(beta * 10)),
        sample_every=sample_every,
        burn_frac=constants["burn_frac"],
    )
    if len(samples) != len(expected):
        raise RuntimeError(
            f"P5 sample count mismatch: {len(samples)} != {len(expected)}"
        )
    rows: list[dict] = []
    for sample_index, (obs, old) in enumerate(zip(samples, expected, strict=True)):
        row = {
            "source": "P5",
            "condition": f"beta={beta:g}",
            "reference_group": "orders_n600",
            "seed": float(seed),
            "beta": beta,
            "sample": float(sample_index),
            "label": float(not crystal),
            "status": old["status"],
            "heldout": old.get("heldout", ""),
            "null_gap": old.get("null_gap", ""),
            "truth": old.get("truth", ""),
            "instrument_margin": _instrument_margin("P5", old),
            "acceptance": acceptance,
            "code_version": git_describe(),
            **obs,
        }
        rows.append(row)
    chain_s = float(np.mean([row["S"] for row in samples]))
    expected_s = float(expected[0]["chain_S"])
    if not math.isclose(chain_s, expected_s, rel_tol=0.0, abs_tol=1e-9):
        raise RuntimeError(f"P5 chain provenance mismatch: {chain_s} != {expected_s}")
    write_rows_csv(OUT / f"p6b_raw_p5_b{beta:g}_s{seed}.csv", rows)


def _average_ranks(values: list[float]) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1) + 1.0
        start = stop
    return ranks


def roc_auc(labels: list[int], scores: list[float]) -> float:
    ranks = _average_ranks(scores)
    positives = np.asarray(labels, dtype=int) == 1
    n_pos = int(positives.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    return float((ranks[positives].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def spearman(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 3:
        return float("nan")
    rx = _average_ranks(xs)
    ry = _average_ranks(ys)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _reference_models(rows: list[dict], quantile: float) -> tuple[dict, dict]:
    refs = [row for row in rows if row["source"] == "reference"]
    models: dict[str, dict] = {}
    distances: dict[str, list[float]] = defaultdict(list)
    for group in sorted({row["reference_group"] for row in refs}):
        subset = [row for row in refs if row["reference_group"] == group]
        vectors = np.array(
            [
                [
                    row["mm_dim"],
                    np.log1p(row["n0"]),
                    np.log1p(row["n1"]),
                    np.log1p(row["n2"]),
                    np.log1p(row["height"]),
                ]
                for row in subset
            ],
            dtype=float,
        )
        models[group] = {
            "mean": vectors.mean(axis=0),
            "std": np.maximum(vectors.std(axis=0, ddof=1), 1e-12),
        }
    for row in refs:
        scored = _score_cheap(row, models[row["reference_group"]])
        for key, value in scored.items():
            distances[f"{row['reference_group']}:{key}"].append(value)
    cutoffs = {
        key: float(np.quantile(values, quantile)) for key, values in distances.items()
    }
    return models, cutoffs


def _score_cheap(row: dict, model: dict) -> dict[str, float]:
    vector = np.array(
        [
            float(row["mm_dim"]),
            np.log1p(float(row["n0"])),
            np.log1p(float(row["n1"])),
            np.log1p(float(row["n2"])),
            np.log1p(float(row["height"])),
        ]
    )
    z = (vector - model["mean"]) / model["std"]
    return {
        "mm_distance": float(abs(z[0])),
        "abundance_distance": float(np.sqrt(np.mean(z[1:4] ** 2))),
        "height_distance": float(abs(z[4])),
    }


def _rank_summary(rank_rows: list[dict], diagnostic: str) -> dict[str, float]:
    subset = [row for row in rank_rows if row["diagnostic"] == diagnostic]
    result: dict[str, float] = {}
    for key in ("rho_epsilon", "rho_truth"):
        values = np.array([row[key] for row in subset], dtype=float)
        result[f"median_{key}"] = float(np.nanmedian(values))
        result[f"min_{key}"] = float(np.nanmin(values))
        result[f"max_{key}"] = float(np.nanmax(values))
    return result


def _write_comparison_deliverables(summary: dict) -> None:
    table_rows = []
    for diagnostic in RANK_DIAGNOSTICS:
        overlap = summary["p1_h_lag_overlap"].get(diagnostic, {})
        table_rows.append(
            {
                "diagnostic": diagnostic,
                "roc_auc": summary["roc_auc"][diagnostic],
                **summary["p1_rank_summary"][diagnostic],
                "h_lag_false_pass_fraction": overlap.get("fraction", ""),
                "h_lag_false_pass_count": overlap.get("n_false_pass", ""),
                "h_lag_cell_count": overlap.get("n_h_lag", ""),
            }
        )
    write_rows_csv(OUT / "p6b_diagnostics_table.csv", table_rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["Instrument", "MM dim.", "Abundance", "Height"]
    colors = ["#1f2937", "#3b82f6", "#10b981", "#f59e0b"]
    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.4))
    axes[0].bar(
        labels,
        [summary["roc_auc"][name] for name in RANK_DIAGNOSTICS],
        color=colors,
    )
    axes[0].axhline(0.5, color="#6b7280", linestyle="--", linewidth=1)
    axes[0].set_ylabel("ROC AUC")
    axes[0].set_ylim(0, 1.03)
    axes[0].set_title("Frozen class labels")

    axes[1].bar(
        labels,
        [
            summary["p1_rank_summary"][name]["median_rho_epsilon"]
            for name in RANK_DIAGNOSTICS
        ],
        color=colors,
    )
    axes[1].axhline(0, color="#6b7280", linewidth=1)
    axes[1].set_ylabel("Median Spearman rho")
    axes[1].set_ylim(-1.03, 1.03)
    axes[1].set_title("P1 response vs epsilon")

    cheap = RANK_DIAGNOSTICS[1:]
    axes[2].bar(
        labels[1:],
        [summary["p1_h_lag_overlap"][name]["fraction"] for name in cheap],
        color=colors[1:],
    )
    axes[2].set_ylabel("False-pass overlap")
    axes[2].set_ylim(0, 1.03)
    axes[2].set_title("P1 H-LAG cells")
    for axis in axes:
        axis.tick_params(axis="x", rotation=25)
        axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(OUT / "p6b_diagnostics_comparison.png", dpi=180)
    figure.savefig(OUT / "p6b_diagnostics_comparison.pdf")
    plt.close(figure)


def aggregate() -> None:
    required = [
        OUT / "p6b_raw_references.csv",
        OUT / "p6b_raw_p1.csv",
        OUT / "p6b_raw_p3.csv",
        OUT / "p6b_raw_p6.csv",
    ]
    p5_constants = json.loads((FROZEN / "p5_test_constants.json").read_text())
    p5_cells = [(beta, seed) for beta in (2.0, 8.0) for seed in (100, 101, 102)]
    p5_cells += [(float(p5_constants["crystal_beta"]), seed) for seed in (100, 101)]
    required += [OUT / f"p6b_raw_p5_b{beta:g}_s{seed}.csv" for beta, seed in p5_cells]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        message = "P6b aggregation missing frozen inputs:\n" + "\n".join(missing)
        raise SystemExit(message)
    rows: list[dict] = []
    for path in required:
        for row in _read_csv(path):
            rows.append(
                {
                    key: (
                        float(value)
                        if key not in STRING_FIELDS and value != ""
                        else value
                    )
                    for key, value in row.items()
                }
            )
    constants = json.loads(CONSTANTS.read_text(encoding="utf-8"))
    models, cutoffs = _reference_models(
        rows, float(constants["reference_pass_quantile"])
    )
    data = [row for row in rows if row["source"] != "reference"]
    for row in data:
        row.update(_score_cheap(row, models[row["reference_group"]]))
    labelled = [row for row in data if row.get("label", "") != ""]
    labels = [int(row["label"]) for row in labelled]
    diagnostics = {
        "instrument_margin": [float(row["instrument_margin"]) for row in labelled],
        "mm_distance": [-row["mm_distance"] for row in labelled],
        "abundance_distance": [-row["abundance_distance"] for row in labelled],
        "height_distance": [-row["height_distance"] for row in labelled],
    }
    auc = {name: roc_auc(labels, scores) for name, scores in diagnostics.items()}

    p1 = [row for row in data if row["source"] == "P1" and row["status"] == "ok"]
    rank_rows: list[dict] = []
    for seed in sorted({int(row["seed"]) for row in p1}):
        subset = sorted(
            (row for row in p1 if int(row["seed"]) == seed),
            key=lambda row: row["epsilon"],
        )
        for name in RANK_DIAGNOSTICS:
            values = [
                -row[name] if name == "instrument_margin" else row[name]
                for row in subset
            ]
            rank_rows.append(
                {
                    "seed": seed,
                    "diagnostic": name,
                    "rho_epsilon": spearman([row["epsilon"] for row in subset], values),
                    "rho_truth": spearman([row["truth"] for row in subset], values),
                    "n_cells": len(subset),
                }
            )
    h_lag = [row for row in p1 if row["heldout"] <= 0.05 and row["truth"] > 0.15]
    overlap = {}
    for name in RANK_DIAGNOSTICS[1:]:
        overlap[name] = {
            "n_h_lag": len(h_lag),
            "n_false_pass": sum(
                row[name] <= cutoffs[f"p1_scene:{name}"] for row in h_lag
            ),
        }
        overlap[name]["fraction"] = (
            overlap[name]["n_false_pass"] / len(h_lag) if h_lag else None
        )

    write_rows_csv(OUT / "p6b_scored_rows.csv", data)
    write_rows_csv(OUT / "p6b_p1_rank_correlations.csv", rank_rows)
    summary = {
        "code_version": git_describe(),
        "n_rows": len(data),
        "n_labelled": len(labelled),
        "class_counts": {
            "geometric": sum(labels),
            "nongeometric": len(labels) - sum(labels),
        },
        "roc_auc": auc,
        "reference_cutoffs": cutoffs,
        "p1_h_lag_overlap": overlap,
        "p1_rank_summary": {
            name: _rank_summary(rank_rows, name) for name in diagnostics
        },
    }
    (OUT / "p6b_diagnostics_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _write_comparison_deliverables(summary)
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=("references", "p1", "p3", "p5", "p6", "aggregate"),
        required=True,
    )
    parser.add_argument("--beta", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument(
        "--accelerated",
        action="store_true",
        help="use optional trajectory-equivalent Numba replay for P5",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.source == "references":
        collect_references()
    elif args.source == "p1":
        collect_p1()
    elif args.source == "p3":
        collect_p3()
    elif args.source == "p5":
        if args.beta is None or args.seed is None:
            raise SystemExit("--source p5 requires --beta and --seed")
        collect_p5(args.beta, args.seed, accelerated=args.accelerated)
    elif args.source == "p6":
        collect_p6()
    else:
        aggregate()


if __name__ == "__main__":
    main()
