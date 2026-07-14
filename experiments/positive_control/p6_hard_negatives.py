"""P6a: chain-rich hard negatives for the frozen P3 geometry discriminator.

Stage A maps the constructed layered family and audits the proposed local-
shuffle family. Stage B is intentionally locked until Stage A results and
expectations have been frozen in ``p6_test_constants.json``.

The local-shuffle audit scores the same fitted embedding against both the
pre-shuffle and post-shuffle exact 2D-order coordinates. This distinguishes
loss of an externally retained coordinate label from loss of intrinsic 2D
order geometry.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable
from pathlib import Path

import numpy as np
from p3_dynamics import _median, analyze_order
from pc_common import git_describe, parse_seed_spec, write_rows_csv

from causal_spacetime_lab.positive_control.two_orders import (
    balanced_layered_perm,
    chain_observables,
    perm_to_causal_matrix,
    windowed_transpositions,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"
FROZEN = ROOT / "docs" / "prereg" / "frozen" / "p6_test_constants.json"

N_ELEMENTS = 600
GATE_HELDOUT = 0.10
GATE_NULL_GAP = 0.10
GATE_TRUTH = 0.40

LAYER_COUNTS = (25, 40, 60)
LAYER_JITTER_MOVES = (20, 100, 500, 2000)
LAYER_JITTER_WINDOW = 60
LOCAL_SHUFFLE_MOVES = (30, 150, 600, 2400)
LOCAL_SHUFFLE_WINDOW = 32


def order_inputs(pi: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return causal order, topological time, and exact spatial coordinate."""
    idx = np.arange(pi.size, dtype=float)
    return perm_to_causal_matrix(pi), idx + pi, idx - pi


def gate_pass(row: dict, truth_key: str = "truth") -> bool:
    """Apply the inherited frozen P3 gates to one successful analysis row."""
    return bool(
        row.get("status") == "ok"
        and row["heldout"] <= GATE_HELDOUT
        and row["null_gap"] >= GATE_NULL_GAP
        and row[truth_key] <= GATE_TRUTH
    )


def _with_observables(row: dict, causal: np.ndarray) -> dict:
    return {**row, **chain_observables(causal)}


def analyze_layered(layer_count: int, moves: int, seed: int) -> dict:
    """Analyze one jittered complete-layer construction."""
    base = balanced_layered_perm(
        N_ELEMENTS, layer_count=layer_count, seed=seed, min_layer_size=6
    )
    pi = windowed_transpositions(
        base, moves=moves, window=LAYER_JITTER_WINDOW, seed=10_000 + seed
    )
    causal, times, coords = order_inputs(pi)
    row = analyze_order(causal, times, coords, seed=seed, want_truth=True)
    row = _with_observables(row, causal)
    row.update(
        {
            "family": "layered",
            "layer_count": float(layer_count),
            "moves": float(moves),
            "window": float(LAYER_JITTER_WINDOW),
            "seed": float(seed),
            "gate_pass": float(gate_pass(row)),
        }
    )
    return row


def analyze_local_shuffle(moves: int, seed: int) -> dict:
    """Audit a local shuffle of a uniform continuum 2D order.

    Uniform permutations are the beta=0 continuum ensemble and are used only
    for the Stage-A construction audit. The exact argument applies equally to
    P5 beta=2/8 states: transposing a permutation produces another exact 2D
    order, so the operation cannot by itself construct a geometry-free order.
    """
    base = np.random.default_rng(seed).permutation(N_ELEMENTS)
    pi = windowed_transpositions(
        base, moves=moves, window=LOCAL_SHUFFLE_WINDOW, seed=20_000 + seed
    )
    causal, times, coords = order_inputs(pi)
    idx = np.arange(N_ELEMENTS, dtype=float)
    original_coords = idx - base
    row = analyze_order(
        causal,
        times,
        original_coords,
        seed=seed,
        want_truth=True,
        extra_truth_coords={"current": coords},
    )
    row = _with_observables(row, causal)
    row.update(
        {
            "family": "local_shuffle_audit",
            "moves": float(moves),
            "window": float(LOCAL_SHUFFLE_WINDOW),
            "seed": float(seed),
            "gate_pass_original": float(gate_pass(row)),
            "gate_pass_current": float(gate_pass(row, "truth_current")),
        }
    )
    return row


def _cell_summary(rows: list[dict], keys: tuple[str, ...]) -> list[dict]:
    cells: list[dict] = []
    identities = sorted({tuple(row[key] for key in keys) for row in rows})
    for identity in identities:
        subset = [row for row in rows if tuple(row[key] for key in keys) == identity]
        ok = [row for row in subset if row["status"] == "ok"]
        cell = {key: value for key, value in zip(keys, identity, strict=True)}
        cell.update(
            {
                "n_total": len(subset),
                "n_valid": len(ok),
                "n_gate_pass": sum(int(row.get("gate_pass", 0)) for row in ok),
                "n_gate_pass_original": sum(
                    int(row.get("gate_pass_original", 0)) for row in ok
                ),
                "n_gate_pass_current": sum(
                    int(row.get("gate_pass_current", 0)) for row in ok
                ),
                "median_heldout": _median([row["heldout"] for row in ok]),
                "median_null_gap": _median([row["null_gap"] for row in ok]),
                "median_truth": _median([row["truth"] for row in ok]),
                "median_truth_current": _median(
                    [row["truth_current"] for row in ok if "truth_current" in row]
                ),
                "median_height": _median([row["height"] for row in subset]),
            }
        )
        cells.append(cell)
    return cells


def stage_a(seeds: Iterable[int]) -> None:
    """Run exploratory calibration and emit a mechanical freeze proposal."""
    seeds = list(seeds)
    OUT.mkdir(parents=True, exist_ok=True)
    code_version = git_describe()
    layered: list[dict] = []
    local: list[dict] = []
    for seed in seeds:
        for layer_count in LAYER_COUNTS:
            for moves in LAYER_JITTER_MOVES:
                row = analyze_layered(layer_count, moves, seed)
                row.update({"stage": "P6-A", "code_version": code_version})
                layered.append(row)
                print(
                    f"layered seed={seed} k={layer_count} moves={moves}: "
                    f"{row['status']} pass={int(row['gate_pass'])}",
                    flush=True,
                )
        for moves in LOCAL_SHUFFLE_MOVES:
            row = analyze_local_shuffle(moves, seed)
            row.update({"stage": "P6-A", "code_version": code_version})
            local.append(row)
            print(
                f"local seed={seed} moves={moves}: {row['status']} "
                f"old={int(row['gate_pass_original'])} "
                f"current={int(row['gate_pass_current'])}",
                flush=True,
            )

    write_rows_csv(OUT / "p6_stage_a_layered.csv", layered)
    write_rows_csv(OUT / "p6_stage_a_local_shuffle_audit.csv", local)
    layered_cells = _cell_summary(layered, ("layer_count", "moves"))
    local_cells = _cell_summary(local, ("moves",))

    # A candidate must exercise the discriminator rather than fail chain
    # extraction: >= 8/10 valid and >= 8/10 valid gate blocks.
    proposed_layered = [
        {"layer_count": int(cell["layer_count"]), "moves": int(cell["moves"])}
        for cell in layered_cells
        if cell["n_valid"] >= 8 and cell["n_valid"] - cell["n_gate_pass"] >= 8
    ]
    summary = {
        "stage": "P6-A",
        "code_version": code_version,
        "n_elements": N_ELEMENTS,
        "seeds": list(seeds),
        "gates_inherited_from": "P3",
        "layered_cells": layered_cells,
        "proposed_layered_confirmatory_cells": proposed_layered,
        "local_shuffle_cells": local_cells,
        "local_shuffle_disposition": (
            "construction audit only: every output is still an exact 2D order; "
            "a current-coordinate pass with an original-coordinate failure is "
            "coordinate remapping, not a geometry-free hard negative"
        ),
    }
    path = OUT / "p6_stage_a_summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    print(f"Freeze selected cells in {FROZEN} before Stage B.", flush=True)


def _load_frozen() -> dict:
    if not FROZEN.exists():
        raise SystemExit(
            f"frozen P6 constants not found at {FROZEN}; run Stage A, inspect "
            "the calibration, and freeze expectations before Stage B"
        )
    return json.loads(FROZEN.read_text(encoding="utf-8"))


def stage_b(seeds: Iterable[int]) -> None:
    """Apply frozen layered expectations on fresh construction seeds."""
    frozen = _load_frozen()
    expected_seeds = parse_seed_spec(frozen["confirmatory_seeds"])
    seeds = list(seeds)
    if seeds != expected_seeds:
        raise SystemExit(
            f"Stage B seeds are frozen as {frozen['confirmatory_seeds']}; got {seeds}"
        )
    OUT.mkdir(parents=True, exist_ok=True)
    code_version = git_describe()
    rows: list[dict] = []
    for seed in seeds:
        for cell in frozen["layered_cells"]:
            row = analyze_layered(cell["layer_count"], cell["moves"], seed)
            row.update({"stage": "P6-B", "code_version": code_version})
            rows.append(row)
            print(
                f"seed={seed} k={cell['layer_count']} moves={cell['moves']}: "
                f"{row['status']} pass={int(row['gate_pass'])}",
                flush=True,
            )
    write_rows_csv(OUT / "p6_stage_b_layered.csv", rows)
    decide_stage_b(rows, frozen, code_version)


def decide_stage_b(rows: list[dict], frozen: dict, code_version: str) -> None:
    cells = _cell_summary(rows, ("layer_count", "moves"))
    denominator = len(parse_seed_spec(frozen["confirmatory_seeds"]))
    minimum = int(frozen["minimum_chain_rich_blocks"])
    for cell in cells:
        cell["expectation"] = "chain-rich block"
        cell["expectation_met"] = bool(
            cell["n_valid"] >= minimum
            and cell["n_valid"] - cell["n_gate_pass"] >= minimum
        )
    registry = {
        "stage": "P6-B",
        "code_version": code_version,
        "confirmatory_denominator": denominator,
        "minimum_chain_rich_blocks": minimum,
        "cells": cells,
        "all_expectations_met": bool(cells)
        and all(cell["expectation_met"] for cell in cells),
    }
    path = OUT / "p6_stage_b_decision_registry.json"
    path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(registry, indent=2), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("a", "b"), required=True)
    parser.add_argument("--seeds", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    default = "0-9" if args.stage == "a" else "100-119"
    seeds = parse_seed_spec(args.seeds or default)
    if args.stage == "a":
        stage_a(seeds)
    else:
        stage_b(seeds)


if __name__ == "__main__":
    main()
