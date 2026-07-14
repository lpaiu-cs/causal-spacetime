"""P7 Stage A validation of the continuous geometry order parameter G."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from pc_common import git_describe, write_rows_csv

from causal_spacetime_lab.positive_control.geometry_score import (
    clipped_gate_margin_score,
    geometry_order_parameter,
)

ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "docs" / "prereg" / "frozen"
OUT = ROOT / "outputs" / "positive_control"


def _read(name: str) -> list[dict[str, str]]:
    with (FROZEN / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], *names: str) -> float | None:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return float(value)
    return None


def _score_row(row: dict[str, str]) -> float:
    status = row.get("status", "")
    if not status and row.get("beta") == "32.0":
        status = row.get("n_targets", "")
    return geometry_order_parameter(
        status=status,
        heldout=_number(row, "heldout", "heldout_violation"),
        null_gap=_number(row, "null_gap"),
        truth_error=_number(row, "truth", "truth_order_error"),
    )


def _rank(values: list[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    order = np.argsort(array, kind="stable")
    ranks = np.empty(array.size, dtype=float)
    start = 0
    while start < array.size:
        stop = start + 1
        while stop < array.size and array[order[stop]] == array[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1)
        start = stop
    return ranks


def _spearman(left: list[float], right: list[float]) -> float:
    x = _rank(left)
    y = _rank(right)
    x -= x.mean()
    y -= y.mean()
    denominator = float(np.sqrt(np.dot(x, x) * np.dot(y, y)))
    return float(np.dot(x, y) / denominator) if denominator else float("nan")


def _record(
    rows: list[dict], source: str, condition: str, raw: dict[str, str]
) -> None:
    status = raw.get("status", "")
    if not status and raw.get("beta") == "32.0":
        status = raw.get("n_targets", "")
    rows.append(
        {
            "source": source,
            "condition": condition,
            "seed": _number(raw, "seed"),
            "status": status,
            "heldout": _number(raw, "heldout", "heldout_violation"),
            "null_gap": _number(raw, "null_gap"),
            "truth_error": _number(raw, "truth", "truth_order_error"),
            "G": _score_row(raw),
            "score_kind": "exact",
            "code_version": git_describe(),
        }
    )


def _summarize(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    return {
        "n": int(array.size),
        "min": float(array.min()),
        "median": float(np.median(array)),
        "mean": float(array.mean()),
        "max": float(array.max()),
        "n_ge_0_5": int(np.sum(array >= 0.5)),
    }


def main() -> None:
    rows: list[dict] = []

    p3 = _read("p3_stage_a_calibration.csv")
    for raw in p3:
        _record(rows, "P3", "sprinkled", raw)
        shuffled_heldout = float(raw["heldout"]) + float(raw["null_gap"])
        upper_bound = clipped_gate_margin_score(
            ((0.10 - shuffled_heldout) / 0.10,)
        )
        rows.append(
            {
                "source": "P3",
                "condition": "column_shuffle",
                "seed": float(raw["seed"]),
                "status": "heldout_upper_bound",
                "heldout": shuffled_heldout,
                "null_gap": "",
                "truth_error": "",
                "G": upper_bound,
                "score_kind": "upper_bound",
                "code_version": git_describe(),
            }
        )

    for raw in _read("p5_stage_a_calibration.csv"):
        _record(rows, "P5", "uniform", raw)
    for raw in _read("p5_stage_b_all.csv"):
        _record(rows, "P5", f"beta={float(raw['beta']):g}", raw)
    for raw in _read("p6_stage_b_layered.csv"):
        _record(rows, "P6", "layered_hard_negative", raw)

    p1_by_seed: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for raw in _read("p1_stage_b_epsilon_sweep.csv"):
        status = raw["status"]
        if status == "ok":
            margins = (
                (0.05 - float(raw["heldout_violation"])) / 0.05,
                (0.15 - float(raw["truth_order_error"])) / 0.15,
            )
            score = clipped_gate_margin_score(margins)
        else:
            score = 0.0
        seed = int(float(raw["seed"]))
        epsilon = float(raw["epsilon"])
        p1_by_seed[seed].append((epsilon, score))
        rows.append(
            {
                "source": "P1",
                "condition": "dilution_analogue",
                "seed": seed,
                "epsilon": epsilon,
                "status": status,
                "heldout": raw.get("heldout_violation", ""),
                "null_gap": "",
                "truth_error": raw.get("truth_order_error", ""),
                "G": score,
                "score_kind": "two_gate_analogue",
                "code_version": git_describe(),
            }
        )

    groups: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        if row["source"] != "P1":
            groups[f"{row['source']}:{row['condition']}"].append(float(row["G"]))
    p1_rho = {}
    for seed, pairs in sorted(p1_by_seed.items()):
        pairs.sort()
        p1_rho[str(seed)] = _spearman(
            [epsilon for epsilon, _ in pairs], [score for _, score in pairs]
        )

    group_summary = {name: _summarize(values) for name, values in groups.items()}
    requirements = {
        "p3_sprinkled_all_ge_0_5": group_summary["P3:sprinkled"]["n_ge_0_5"]
        == group_summary["P3:sprinkled"]["n"],
        "p5_uniform_all_ge_0_5": group_summary["P5:uniform"]["n_ge_0_5"]
        == group_summary["P5:uniform"]["n"],
        "column_shuffle_upper_bound_zero": group_summary[
            "P3:column_shuffle"
        ]["max"]
        == 0.0,
        "crystal_all_zero": group_summary["P5:beta=32"]["max"] == 0.0,
        "p1_all_rho_le_minus_0_8": all(value <= -0.8 for value in p1_rho.values()),
    }
    result = {
        "stage": "P7-A-G",
        "code_version": git_describe(),
        "definition": {
            "normalized_margins": [
                "(0.10 - heldout) / 0.10",
                "(null_gap - 0.10) / 0.10",
                "(0.40 - truth_error) / 0.40",
            ],
            "G": "clip(0.5 + min(normalized_margins), 0, 1)",
            "structural_block": 0.0,
            "gate_equivalent_threshold": 0.5,
        },
        "groups": group_summary,
        "p1_dilution_rho_by_seed": p1_rho,
        "p1_dilution_median_rho": float(np.median(list(p1_rho.values()))),
        "requirements": requirements,
        "all_requirements_met": all(requirements.values()),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    write_rows_csv(OUT / "p7_stage_a_geometry_score.csv", rows)
    (OUT / "p7_stage_a_geometry_score_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
