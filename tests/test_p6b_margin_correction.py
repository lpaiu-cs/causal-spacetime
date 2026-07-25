"""Regression tests for the analysis-only correction to the frozen P6b margin."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

import pytest

from causal_spacetime_lab.positive_control.p6b_margin import p6b_instrument_margin

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "docs" / "prereg" / "frozen"
CORRECTION = (
    ROOT
    / "docs"
    / "paper"
    / "paper_b"
    / "figures"
    / "p6b_margin_correction.json"
)


def _rows(name: str) -> list[dict[str, str]]:
    with (FROZEN / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _rank_auc(labels: list[int], scores: list[float]) -> float:
    order = sorted(range(len(scores)), key=scores.__getitem__)
    ranks = [0.0] * len(scores)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and scores[order[stop]] == scores[order[start]]:
            stop += 1
        average_rank = 0.5 * (start + 1 + stop)
        for index in order[start:stop]:
            ranks[index] = average_rank
        start = stop
    positive_ranks = sum(
        rank for rank, label in zip(ranks, labels, strict=True) if label == 1
    )
    n_positive = sum(labels)
    n_negative = len(labels) - n_positive
    return (
        positive_ranks - n_positive * (n_positive + 1) / 2
    ) / (n_positive * n_negative)


def _all_corrected_rows() -> list[dict[str, str]]:
    stability = {
        (int(float(row["seed"])), float(row["epsilon"])): row[
            "restart_order_disagreement"
        ]
        for row in _rows("p1_stage_b_epsilon_sweep.csv")
    }
    rows = _rows("p6b_scored_rows.csv")
    for row in rows:
        if row["source"] == "P1":
            row["restart_order_disagreement"] = stability[
                (int(float(row["seed"])), float(row["epsilon"]))
            ]
    return rows


def _corrected_rows() -> list[dict[str, str]]:
    return [row for row in _all_corrected_rows() if row["label"] != ""]


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        average_rank = 0.5 * (start + 1 + stop)
        for index in order[start:stop]:
            ranks[index] = average_rank
        start = stop
    return ranks


def _spearman(xs: list[float], ys: list[float]) -> float:
    rx = _average_ranks(xs)
    ry = _average_ranks(ys)
    mx = statistics.mean(rx)
    my = statistics.mean(ry)
    covariance = sum((x - mx) * (y - my) for x, y in zip(rx, ry, strict=True))
    scale_x = sum((x - mx) ** 2 for x in rx)
    scale_y = sum((y - my) ** 2 for y in ry)
    return covariance / (scale_x * scale_y) ** 0.5


def _historical_incomplete_order_only(row: dict[str, str]) -> float:
    if row["status"] != "ok":
        return -5.0
    heldout = float(row["heldout"])
    heldout_max = 0.05 if row["source"] == "P1" else 0.10
    margins = [(heldout_max - heldout) / heldout_max]
    if row["source"] != "P1":
        margins.append((float(row["null_gap"]) - 0.10) / 0.10)
    return min(margins)


def test_corrected_frozen_p6b_auc_values():
    rows = _corrected_rows()
    labels = [int(float(row["label"])) for row in rows]
    full = [p6b_instrument_margin(row["source"], row) for row in rows]
    order_only = [
        p6b_instrument_margin(row["source"], row, include_truth=False)
        for row in rows
    ]

    assert _rank_auc(labels, full) == pytest.approx(0.9897629310344828)
    assert _rank_auc(labels, order_only) == pytest.approx(0.9679559891107078)


def test_correction_audit_matches_every_cited_derived_value():
    audit = json.loads(CORRECTION.read_text(encoding="utf-8"))
    rows = _corrected_rows()
    labels = [int(float(row["label"])) for row in rows]
    full = [p6b_instrument_margin(row["source"], row) for row in rows]
    order_only = [
        p6b_instrument_margin(row["source"], row, include_truth=False)
        for row in rows
    ]
    score_series = {
        "historical_frozen_proxy": [
            float(row["instrument_margin"]) for row in rows
        ],
        "gate_complete_full": full,
        "historical_incomplete_order_only": [
            _historical_incomplete_order_only(row) for row in rows
        ],
        "gate_complete_order_only": order_only,
        "height_distance": [-float(row["height_distance"]) for row in rows],
        "mm_distance": [-float(row["mm_distance"]) for row in rows],
        "abundance_distance": [-float(row["abundance_distance"]) for row in rows],
    }
    for name, scores in score_series.items():
        assert _rank_auc(labels, scores) == pytest.approx(audit["roc_auc"][name])

    all_p1 = [
        row
        for row in _all_corrected_rows()
        if row["source"] == "P1" and row["status"] == "ok"
    ]
    full_rho: list[float] = []
    order_rho: list[float] = []
    for seed in sorted({int(float(row["seed"])) for row in all_p1}):
        subset = sorted(
            (row for row in all_p1 if int(float(row["seed"])) == seed),
            key=lambda row: float(row["epsilon"]),
        )
        epsilon = [float(row["epsilon"]) for row in subset]
        full_rho.append(
            _spearman(
                epsilon,
                [-p6b_instrument_margin("P1", row) for row in subset],
            )
        )
        order_rho.append(
            _spearman(
                epsilon,
                [
                    -p6b_instrument_margin("P1", row, include_truth=False)
                    for row in subset
                ],
            )
        )
    p1_audit = audit["p1_rank_and_h_lag"]
    assert statistics.median(full_rho) == pytest.approx(
        p1_audit["gate_complete_full_median_rho_epsilon"]
    )
    assert statistics.median(order_rho) == pytest.approx(
        p1_audit["gate_complete_order_only_median_rho_epsilon"]
    )
    h_lag = [
        row
        for row in all_p1
        if float(row["heldout"]) <= 0.05 and float(row["truth"]) > 0.15
    ]
    assert len(h_lag) == p1_audit["h_lag_count"]
    assert (
        sum(
            p6b_instrument_margin("P1", row, include_truth=False) >= 0.0
            for row in h_lag
        )
        == p1_audit["gate_complete_order_only_pass_count"]
    )

    def pair_auc(
        positives: list[dict[str, str]],
        negatives: list[dict[str, str]],
        *,
        order_score: bool,
    ) -> float:
        pair = [*positives, *negatives]
        pair_labels = [1] * len(positives) + [0] * len(negatives)
        if order_score:
            scores = [
                p6b_instrument_margin(
                    row["source"], row, include_truth=False
                )
                for row in pair
            ]
        else:
            scores = [-float(row["height_distance"]) for row in pair]
        return _rank_auc(pair_labels, scores)

    p1_positive = [
        row for row in rows if row["source"] == "P1" and row["label"] == "1.0"
    ]
    p1_negative = [
        row for row in rows if row["source"] == "P1" and row["label"] == "0.0"
    ]
    p5_positive = [
        row for row in rows if row["source"] == "P5" and row["label"] == "1.0"
    ]
    p6_negative = [row for row in rows if row["source"] == "P6"]
    pairs = {
        "p1_positive_vs_p1_negative": (p1_positive, p1_negative),
        "p1_positive_vs_p6_negative": (p1_positive, p6_negative),
        "p5_positive_vs_p6_negative": (p5_positive, p6_negative),
    }
    for name, (positives, negatives) in pairs.items():
        expected = audit["selected_family_pair_auc"][name]
        assert pair_auc(positives, negatives, order_score=True) == pytest.approx(
            expected["gate_complete_order_only"]
        )
        assert pair_auc(positives, negatives, order_score=False) == pytest.approx(
            expected["height_distance"]
        )

    for group in ("p1_scene", "orders_n600"):
        subset = [row for row in rows if row["reference_group"] == group]
        group_labels = [int(float(row["label"])) for row in subset]
        expected = audit["within_reference_group_auc"][group]
        assert len(subset) == expected["n"]
        assert _rank_auc(
            group_labels,
            [p6b_instrument_margin(row["source"], row) for row in subset],
        ) == pytest.approx(expected["gate_complete_full"])
        assert _rank_auc(
            group_labels,
            [
                p6b_instrument_margin(
                    row["source"], row, include_truth=False
                )
                for row in subset
            ],
        ) == pytest.approx(expected["gate_complete_order_only"])
        assert _rank_auc(
            group_labels,
            [-float(row["height_distance"]) for row in subset],
        ) == pytest.approx(expected["height_distance"])
    p3_rows = [row for row in rows if row["reference_group"] == "p3_n1500"]
    assert len(p3_rows) == audit["within_reference_group_auc"]["p3_n1500"]["n"]
    assert {row["label"] for row in p3_rows} == {"0.0"}


def test_correction_changes_ranks_but_not_frozen_pass_block_signs():
    rows = _corrected_rows()
    corrected = [p6b_instrument_margin(row["source"], row) for row in rows]
    historical = [float(row["instrument_margin"]) for row in rows]

    assert corrected != historical
    assert [value >= 0.0 for value in corrected] == [
        value >= 0.0 for value in historical
    ]


def test_corrected_p1_margin_sign_matches_all_three_declared_gates():
    rows = [row for row in _corrected_rows() if row["source"] == "P1"]
    for row in rows:
        margin_pass = p6b_instrument_margin("P1", row) >= 0.0
        declared_pass = (
            row["status"] == "ok"
            and float(row["heldout"]) <= 0.05
            and float(row["truth"]) <= 0.15
            and float(row["restart_order_disagreement"]) <= 0.15
        )
        assert margin_pass == declared_pass
