"""Focused tests for P7 N=600 summary statistics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p7_n600_stage_a import (  # noqa: E402
    _binder,
    _group_chains_by_start,
    _seed,
    _validate_shard_contents,
)


def test_seed_formula_separates_beta_and_chain():
    assert _seed(900000, 14.0, 0) == 914000
    assert _seed(900000, 14.0, 1) == 914001
    assert _seed(900000, 18.0, 0) != _seed(900000, 14.0, 0)


def test_centered_binder_is_finite_for_nonconstant_series():
    value = _binder([1.0, 2.0, 3.0, 4.0])
    assert value == pytest.approx(0.45333333333333337)


def test_start_grouping_does_not_depend_on_constant_order_or_chain_ids():
    chains = [
        {"chain": 9, "start": "bipartite", "phase": "intermediate"},
        {"chain": 7, "start": "random", "phase": "continuum"},
        {"chain": 3, "start": "random", "phase": "continuum"},
    ]
    grouped = _group_chains_by_start(chains)
    assert [row["chain"] for row in grouped["random"]] == [3, 7]
    assert [row["chain"] for row in grouped["bipartite"]] == [9]


def test_start_grouping_rejects_wrong_multiplicity():
    with pytest.raises(SystemExit, match="exactly two random"):
        _group_chains_by_start(
            [
                {"chain": 0, "start": "random"},
                {"chain": 2, "start": "bipartite"},
            ]
        )


def _shard_rows(count: int, *, instrument: bool = False):
    indices = [0, count // 2, count - 1] if instrument else range(count)
    return [
        {"beta": "14", "chain": "0", "start": "random", "sample": str(index)}
        for index in indices
    ]


def test_shard_validation_accepts_complete_chain_and_instruments():
    _validate_shard_contents(
        beta=14.0,
        chain=0,
        start="random",
        chain_rows=_shard_rows(45),
        instrument_rows=_shard_rows(45, instrument=True),
        minimum_samples=45,
        nominal_samples=48,
    )


def test_shard_validation_rejects_truncated_instruments():
    with pytest.raises(SystemExit, match="instrument.*sample indices"):
        _validate_shard_contents(
            beta=14.0,
            chain=0,
            start="random",
            chain_rows=_shard_rows(45),
            instrument_rows=_shard_rows(45, instrument=True)[:-1],
            minimum_samples=45,
            nominal_samples=48,
        )


def test_shard_validation_rejects_wrong_identity():
    rows = _shard_rows(45)
    rows[-1]["beta"] = "16"
    with pytest.raises(SystemExit, match="identity mismatch"):
        _validate_shard_contents(
            beta=14.0,
            chain=0,
            start="random",
            chain_rows=rows,
            instrument_rows=_shard_rows(45, instrument=True),
            minimum_samples=45,
            nominal_samples=48,
        )
