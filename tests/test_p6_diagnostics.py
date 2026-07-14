"""Unit tests for frozen P6b comparison statistics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p6_diagnostics import (  # noqa: E402
    _instrument_margin,
    _normalize_p5_row,
    roc_auc,
    spearman,
)
from pc_common import git_describe, write_rows_csv  # noqa: E402


def test_roc_auc_handles_perfect_reverse_and_ties():
    labels = [0, 0, 1, 1]
    assert roc_auc(labels, [0.0, 1.0, 2.0, 3.0]) == 1.0
    assert roc_auc(labels, [3.0, 2.0, 1.0, 0.0]) == 0.0
    assert roc_auc(labels, [1.0, 1.0, 1.0, 1.0]) == 0.5


def test_spearman_uses_average_tie_ranks():
    assert spearman([0.0, 1.0, 2.0], [2.0, 1.0, 0.0]) == pytest.approx(-1.0)
    assert spearman([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 1.0, 2.0]) > 0.9


def test_p1_instrument_margin_accepts_frozen_column_names():
    row = {
        "status": "ok",
        "heldout_violation": "0.025",
        "truth_order_error": "0.075",
    }
    assert _instrument_margin("P1", row) == pytest.approx(0.5)


def test_legacy_short_p5_crystal_row_is_normalized():
    row = {
        "beta": "32.0",
        "heldout": "1",
        "min_chain_len": "101",
        "n_targets": "structural_block: only 0 chains",
        "sample": None,
        "seed": None,
        "status": None,
    }
    normalized = _normalize_p5_row(row)
    assert normalized["sample"] == "1"
    assert normalized["seed"] == "101"
    assert normalized["status"] == "structural_block: only 0 chains"


def test_shared_csv_writer_uses_lf_line_endings(tmp_path):
    path = tmp_path / "rows.csv"
    write_rows_csv(path, [{"seed": 1.0, "status": "ok"}])
    payload = path.read_bytes()
    assert b"\r\n" not in payload
    assert payload.count(b"\n") == 2


def test_git_describe_is_cached(monkeypatch):
    calls = 0

    def fake_run(*args, **kwargs):
        nonlocal calls
        calls += 1
        return type("Result", (), {"stdout": "abc123\n"})()

    git_describe.cache_clear()
    monkeypatch.setattr("pc_common.subprocess.run", fake_run)
    try:
        assert git_describe() == "abc123"
        assert git_describe() == "abc123"
        assert calls == 1
    finally:
        git_describe.cache_clear()
