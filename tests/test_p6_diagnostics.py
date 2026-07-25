"""Unit tests for frozen P6b comparison statistics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p1_epsilon_sweep import _p1_gate_pass  # noqa: E402
from p6_diagnostics import (  # noqa: E402
    _instrument_margin,
    _normalize_p5_row,
    _p5_row_key,
    _validate_expected_rows,
    roc_auc,
    spearman,
)
from pc_common import git_describe, write_rows_csv  # noqa: E402

from causal_spacetime_lab.positive_control.p6b_margin import (  # noqa: E402
    p6b_instrument_margin,
)


def test_roc_auc_handles_perfect_reverse_and_ties():
    labels = [0, 0, 1, 1]
    assert roc_auc(labels, [0.0, 1.0, 2.0, 3.0]) == 1.0
    assert roc_auc(labels, [3.0, 2.0, 1.0, 0.0]) == 0.0
    assert roc_auc(labels, [1.0, 1.0, 1.0, 1.0]) == 0.5


def test_spearman_uses_average_tie_ranks():
    assert spearman([0.0, 1.0, 2.0], [2.0, 1.0, 0.0]) == pytest.approx(-1.0)
    assert spearman([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 1.0, 2.0]) > 0.9


def test_historical_p1_proxy_remains_a_faithful_frozen_replayer():
    row = {
        "status": "ok",
        "heldout_violation": "0.025",
        "truth_order_error": "0.075",
        "restart_order_disagreement": "0.12",
    }
    assert _instrument_margin("P1", row) == pytest.approx(0.5)
    assert p6b_instrument_margin("P1", row) == pytest.approx(0.2)


def test_p1_endpoint_gate_includes_restart_stability():
    common = {
        "truth": 0.05,
        "heldout": 0.01,
        "truth_max": 0.15,
        "heldout_max": 0.05,
        "stability_max": 0.15,
    }
    assert _p1_gate_pass(stability=0.15, **common)
    assert not _p1_gate_pass(stability=0.151, **common)


def test_p1_order_only_margin_retains_restart_stability():
    row = {
        "status": "ok",
        "heldout_violation": "0.01",
        "truth_order_error": "0.75",
        "restart_order_disagreement": "0.12",
    }
    assert p6b_instrument_margin(
        "P1", row, include_truth=False
    ) == pytest.approx(0.2)


def test_p1_margin_rejects_missing_declared_stability_gate():
    row = {
        "status": "ok",
        "heldout_violation": "0.025",
        "truth_order_error": "0.075",
    }
    with pytest.raises(ValueError, match="restart-stability"):
        p6b_instrument_margin("P1", row)


def test_p5_full_margin_requires_truth_but_order_only_does_not():
    row = {
        "status": "ok",
        "heldout": "0.05",
        "null_gap": "0.20",
    }
    with pytest.raises(ValueError, match="P5 row has no truth-order"):
        p6b_instrument_margin("P5", row)
    assert p6b_instrument_margin(
        "P5", row, include_truth=False
    ) == pytest.approx(0.5)


def test_p3_has_no_truth_gate_in_either_variant():
    row = {
        "status": "ok",
        "heldout": "0.05",
        "null_gap": "0.20",
    }
    assert p6b_instrument_margin("P3", row) == pytest.approx(
        p6b_instrument_margin("P3", row, include_truth=False)
    )


def test_valid_p3_row_requires_its_declared_null_gap():
    row = {"status": "ok", "heldout": "0.05"}
    with pytest.raises(ValueError, match="P3 row has no null-gap"):
        p6b_instrument_margin("P3", row)


def test_p6b_margin_rejects_unknown_source():
    row = {"status": "ok", "heldout": "0.01"}
    with pytest.raises(ValueError, match="unsupported P6b source"):
        p6b_instrument_margin("unknown", row)


def test_p6b_margin_validates_source_before_structural_block():
    row = {"status": "structural_block: no chains"}
    with pytest.raises(ValueError, match="unsupported P6b source"):
        p6b_instrument_margin("unknown", row)
    assert p6b_instrument_margin("P3", row) == -5.0


def test_p6b_margin_rejects_missing_status():
    with pytest.raises(ValueError, match="P3 row has no valid status"):
        p6b_instrument_margin("P3", {"heldout": "0.01", "null_gap": "0.20"})


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


def _p5_rows(samples=(0, 1, 2)):
    return [
        {
            "source": "P5",
            "beta": "2",
            "seed": "100",
            "sample": str(sample),
        }
        for sample in samples
    ]


def test_p6b_row_validation_accepts_exact_frozen_keys():
    expected = {("P5", 2.0, 100, sample) for sample in (0, 1, 2)}
    _validate_expected_rows("P5", _p5_rows(), expected, _p5_row_key)


def test_p6b_row_validation_rejects_truncated_shard():
    expected = {("P5", 2.0, 100, sample) for sample in (0, 1, 2)}
    with pytest.raises(SystemExit, match="missing=.*2"):
        _validate_expected_rows("P5", _p5_rows((0, 1)), expected, _p5_row_key)


def test_p6b_row_validation_rejects_duplicate_key():
    expected = {("P5", 2.0, 100, sample) for sample in (0, 1, 2)}
    rows = _p5_rows()
    rows[-1] = dict(rows[0])
    with pytest.raises(SystemExit, match="duplicates="):
        _validate_expected_rows("P5", rows, expected, _p5_row_key)
