"""O4b checkpoints: progress that survives, and never reads as a
result.

O4 lost twelve hours of statistics because nothing intermediate was
written. These tests hold the two properties that make the fix real --
the write is atomic, and the record cannot be mistaken for a verdict.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))

import o4b_checkpoint as cp  # noqa: E402


def _payload(**over) -> dict:
    base = {
        "freeze_sha": "0" * 40,
        "manifest_digest": "a" * 64,
        "seed": 40_000_401,
        "rng_position": 8_192,
        "samples": 8_192,
        "statistics": {"n": 8_192, "mean_z": 0.25, "var_z": 0.03},
        "budget": {"calls": 16_384, "wall_s": 12.5},
    }
    base.update(over)
    return base


def test_a_checkpoint_says_which_freeze_which_stream_and_how_far(
        tmp_path):
    """Resuming or auditing needs all three. A checkpoint that has the
    numbers but cannot say which run produced them is not evidence."""

    path = cp.write(tmp_path / "ck.json", "g1_chunk", _payload())
    record = json.loads(path.read_text(encoding="utf-8"))
    for key in cp.REQUIRED:
        assert key in record
    assert record["rng_position"] == 8_192


def test_a_missing_required_field_fails_at_the_write(tmp_path):
    """Not at the recovery -- by then the run that could have supplied
    it is gone."""

    for key in cp.REQUIRED:
        payload = _payload()
        del payload[key]
        with pytest.raises(ValueError, match=key):
            cp.write(tmp_path / "ck.json", "g1_chunk", payload)


def test_the_stages_are_the_four_the_freeze_names(tmp_path):
    assert cp.STAGES == ("g3b", "g1_chunk", "g1_complete",
                         "g2_complete")
    with pytest.raises(ValueError, match="frozen checkpoint stages"):
        cp.write(tmp_path / "ck.json", "g3a", _payload())


def test_every_checkpoint_is_stamped_partial_and_non_verdict(tmp_path):
    """Written by the module, not passed in, so a caller can neither
    forget the stamp nor claim otherwise."""

    for stage in cp.STAGES:
        path = cp.write(tmp_path / f"{stage}.json", stage, _payload())
        record = json.loads(path.read_text(encoding="utf-8"))
        assert record["partial"] is True
        assert record["non_verdict"] is True
        assert record["kind"] == "checkpoint"
        assert cp.is_verdict(record) is False


@pytest.mark.parametrize("key", cp.RESERVED)
def test_a_payload_may_not_supply_the_keys_this_module_writes(
        tmp_path, key):
    """Refused, not overridden. Ordering the spread so the stamps win
    would protect `partial`, but a payload carrying `stage` is a
    caller confusion -- a reused run-state dict whose `stage` is the
    run stage `g1`, not the checkpoint point `g1_chunk`. Winning hides
    it; losing writes a `stage` outside STAGES and a resume continues
    from the wrong place."""

    with pytest.raises(ValueError, match="which this module writes"):
        cp.write(tmp_path / "ck.json", "g1_complete",
                 _payload(**{key: "g1"}))
    assert not (tmp_path / "ck.json").exists()


def test_a_caller_cannot_stamp_a_checkpoint_as_final(tmp_path):
    """The specific case the stamps exist for: a payload claiming the
    record is a result does not get to write one."""

    with pytest.raises(ValueError, match="reading as a result"):
        cp.write(tmp_path / "ck.json", "g1_complete",
                 _payload(partial=False, non_verdict=False))

    record = json.loads(cp.write(tmp_path / "ck.json", "g1_complete",
                                 _payload()).read_text(
                                     encoding="utf-8"))
    assert record["partial"] is True
    assert record["non_verdict"] is True
    assert record["stage"] == "g1_complete"


def test_a_failed_write_leaves_the_previous_checkpoint_intact(
        tmp_path, monkeypatch):
    """The whole point. A run dies at an awkward moment; the last good
    checkpoint has to still be there."""

    path = tmp_path / "ck.json"
    cp.write(path, "g1_chunk", _payload(samples=1_000))
    good = json.loads(path.read_text(encoding="utf-8"))

    def exploding(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(cp.os, "replace", exploding)
    with pytest.raises(OSError):
        cp.write(path, "g1_chunk", _payload(samples=2_000))

    assert json.loads(path.read_text(encoding="utf-8")) == good
    # and nothing was left behind to be mistaken for a checkpoint
    assert [p.name for p in tmp_path.iterdir()] == ["ck.json"]


def test_a_checkpoint_is_meant_to_be_overwritten(tmp_path):
    """Unlike the incident artifacts, which are write-once. Every G1
    chunk replaces the last -- an accumulating pile would be a second
    thing to reconcile."""

    path = tmp_path / "ck.json"
    cp.write(path, "g1_chunk", _payload(samples=1_000))
    cp.write(path, "g1_chunk", _payload(samples=2_000))
    assert cp.read(path)["samples"] == 2_000


def test_no_checkpoint_yet_is_not_an_error(tmp_path):
    """The run may have died before the first chunk, and that has to
    be distinguishable from a corrupt file, which raises."""

    assert cp.read(tmp_path / "absent.json") == {}
    (tmp_path / "broken.json").write_text("{not json",
                                          encoding="utf-8")
    with pytest.raises(json.JSONDecodeError):
        cp.read(tmp_path / "broken.json")


def test_the_document_freezes_the_four_points_and_their_contents():
    doc = (_REPO / "docs" / "prereg"
           / "p14_o4_g3_prereg_reopen.md").read_text(encoding="utf-8")
    assert "### 5.2 원자적 체크포인트" in doc
    assert "G3b 종료 · 각 G1 청크 · G1 완료 · G2 완료" in doc
    assert "RNG 위치" in doc
    assert "os.replace" in doc
    assert "`partial` · `non-verdict`" in doc
