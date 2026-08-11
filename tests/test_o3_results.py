"""O3 campaign-result contract tests.

The frozen (12, 18, 8.5)M volume was executed ONCE, from the clean
exact checkout of the freeze merge commit, under the frozen caps. The
committed artifact is that single observation; these tests pin its
lineage, its termination provenance, and the target arithmetic to the
ruling, and verify the executed freeze surface against the historical
blobs at the freeze commit -- so the result stays checkable long after
the CURRENT manifest is re-pinned for later work.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_FREEZE = "785148ecf8be8b7b1baaa1f3866d5d827b8dfdf7"
_ART = _REPO / "docs" / "prereg" / "p14_o3_volume.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_o3_executed_freeze_manifest.json")


def _artifact() -> dict:
    return json.loads(_ART.read_text(encoding="utf-8"))


def _blob(rel: str) -> bytes:
    return subprocess.run(
        ["git", "cat-file", "blob", f"{_FREEZE}:{rel}"], cwd=_REPO,
        capture_output=True, check=True).stdout


def test_lineage_is_the_freeze_commit_clean_at_entry_and_exit():
    """The approval's scope: one run, from the clean exact checkout
    of the freeze merge -- start and end revs must both be that
    commit and both clean."""

    code = _artifact()["code"]
    assert code["start"] == {"rev": _FREEZE, "dirty": False}
    assert code["end"] == {"rev": _FREEZE, "dirty": False}


def test_termination_is_target_met_with_full_provenance():
    """The ruling's provenance fields, with the values the run
    produced: target-met, depth headroom visible (12 of 18), no
    uncosted cells, and the L6d intersection never active -- which is
    what makes the recorded mode decomposition exact rather than an
    upper bound."""

    r = _artifact()["result"]
    assert r["status"] == "target-met"
    assert r["termination_reason"] == "target-met"
    assert r["max_depth_reached"] == 12
    assert r["cells_at_max_depth"] == 0
    assert r["uncosted_cells"] == 0
    assert r["intersection_active"] is False


def test_caps_were_respected_not_raised():
    """The frozen caps bound the run they governed: the artifact's
    frozen_config must equal the ruled freeze, and the actual spend
    must sit inside every cap."""

    art = _artifact()
    import sys
    sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
    import o3_frozen_volume as o3
    assert art["frozen_config"] == o3.FROZEN
    r = art["result"]
    assert r["calls"] <= o3.FROZEN["max_calls"]
    assert art["result"]["wall_s"] <= o3.FROZEN["max_wall_s"]
    assert r["max_depth_reached"] <= o3.FROZEN["max_depth"]


def test_target_arithmetic_holds_on_the_serialized_endpoints():
    """The frozen target on the OUTWARD binary64 endpoints exactly as
    committed: V_lo > 0 (the ratio's precondition) and
    (V_hi - V_lo)/(V_hi + V_lo) <= 0.01 -- outward serialization
    widens the interval, so the target must survive the widening for
    the artifact to state target-met."""

    r = _artifact()["result"]
    assert 0.0 < r["v_lo"] < r["v_hi"]
    assert (r["v_hi"] - r["v_lo"]) / (r["v_hi"] + r["v_lo"]) <= 0.01


def test_environment_matches_the_executed_lock():
    """The instrument that ran is the instrument that was frozen:
    the artifact's recorded environment must equal the executed
    manifest's lock, MPFR included."""

    art = _artifact()
    snap = json.loads(_SNAP.read_text(encoding="utf-8"))
    assert art["environment"] == snap["environment"]


def test_executed_freeze_snapshot_matches_the_historical_blobs():
    """The S4/S5 pattern: the snapshot must be byte-identical to the
    freeze-commit blob of the manifest, and every digest it froze
    must equal the blob it pinned at that commit -- so the executed
    protocol surface stays verifiable even after the CURRENT manifest
    is re-pinned for later work."""

    hist = _blob("docs/prereg/p14_o3_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    assert len(m["files"]) == 8
    for rel, want in m["files"].items():
        assert hashlib.sha256(_blob(rel)).hexdigest() == want, rel


def test_curve_endpoints_enclose_the_final_interval():
    """The convergence curve's last sample was serialized with the
    same outward discipline as the final result: its enclosure must
    contain the final certified interval, and the curve must be
    monotone in calls."""

    art = _artifact()
    curve = art["curve"]
    r = art["result"]
    calls = [s["calls"] for s in curve]
    assert calls == sorted(calls)
    last = curve[-1]
    assert last["v_lo"] <= r["v_lo"] + 1e-12
    assert last["v_hi"] >= r["v_hi"] - 1e-12
