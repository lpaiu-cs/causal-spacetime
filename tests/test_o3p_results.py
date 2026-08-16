"""O3' campaign-result contract tests.

The re-certification at target ratio 0.005 was executed ONCE, from the
clean exact checkout of the freeze BRANCH HEAD `30d286a` -- the second
parent of the merge commit, per the execution ruling -- under the
frozen caps. The committed artifact is that single observation; these
tests pin its lineage, its termination provenance, the target
arithmetic, the O3-intersection consistency, and the gated V_ref'
recommendation, and verify the executed freeze surface against the
historical blobs at the freeze head -- so the result stays checkable
long after the CURRENT manifest is re-pinned for later work.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
#: The approved execution SHA: the freeze branch head, NOT the merge.
_FREEZE = "30d286a52790bd169d07a8eab181b4bee5551144"
_ART = _REPO / "docs" / "prereg" / "p14_o3p_volume.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_o3p_executed_freeze_manifest.json")

#: The O3 base interval the consistency check ran against.
_O3_LO, _O3_HI = 56.21273686780051, 57.348018526434714


def _artifact() -> dict:
    return json.loads(_ART.read_text(encoding="utf-8"))


def _blob(rel: str) -> bytes:
    return subprocess.run(
        ["git", "cat-file", "blob", f"{_FREEZE}:{rel}"], cwd=_REPO,
        capture_output=True, check=True).stdout


def test_lineage_is_the_freeze_branch_head_clean_at_entry_and_exit():
    """The approval's scope: one run, from the clean exact checkout of
    the freeze BRANCH HEAD (the merge commit's second parent, never
    the merge commit) -- start and end revs both that commit, clean,
    and the exit rev compared to the approved SHA by the runner."""

    code = _artifact()["code"]
    assert code["start"] == {"rev": _FREEZE, "dirty": False}
    assert code["end"] == {"rev": _FREEZE, "dirty": False}
    merge_parent2 = subprocess.run(
        ["git", "rev-parse",
         "3b869c3a9cb70f6d1d1354c190260c3f8358892f^2"],
        cwd=_REPO, capture_output=True, text=True,
        check=True).stdout.strip()
    assert merge_parent2 == _FREEZE


def test_termination_is_target_met_with_full_provenance():
    r = _artifact()["result"]
    assert r["status"] == "target-met"
    assert r["termination_reason"] == "target-met"
    assert r["max_depth_reached"] == 14
    assert r["cells_at_max_depth"] == 0
    assert r["uncosted_cells"] == 0
    assert r["intersection_active"] is False


def test_the_certified_interval_is_the_published_one():
    """The observation, pinned to the exact serialized endpoints."""

    r = _artifact()["result"]
    assert r["v_lo"] == 56.49295916680885
    assert r["v_hi"] == 57.06066656647874
    assert r["calls"] == 565_478
    assert r["cells"] == 141_488


def test_caps_were_respected_not_raised():
    art = _artifact()
    import sys
    sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
    import o3p_frozen_volume as o3p
    assert art["frozen_config"] == o3p.FROZEN
    r = art["result"]
    assert r["calls"] <= o3p.FROZEN["max_calls"]
    assert r["wall_s"] <= o3p.FROZEN["max_wall_s"]
    assert r["max_depth_reached"] <= o3p.FROZEN["max_depth"]


def test_target_arithmetic_holds_on_the_serialized_endpoints():
    """(V_hi - V_lo)/(V_hi + V_lo) <= 0.005 on the OUTWARD binary64
    endpoints exactly as committed -- outward serialization widens, so
    target-met must survive the widening."""

    r = _artifact()["result"]
    assert 0.0 < r["v_lo"] < r["v_hi"]
    assert (r["v_hi"] - r["v_lo"]) / (r["v_hi"] + r["v_lo"]) <= 0.005


def test_the_o3_intersection_is_consistent_and_contained():
    """The consistency check the freeze added: both certifications
    enclose the same true volume, and this run's interval lies
    strictly inside the O3 interval, so the intersection equals the
    standalone interval and no inconsistency is recorded."""

    art = _artifact()
    x = art["intersection_with_o3"]
    r = art["result"]
    assert x["consistent"] is True
    assert x["certification_inconsistency"] is False
    assert x["o3_base"]["v_lo"] == _O3_LO
    assert x["o3_base"]["v_hi"] == _O3_HI
    assert _O3_LO < r["v_lo"] < r["v_hi"] < _O3_HI
    assert (x["v_lo"], x["v_hi"]) == (r["v_lo"], r["v_hi"])
    assert "STANDALONE" in x["note"]


def test_v_ref_prime_is_the_standalone_midpoint_and_a_recommendation():
    """Emitted because the check was consistent; the value is the
    standalone midpoint, and its adoption (plus every O5 re-sizing)
    belongs to the O5 freeze."""

    art = _artifact()
    v = art["v_ref_prime_recommendation"]
    r = art["result"]
    assert v["value"] == 0.5 * (r["v_lo"] + r["v_hi"])
    assert v["definition"] == "midpoint of the STANDALONE O3' interval"
    assert "recommendation only" in v["status"]


def test_environment_matches_the_executed_lock():
    art = _artifact()
    snap = json.loads(_SNAP.read_text(encoding="utf-8"))
    assert art["environment"] == snap["environment"]


def test_executed_freeze_snapshot_matches_the_historical_blobs():
    """The S4/S5/O3 pattern: the snapshot must be byte-identical to
    the freeze-head blob of the manifest, and every digest it froze
    must equal the blob it pinned at that commit."""

    hist = _blob("docs/prereg/p14_o3p_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    assert len(m["files"]) == 11
    for rel, want in m["files"].items():
        assert hashlib.sha256(_blob(rel)).hexdigest() == want, rel


def test_curve_endpoints_enclose_the_final_interval():
    art = _artifact()
    curve = art["curve"]
    r = art["result"]
    calls = [s["calls"] for s in curve]
    assert calls == sorted(calls)
    last = curve[-1]
    assert last["v_lo"] <= r["v_lo"] + 1e-12
    assert last["v_hi"] >= r["v_hi"] - 1e-12


def test_the_actual_cost_sits_inside_the_projected_range():
    """The freeze's two-curve projection said ~467k-750k calls; the
    run spent 565,478 -- inside the range, 2.5% off the central p=2
    projection. Recorded as a fact about the projection method, not
    as a certification of it."""

    import sys
    sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
    import o3p_projection as proj
    t = proj.table(0.005)
    lo, hi = t["range_calls"]
    assert lo <= _artifact()["result"]["calls"] <= hi
