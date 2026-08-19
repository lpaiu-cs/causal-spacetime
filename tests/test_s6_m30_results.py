"""S6 rung M = 3.0 oracle-result contract tests.

The rung's certified volume was executed ONCE, from the clean exact
checkout of the freeze BRANCH HEAD `f980833` (the second parent of
merge `908f6cc`, PR #106), under the frozen caps. These tests pin the
lineage, the target arithmetic on the serialized endpoints, the cap
discipline, the v_ref_rung recommendation, and the executed freeze
surface against the historical blobs at the freeze head."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_FREEZE = "f9808339c5e411fce1045be9e90a18f25133c7ee"
_MERGE = "908f6cc"
_ART = _REPO / "docs" / "prereg" / "p14_s6_m30_volume.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_s6_m30_executed_freeze_manifest.json")

sys.path.insert(0, str(_REPO / "experiments" / "oracle"))

import s6_m30_frozen_volume as run  # noqa: E402
import s6_rungs as s6  # noqa: E402


def _artifact() -> dict:
    return json.loads(_ART.read_text(encoding="utf-8"))


def _blob(rel: str) -> bytes:
    return subprocess.run(
        ["git", "cat-file", "blob", f"{_FREEZE}:{rel}"], cwd=_REPO,
        capture_output=True, check=True).stdout


def test_lineage_is_the_freeze_branch_head_clean_at_entry_and_exit():
    code = _artifact()["code"]
    assert code["start"] == {"rev": _FREEZE, "dirty": False}
    assert code["end"] == {"rev": _FREEZE, "dirty": False}
    merge_parent2 = subprocess.run(
        ["git", "rev-parse", f"{_MERGE}^2"], cwd=_REPO,
        capture_output=True, text=True, check=True).stdout.strip()
    assert merge_parent2 == _FREEZE


def test_termination_is_target_met_with_full_provenance():
    r = _artifact()["result"]
    assert r["status"] == "target-met"
    assert r["termination_reason"] == "target-met"
    assert r["uncosted_cells"] == 0
    assert r["max_depth_reached"] == 15
    assert r["cells_at_max_depth"] == 0


def test_the_certified_interval_is_the_published_one():
    r = _artifact()["result"]
    assert r["v_lo"] == 121.22502058308072
    assert r["v_hi"] == 122.44328048247328
    assert r["calls"] == 455_520


def test_caps_were_respected_not_raised():
    art = _artifact()
    assert art["frozen_config"] == run.FROZEN
    r = art["result"]
    assert r["calls"] <= run.FROZEN["max_calls"]
    assert r["wall_s"] <= run.FROZEN["max_wall_s"]
    assert r["max_depth_reached"] <= run.FROZEN["max_depth"]


def test_target_arithmetic_holds_on_the_serialized_endpoints():
    r = _artifact()["result"]
    assert 0.0 < r["v_lo"] < r["v_hi"]
    assert (r["v_hi"] - r["v_lo"]) / (r["v_hi"] + r["v_lo"]) <= 0.005


def test_the_rung_identity_is_the_ladder_constant():
    art = _artifact()
    assert art["rung"]["m"] == 3.0
    assert art["rung"]["mu"] == s6.RUNG_CONSTANTS[3.0]["mu"]
    assert art["rung"]["scale"] == s6.RUNG_CONSTANTS[3.0]["scale"]
    assert art["frozen_config"]["dt"] == s6.RUNG_CONSTANTS[3.0]["dt"]
    # the lemma-table snapshot the import gate verified is embedded
    assert all(row["certified"]
               for row in art["lemma_table"].values())


def test_v_ref_rung_is_the_midpoint_and_a_recommendation():
    art = _artifact()
    v = art["v_ref_rung_recommendation"]
    r = art["result"]
    assert v["value"] == 0.5 * (r["v_lo"] + r["v_hi"])
    assert "recommendation only" in v["status"]


def test_environment_matches_the_executed_lock():
    art = _artifact()
    snap = json.loads(_SNAP.read_text(encoding="utf-8"))
    assert art["environment"] == snap["environment"]


def test_executed_freeze_snapshot_matches_the_historical_blobs():
    hist = _blob("docs/prereg/p14_s6_m30_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    assert set(m["files"]) == set(run.PROTOCOL_SURFACE)
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
