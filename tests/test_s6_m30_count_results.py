"""S6 M=3.0 count-campaign result contract tests.

Executed ONCE from freeze head `ed353d3` (second parent of merge
`f9d6067`, PR #112). CONCORDANT: K_certain = 53,606 (inside the
frozen acceptance [52,986, 54,220]) with ZERO ambiguous points,
D = [-1.017, +0.987] contained in the tau band. These tests pin the
lineage, the accounting, the bit-exact decision reproduction, the
reservation provenance, the retired seed, and the executed
surface."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_FREEZE = "ed353d3fd18c8f6cfd12989b14bf55ef9f7186fd"
_MERGE = "f9d6067"
_ART = _REPO / "docs" / "prereg" / "p14_s6_m30_count.json"
_CK = _REPO / "docs" / "prereg" / "p14_s6_m30_count_checkpoint.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_s6_m30_count_executed_freeze_manifest.json")

sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import probe_seed_ledger as ledger  # noqa: E402
import s6_m30_count as camp  # noqa: E402


def _artifact() -> dict:
    return json.loads(_ART.read_text(encoding="utf-8"))


def _blob(rel: str) -> bytes:
    return subprocess.run(
        ["git", "cat-file", "blob", f"{_FREEZE}:{rel}"], cwd=_REPO,
        capture_output=True, check=True).stdout


def test_lineage_is_the_freeze_branch_head_clean_at_entry_and_exit():
    code = _artifact()["code"]
    assert code["start"]["rev"] == _FREEZE
    assert code["end"]["rev"] == _FREEZE
    assert not code["start"]["dirty"] and not code["end"]["dirty"]
    merge_parent2 = subprocess.run(
        ["git", "rev-parse", f"{_MERGE}^2"], cwd=_REPO,
        capture_output=True, text=True, check=True).stdout.strip()
    assert merge_parent2 == _FREEZE
    assert _artifact()["freeze_sha"] == _FREEZE


def test_the_scan_accounting_is_consistent():
    art = _artifact()
    scan = art["scan"]
    assert scan["n_total"] == 4_605_671
    assert scan["points_done"] == scan["n_total"]
    assert scan["k_certain"] == 53_606
    assert scan["u_ambiguous"] == 0
    assert scan["calls"] == 6_083_193
    assert scan["n_total"] <= scan["calls"] <= 2 * scan["n_total"]
    b = art["budget"]
    assert b["reserved"] == scan["calls"] + camp.G3A_PREFLIGHT_CALLS
    assert b["reserved"] <= camp.FROZEN["max_calls"]
    assert art["frozen_config"] == camp.FROZEN


def test_the_published_decision_is_the_frozen_rule_bit_exactly():
    art = _artifact()
    assert art["decision"] == camp.decide(53_606, 0)
    assert art["decision"]["verdict"] == "CONCORDANT"
    k_lo, k_hi = camp.ACCEPTANCE
    assert k_lo <= 53_606 and 53_606 <= k_hi


def test_containment_holds_on_the_serialized_values():
    d = _artifact()["decision"]
    assert d["d_lo"] >= -d["band"]
    assert d["d_hi"] <= d["band"]
    assert d["c_lo_volume"] < d["c_hi_volume"]


def test_the_reservation_is_recorded_and_the_ref_retained():
    res = _artifact()["reservation"]
    assert res["ref"] == "refs/s6m30/reservation"
    assert res["object"] == "05ad01058677c96d57433c15c543c3010900c48c"
    assert res["seed_spent"] is True
    assert res["verified_at_exit"] is True


def test_the_seed_is_the_retired_stream_and_only_that():
    art = _artifact()
    assert art["seed"] == {"s6_m30_count": 40_000_511}
    assert ledger.OBSERVED_PROBE_SCALARS["s6_m30_count"] == 40_000_511
    assert "s6_m30_count" not in ledger.FRESH_PROBE_SCALARS


def test_the_final_checkpoint_is_a_non_verdict_partial_at_full_n():
    ck = json.loads(_CK.read_text(encoding="utf-8"))
    art = _artifact()
    assert ck["partial"] is True and ck["non_verdict"] is True
    assert ck["points_done"] == art["scan"]["points_done"]
    assert ck["k_certain"] == art["scan"]["k_certain"]
    assert ck["u_ambiguous"] == art["scan"]["u_ambiguous"]
    assert ck["rng_stream"] == "s6_m30_count"
    assert ck["freeze_sha"] == _FREEZE
    assert ck["manifest_digest"] == art["manifest_digest"]


def test_the_artifact_names_the_executed_manifest_digest():
    art = _artifact()
    snap_digest = hashlib.sha256(
        _SNAP.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    assert art["manifest_digest"] == snap_digest
    snap = json.loads(_SNAP.read_text(encoding="utf-8"))
    assert art["environment"] == snap["environment"]


def test_executed_freeze_snapshot_matches_the_historical_blobs():
    hist = _blob("docs/prereg/p14_s6_m30_count_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    for rel, want in m["files"].items():
        assert hashlib.sha256(_blob(rel)).hexdigest() == want, rel


def test_the_freeze_commit_carried_the_count_scalar_as_fresh():
    blob = subprocess.run(
        ["git", "show", f"{_FREEZE}:experiments/"
         "positive_control/probe_seed_ledger.py"],
        cwd=_REPO, capture_output=True, text=True, check=True).stdout
    frozen = blob.split("FRESH_PROBE_SCALARS")[1].split("}")[0]
    assert '"s6_m30_count": 40_000_511' in frozen
