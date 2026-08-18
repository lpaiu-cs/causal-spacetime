"""S6 M=1.8 ambiguity-pilot result contract tests.

Executed ONCE from freeze head `994ec89` (second parent of merge
`8606bde`, PR #99). FEASIBLE at k = 0 of 3,940,846 -- the deepest
rung produced ZERO ambiguous events, calls exactly 2n. These tests
pin the lineage, the call accounting, the bit-exact decision
reproduction, the reservation provenance, the retired seed, and the
executed surface."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_FREEZE = "994ec895ad8a6a93890da7b095c20d5aec094300"
_MERGE = "8606bde"
_ART = _REPO / "docs" / "prereg" / "p14_s6_m18_pilot.json"
_CK = _REPO / "docs" / "prereg" / "p14_s6_m18_pilot_checkpoint.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_s6_m18_pilot_executed_freeze_manifest.json")

sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import probe_seed_ledger as ledger  # noqa: E402
import s6_m18_amb_pilot as pilot  # noqa: E402


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


def test_the_scan_accounting_is_exactly_two_n():
    """k = 0: every event decided both legs, calls exactly 2n and the
    budget closed exactly at max_calls."""

    art = _artifact()
    scan = art["scan"]
    n = pilot.FROZEN["n_events"]
    assert scan["events_done"] == n == 3_940_846
    assert scan["k_ambiguous"] == 0
    assert scan["calls"] == 2 * n == 7_881_692
    b = art["budget"]
    assert b["reserved"] == scan["calls"] + pilot.G3A_PREFLIGHT_CALLS
    assert b["reserved"] == pilot.FROZEN["max_calls"]
    assert art["frozen_config"] == pilot.FROZEN


def test_the_published_decision_is_the_frozen_rule_bit_exactly():
    """The boundary case the sizing was built for: k = 1 FEASIBLE by
    the pre-frozen margin, reproduced bit-for-bit by decide(1)."""

    art = _artifact()
    assert art["decision"] == pilot.decide(0)
    assert art["decision"]["verdict"] == "FEASIBLE"
    tail = art["decision"]["tail_p_u_gt_u_max"]
    assert tail <= pilot.FROZEN["tail_budget"]


def test_the_reservation_is_recorded_and_the_ref_retained():
    res = _artifact()["reservation"]
    assert res["ref"] == "refs/s6m18pilot/reservation"
    assert res["object"] == "f695f96c03e5a0781428f4a2ddb77097ff16e94d"
    assert res["seed_spent"] is True
    assert res["verified_at_exit"] is True


def test_the_seed_is_the_retired_stream_and_only_that():
    art = _artifact()
    assert art["seed"] == {"s6_m18_pilot": 40_000_481}
    assert ledger.OBSERVED_PROBE_SCALARS["s6_m18_pilot"] == 40_000_481
    assert "s6_m18_pilot" not in ledger.FRESH_PROBE_SCALARS


def test_the_final_checkpoint_is_a_non_verdict_partial_at_full_n():
    ck = json.loads(_CK.read_text(encoding="utf-8"))
    art = _artifact()
    assert ck["partial"] is True and ck["non_verdict"] is True
    assert ck["events_done"] == art["scan"]["events_done"]
    assert ck["k_ambiguous"] == art["scan"]["k_ambiguous"]
    assert ck["calls"] == art["scan"]["calls"]
    assert ck["rng_stream"] == "s6_m18_pilot"
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
    hist = _blob("docs/prereg/p14_s6_m18_pilot_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    for rel, want in m["files"].items():
        assert hashlib.sha256(_blob(rel)).hexdigest() == want, rel


def test_the_freeze_commit_carried_the_pilot_scalar_as_fresh():
    blob = subprocess.run(
        ["git", "show", f"{_FREEZE}:experiments/"
         "positive_control/probe_seed_ledger.py"],
        cwd=_REPO, capture_output=True, text=True, check=True).stdout
    frozen = blob.split("FRESH_PROBE_SCALARS")[1].split("}")[0]
    assert '"s6_m18_pilot": 40_000_481' in frozen
