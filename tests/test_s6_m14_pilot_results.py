"""S6 M=1.4 ambiguity-pilot result contract tests.

Executed ONCE from freeze head `69c79ec` (second parent of merge
`0fcb2c6`, PR #95). FEASIBLE at k = 1 of 6,024,777 -- the pre-sized
boundary case: the deep rung DID produce an ambiguous event, and the
frozen n absorbed it by the designed +3.97e-9 tail margin. These
tests pin the lineage, the call accounting (2n minus the one leg-1
early exit), the bit-exact decision reproduction, the reservation
provenance, the retired seed, and the executed surface."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_FREEZE = "69c79ec29539135d41cab769cd69a4f58bb3041c"
_MERGE = "0fcb2c6"
_ART = _REPO / "docs" / "prereg" / "p14_s6_m14_pilot.json"
_CK = _REPO / "docs" / "prereg" / "p14_s6_m14_pilot_checkpoint.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_s6_m14_pilot_executed_freeze_manifest.json")

sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import probe_seed_ledger as ledger  # noqa: E402
import s6_m14_amb_pilot as pilot  # noqa: E402


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


def test_the_scan_accounting_names_the_one_early_exit():
    """k = 1 with exactly 2n - 1 calls: the single ambiguous event
    exited on leg 1, every other event decided both legs."""

    art = _artifact()
    scan = art["scan"]
    n = pilot.FROZEN["n_events"]
    assert scan["events_done"] == n == 6_024_777
    assert scan["k_ambiguous"] == 1
    assert scan["calls"] == 2 * n - 1 == 12_049_553
    b = art["budget"]
    assert b["reserved"] == scan["calls"] + pilot.G3A_PREFLIGHT_CALLS
    assert b["reserved"] == pilot.FROZEN["max_calls"] - 1
    assert art["frozen_config"] == pilot.FROZEN


def test_the_published_decision_is_the_frozen_rule_bit_exactly():
    """The boundary case the sizing was built for: k = 1 FEASIBLE by
    the pre-frozen margin, reproduced bit-for-bit by decide(1)."""

    art = _artifact()
    assert art["decision"] == pilot.decide(1)
    assert art["decision"]["verdict"] == "FEASIBLE"
    tail = art["decision"]["tail_p_u_gt_u_max"]
    assert tail <= pilot.FROZEN["tail_budget"]
    assert 0 < pilot.FROZEN["tail_budget"] - tail < 5e-9


def test_the_reservation_is_recorded_and_the_ref_retained():
    res = _artifact()["reservation"]
    assert res["ref"] == "refs/s6m14pilot/reservation"
    assert res["object"] == "3fc1b28bdf8014245078f2ac3fd08faa17ea1f5b"
    assert res["seed_spent"] is True
    assert res["verified_at_exit"] is True


def test_the_seed_is_the_retired_stream_and_only_that():
    art = _artifact()
    assert art["seed"] == {"s6_m14_pilot": 40_000_461}
    assert ledger.OBSERVED_PROBE_SCALARS["s6_m14_pilot"] == 40_000_461
    assert "s6_m14_pilot" not in ledger.FRESH_PROBE_SCALARS


def test_the_final_checkpoint_is_a_non_verdict_partial_at_full_n():
    ck = json.loads(_CK.read_text(encoding="utf-8"))
    art = _artifact()
    assert ck["partial"] is True and ck["non_verdict"] is True
    assert ck["events_done"] == art["scan"]["events_done"]
    assert ck["k_ambiguous"] == art["scan"]["k_ambiguous"]
    assert ck["calls"] == art["scan"]["calls"]
    assert ck["rng_stream"] == "s6_m14_pilot"
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
    hist = _blob("docs/prereg/p14_s6_m14_pilot_freeze_manifest.json")
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
    assert '"s6_m14_pilot": 40_000_461' in frozen
