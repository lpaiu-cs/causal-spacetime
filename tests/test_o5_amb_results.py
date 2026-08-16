"""O5 ambiguity-pilot result contract tests.

The pilot was executed ONCE, from the clean exact checkout of the
freeze BRANCH HEAD `e9e7e15` -- the second parent of the merge commit,
per the delegated execution ruling -- under the frozen caps. The
committed artifact is that single observation: FEASIBLE at
k_ambiguous = 0 of 10,736,965 events. These tests pin its lineage, the
call accounting, the bit-exact reproducibility of the published
decision by the frozen general-in-k rule, the reservation provenance,
the retired seed, and the executed freeze surface against the
historical blobs at the freeze head -- so the result stays checkable
long after the CURRENT manifest is re-pinned for later work.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
#: The approved execution SHA: the freeze branch head, NOT the merge.
_FREEZE = "e9e7e1592319989233f11b479586c603e8bad958"
_MERGE = "8e1d1c36be873c92791366c176c37da1a24ff938"
_ART = _REPO / "docs" / "prereg" / "p14_o5_amb_pilot.json"
_CK = _REPO / "docs" / "prereg" / "p14_o5_amb_pilot_checkpoint.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_o5_amb_pilot_executed_freeze_manifest.json")

sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o5_amb_pilot as pilot  # noqa: E402
import probe_seed_ledger as ledger  # noqa: E402


def _artifact() -> dict:
    return json.loads(_ART.read_text(encoding="utf-8"))


def _blob(rel: str) -> bytes:
    return subprocess.run(
        ["git", "cat-file", "blob", f"{_FREEZE}:{rel}"], cwd=_REPO,
        capture_output=True, check=True).stdout


def test_lineage_is_the_freeze_branch_head_clean_at_entry_and_exit():
    """One run, from the clean exact checkout of the freeze BRANCH
    HEAD (the merge commit's second parent, never the merge commit) --
    start and end revs both that commit, clean, and the exit rev
    compared to the approved SHA by the runner itself."""

    code = _artifact()["code"]
    assert code["start"]["rev"] == _FREEZE
    assert code["end"]["rev"] == _FREEZE
    assert code["start"]["dirty"] is False
    assert code["end"]["dirty"] is False
    merge_parent2 = subprocess.run(
        ["git", "rev-parse", f"{_MERGE}^2"], cwd=_REPO,
        capture_output=True, text=True, check=True).stdout.strip()
    assert merge_parent2 == _FREEZE
    art = _artifact()
    assert art["freeze_sha"] == _FREEZE
    assert art["run_kind"] == "o5_amb_pilot"


def test_the_scan_ran_the_full_fixed_n_with_exact_call_accounting():
    """Fixed n, never extended, never cut short: every event decided
    both legs (zero ambiguous means zero early exits), so the scan's
    metered calls are exactly 2n, and with G3a's frozen preflight the
    budget closed exactly at max_calls -- reserved, not overshot."""

    art = _artifact()
    scan = art["scan"]
    n = pilot.FROZEN["n_events"]
    assert scan["events_done"] == n == 10_736_965
    assert scan["k_ambiguous"] == 0
    assert scan["calls"] == 2 * n == 21_473_930
    assert art["g3a_preflight"] == {"passed": True,
                                    "metered_calls": 6_537}
    b = art["budget"]
    assert b["reserved"] == b["completed"] == pilot.FROZEN["max_calls"]
    assert b["reserved"] == 2 * n + 6_537
    assert b["wall_s"] <= pilot.FROZEN["max_wall_s"]
    assert art["frozen_config"] == pilot.FROZEN


def test_the_published_decision_is_the_frozen_rule_bit_exactly():
    """Re-applying the frozen general-in-k rule to the observed k
    reproduces the published decision block verbatim -- the same
    recoverability the O4b audit had to prove after the fact, proved
    here while the artifact is fresh."""

    art = _artifact()
    assert art["decision"] == pilot.decide(0)
    assert art["decision"]["verdict"] == "FEASIBLE"
    assert art["decision"]["tail_p_u_gt_u_max"] <= (
        pilot.FROZEN["tail_budget"])


def test_the_reservation_is_recorded_and_the_ref_retained():
    """The claim object is the artifact's provenance nonce; the ref
    stays retained (never deleted), like refs/o4 and refs/o4b."""

    res = _artifact()["reservation"]
    assert res["ref"] == "refs/o5pilot/reservation"
    assert res["object"] == "f5223636bc7c60603b18ecae3f474dd3dc37146f"
    assert res["seed_spent"] is True
    assert res["verified_at_exit"] is True
    assert len(res["object"]) == 40


def test_the_seed_is_the_retired_stream_and_only_that():
    art = _artifact()
    assert art["seed"] == {"o5_amb_pilot": 40_000_441}
    assert ledger.OBSERVED_PROBE_SCALARS["o5_amb_pilot"] == 40_000_441
    assert "o5_amb_pilot" not in ledger.FRESH_PROBE_SCALARS


def test_the_final_checkpoint_is_a_non_verdict_partial_at_full_n():
    """The last chunk's atomic checkpoint is preserved as committed
    provenance: writer-stamped partial/non_verdict (no count in it is
    a verdict -- the artifact is), at exactly the full n and the same
    freeze lineage as the artifact."""

    ck = json.loads(_CK.read_text(encoding="utf-8"))
    art = _artifact()
    assert ck["partial"] is True and ck["non_verdict"] is True
    assert ck["stage"] == "chunk"
    assert ck["events_done"] == art["scan"]["events_done"]
    assert ck["k_ambiguous"] == art["scan"]["k_ambiguous"]
    assert ck["calls"] == art["scan"]["calls"]
    assert ck["rng_stream"] == "o5_amb_pilot"
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
    """The S4/S5/O3/O3' pattern: the snapshot must be byte-identical
    to the freeze-head blob of the manifest, and every digest it froze
    must equal the blob it pinned at that commit -- including the
    ledger blob that carried the stream as FRESH."""

    hist = _blob("docs/prereg/p14_o5_amb_pilot_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    # the EXECUTED surface, pinned literally: the live
    # PROTOCOL_SURFACE may grow (it has -- p14_o3_volume.json was
    # added post-run), but what this run depended on is history
    assert set(m["files"]) == {
        "experiments/oracle/o5_amb_pilot.py",
        "experiments/oracle/o5_pilot_reservation.py",
        "experiments/oracle/o4b_reservation.py",
        "experiments/oracle/o4b_budget.py",
        "experiments/oracle/o4b_meter.py",
        "experiments/oracle/o4b_g3a.py",
        "experiments/oracle/o4_g3_redesign.py",
        "experiments/oracle/o4_sizing.py",
        "experiments/oracle/o4_volume_audit.py",
        "experiments/positive_control/s1_schwarzschild_cost.py",
        "experiments/positive_control/probe_seed_ledger.py",
        "docs/prereg/p14_o5_amb_pilot.md",
        "docs/prereg/p14_o3p_volume.json",
        "pyproject.toml",
    }
    assert set(m["files"]) <= set(pilot.PROTOCOL_SURFACE)
    for rel, want in m["files"].items():
        assert hashlib.sha256(_blob(rel)).hexdigest() == want, rel
