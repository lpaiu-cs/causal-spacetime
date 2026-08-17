"""O5 Poisson-count campaign result contract tests.

The campaign was executed ONCE, under the integrated execution
approval, from the clean exact checkout of the freeze BRANCH HEAD
`c6eb85e` -- the second parent of the merge commit. The committed
artifact is that single observation: CONCORDANT at K_certain = 53,285
with zero ambiguous points out of N = 26,831,117 sprinkled proposals.
These tests pin its lineage, the call accounting, the bit-exact
reproducibility of the published decision by the frozen general rule,
the containment arithmetic on the serialized values, the reservation
provenance, the retired seed, and the executed freeze surface against
the historical blobs at the freeze head."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
#: The approved execution SHA: the freeze branch head, NOT the merge.
_FREEZE = "c6eb85eea1affd33a1246627ea075dca96a859bb"
_MERGE = "26f8fcd0afc0b56b2ca7ab75e2c36f814de99fbb"
_ART = _REPO / "docs" / "prereg" / "p14_o5_count.json"
_CK = _REPO / "docs" / "prereg" / "p14_o5_count_checkpoint.json"
_SNAP = (_REPO / "docs" / "prereg"
         / "p14_o5_count_executed_freeze_manifest.json")

sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o5_count_campaign as camp  # noqa: E402
import probe_seed_ledger as ledger  # noqa: E402


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
    assert code["start"]["dirty"] is False
    assert code["end"]["dirty"] is False
    merge_parent2 = subprocess.run(
        ["git", "rev-parse", f"{_MERGE}^2"], cwd=_REPO,
        capture_output=True, text=True, check=True).stdout.strip()
    assert merge_parent2 == _FREEZE
    art = _artifact()
    assert art["freeze_sha"] == _FREEZE
    assert art["run_kind"] == "o5_count"
    assert art["frozen_config"] == camp.FROZEN


def test_the_scan_ran_the_whole_sprinkling_within_the_caps():
    """N drawn once, every proposal decided, zero ambiguous -- and
    the call accounting closes: between 1 and 2 calls per point, the
    budget reserved exactly the scan calls plus the frozen G3a
    preflight, and nothing was raised."""

    art = _artifact()
    scan = art["scan"]
    assert scan["n_total"] == 26_831_117
    assert scan["points_done"] == scan["n_total"]
    assert scan["k_certain"] == 53_285
    assert scan["u_ambiguous"] == 0
    assert scan["calls"] == 29_335_028
    assert scan["n_total"] <= scan["calls"] <= 2 * scan["n_total"]
    assert art["g3a_preflight"] == {"passed": True,
                                    "metered_calls": 6_537}
    b = art["budget"]
    assert b["reserved"] == scan["calls"] + 6_537 == 29_341_565
    assert b["reserved"] <= camp.FROZEN["max_calls"]
    assert b["wall_s"] <= camp.FROZEN["max_wall_s"]


def test_the_published_decision_is_the_frozen_rule_bit_exactly():
    art = _artifact()
    assert art["decision"] == camp.decide(53_285, 0)
    assert art["decision"]["verdict"] == "CONCORDANT"
    k_lo, k_hi = camp.ACCEPTANCE
    assert k_lo <= art["scan"]["k_certain"] <= k_hi


def test_the_containment_holds_on_the_serialized_values():
    """CONCORDANT means D is contained in [-B, B] -- re-checked on
    the exact committed numbers, not trusted from the label."""

    d = _artifact()["decision"]
    assert d["d_lo"] >= -d["band"]
    assert d["d_hi"] <= d["band"]
    assert d["d_lo"] < d["d_hi"]
    # and the interval arithmetic is self-consistent -- the
    # subtraction ran at 96 bits before serialization, so the double
    # recomputation agrees only to 1 ulp
    import math
    assert math.isclose(d["d_lo"],
                        d["c_lo_volume"] - d["oracle"]["v_hi"],
                        rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(d["d_hi"],
                        d["c_hi_volume"] - d["oracle"]["v_lo"],
                        rel_tol=0.0, abs_tol=1e-12)


def test_the_reservation_is_recorded_and_the_ref_retained():
    res = _artifact()["reservation"]
    assert res["ref"] == "refs/o5/reservation"
    assert res["object"] == "b5d9d6131644e7c8c07477c0b676696aaf2d1ce1"
    assert res["seed_spent"] is True
    assert res["verified_at_exit"] is True
    assert len(res["object"]) == 40


def test_the_seed_is_the_retired_stream_and_only_that():
    art = _artifact()
    assert art["seed"] == {"o5_campaign": 40_000_451}
    assert ledger.OBSERVED_PROBE_SCALARS["o5_campaign"] == 40_000_451
    assert "o5_campaign" not in ledger.FRESH_PROBE_SCALARS
    assert ledger.replay_scalar("o5_campaign") == 40_000_451


def test_the_final_checkpoint_is_a_non_verdict_partial_at_full_n():
    ck = json.loads(_CK.read_text(encoding="utf-8"))
    art = _artifact()
    assert ck["partial"] is True and ck["non_verdict"] is True
    assert ck["stage"] == "chunk"
    assert ck["n_total"] == art["scan"]["n_total"]
    assert ck["points_done"] == art["scan"]["n_total"]
    assert ck["k_certain"] == art["scan"]["k_certain"]
    assert ck["u_ambiguous"] == art["scan"]["u_ambiguous"]
    assert ck["calls"] == art["scan"]["calls"]
    assert ck["rng_stream"] == "o5_campaign"
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
    """The S4/S5/O3/O3'/pilot pattern: the snapshot must be
    byte-identical to the freeze-head blob of the manifest, and every
    digest it froze must equal the blob it pinned at that commit --
    including the ledger blob that carried the stream as FRESH and
    the three quoted artifacts (o3p, o3, pilot)."""

    hist = _blob("docs/prereg/p14_o5_count_freeze_manifest.json")
    assert _SNAP.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    assert set(m["files"]) == {
        "experiments/oracle/o5_count_campaign.py",
        "experiments/oracle/o5_reservation.py",
        "experiments/oracle/o4b_reservation.py",
        "experiments/oracle/o4b_budget.py",
        "experiments/oracle/o4b_meter.py",
        "experiments/oracle/o4b_g3a.py",
        "experiments/oracle/o4_g3_redesign.py",
        "experiments/oracle/o4_sizing.py",
        "experiments/oracle/o4_volume_audit.py",
        "experiments/positive_control/s1_schwarzschild_cost.py",
        "experiments/positive_control/probe_seed_ledger.py",
        "docs/prereg/p14_o5_count.md",
        "docs/prereg/p14_o3p_volume.json",
        "docs/prereg/p14_o3_volume.json",
        "docs/prereg/p14_o5_amb_pilot.json",
        "pyproject.toml",
    }
    assert set(m["files"]) == set(camp.PROTOCOL_SURFACE)
    for rel, want in m["files"].items():
        assert hashlib.sha256(_blob(rel)).hexdigest() == want, rel
