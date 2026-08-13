"""O4 abort-record contract tests.

An incident record is a dangerous document: it sits where a result
would sit, so the tests here exist mostly to pin what it must NOT say.
The campaign completed G1's and G2's execution and then aborted in G3
without persisting any statistic, so the record must carry a null
verdict, no gate status, and no numeric estimate -- and it must remain
lineage-verifiable against the freeze commit and the reservation the
run actually claimed.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

_PREREG = _REPO / "docs" / "prereg"
_INCIDENT = _PREREG / "p14_o4_incident.json"
_RESERVATION = _PREREG / "p14_o4_reservation.json"
_EXECUTED = _PREREG / "p14_o4_executed_freeze_manifest.json"

_FREEZE_REV = "1eb9461f1af739968403c633a7a44cbba9a9f948"
_RESERVATION_OBJ = "c4da1626463e6a6505374813cf3f56d6b429c209"


def _incident() -> dict:
    return json.loads(_INCIDENT.read_text(encoding="utf-8"))


# --------------------------------------------- the negative contract

def test_the_verdict_is_null_and_no_gate_has_a_status():
    """The whole point: an abort publishes no scientific verdict."""

    rec = _incident()
    assert rec["verdict"] is None
    assert rec["termination_reason"] == "g3-undecided"
    assert rec["exit_code"] == 1
    assert rec["results_artifact"] is None
    for gate in ("g1", "g2", "g3"):
        assert rec["gates"][gate]["status"] == "unavailable", gate


def test_g1_and_g2_completed_execution_but_carry_no_statistics():
    """Function completion is not gate passage. The distinction is the
    one an audit of this kind exists to respect, and an earlier working
    report blurred it."""

    gates = _incident()["gates"]
    assert gates["g1"]["execution"] == "completed"
    assert gates["g2"]["execution"] == "completed"
    assert gates["g3"]["execution"] == "aborted"
    for gate in ("g1", "g2"):
        assert "unavailable" in gates[gate]["statistics"]


def test_the_record_states_no_estimate_of_anything_measured():
    """A structural check rather than a wording check: no key anywhere
    in the record may carry a G1/G2/G3 estimate. If a later edit adds
    `mean_z` or `v_s1_lo` or a `leak_upper_abs`, this fails."""

    forbidden = ("mean", "var", "half_width", "v_s1", "identified",
                 "leak", "cp_upper", "concordant", "discordant",
                 "power", "band")
    seen: list[str] = []

    def walk(node, path=""):
        if isinstance(node, dict):
            for k, v in node.items():
                low = k.lower()
                if any(f in low for f in forbidden):
                    seen.append(f"{path}.{k}")
                walk(v, f"{path}.{k}")
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")

    walk(_incident())
    assert not seen, f"the incident record carries estimates: {seen}"


def test_the_record_says_out_loud_what_it_does_not_claim():
    claims = _incident()["what_is_not_claimed"]
    joined = " ".join(claims).lower()
    assert "no gate passed" in joined or "that any gate passed" in joined
    assert "concordant" in joined
    assert any("cap" in c for c in claims)


# ------------------------------------------------------- the lineage

def test_lineage_ties_the_record_to_the_freeze_and_the_reservation():
    rec = _incident()["lineage"]
    assert rec["freeze_rev"] == _FREEZE_REV
    assert rec["reservation_object"] == _RESERVATION_OBJ
    probe = subprocess.run(["git", "cat-file", "-e", _FREEZE_REV],
                           cwd=_REPO, capture_output=True)
    assert probe.returncode == 0, "the freeze commit must be reachable"


def test_the_reservation_record_matches_the_claim_the_run_made():
    res = json.loads(_RESERVATION.read_text(encoding="utf-8"))
    assert res["freeze_rev"] == _FREEZE_REV
    assert res["reservation_object"] == _RESERVATION_OBJ
    assert res["seeds"] == {"g1": 40_000_281, "g2": 40_000_291}
    assert "g3" not in res["seeds"]


def test_the_reservation_ref_is_retained_not_deleted():
    """The claim is permanent. Deleting the ref would let the retired
    streams read as free from a fresh clone."""

    lineage = _incident()["lineage"]
    assert lineage["reservation_ref"] == "refs/o4/reservation"
    assert lineage["reservation_ref_retained"] is True


def test_the_executed_freeze_manifest_is_the_frozen_blob_verbatim():
    """What the campaign actually verified against, preserved even as
    the current manifest moves on (the S4/S5 executed-snapshot split)."""

    blob = subprocess.run(
        ["git", "show", f"{_FREEZE_REV}:docs/prereg/"
         f"p14_o4_freeze_manifest.json"],
        cwd=_REPO, capture_output=True, check=True).stdout
    assert _EXECUTED.read_bytes() == blob


def test_the_result_artifact_does_not_exist():
    assert not (_PREREG / "p14_o4_results.json").exists()


# ---------------------------------------------------------- the seeds

def test_the_drawn_scalars_are_retired_under_abort_revealing_names():
    from probe_seed_ledger import (
        FRESH_PROBE_SCALARS,
        OBSERVED_PROBE_SCALARS,
        spent_scalars,
    )
    assert OBSERVED_PROBE_SCALARS["o4_aborted_g1"] == 40_000_281
    assert OBSERVED_PROBE_SCALARS["o4_aborted_g2"] == 40_000_291
    # the names must not read as the provenance of a verdict
    assert "o4_campaign" not in OBSERVED_PROBE_SCALARS
    assert "g1_audit" not in FRESH_PROBE_SCALARS
    assert "g2_leakage" not in FRESH_PROBE_SCALARS
    assert {40_000_281, 40_000_291} <= spent_scalars()
    assert _incident()["seeds"]["spent"] == {
        "o4_aborted_g1": 40_000_281, "o4_aborted_g2": 40_000_291}


def test_a_retired_stream_cannot_be_re_entered():
    """`assert_fresh_scalar` must refuse by name, with the remedy in
    the message -- a new allocation, or an explicitly labelled replay."""

    import pytest
    from probe_seed_ledger import assert_fresh_scalar, replay_scalar
    for name in ("g1_audit", "g2_leakage"):
        with pytest.raises(KeyError, match="never re-entered"):
            assert_fresh_scalar(name)
    assert replay_scalar("o4_aborted_g1") == 40_000_281


def test_the_never_drawn_scalar_stays_unspent():
    from probe_seed_ledger import OBSERVED_PROBE_SCALARS, spent_scalars
    assert 40_000_301 not in spent_scalars()
    assert 40_000_301 not in OBSERVED_PROBE_SCALARS.values()
    assert "40000301" in _incident()["seeds"]["not_spent"]


# ------------------------------------------- observability and timing

def test_the_timing_separates_machine_stamps_from_operator_records():
    """The runner published nothing, so the only absolute times come
    from git and from the operator's launch wrapper. The record must
    not present the latter as the runner's own clock."""

    t = _incident()["timing"]
    assert t["git_stamped"]["reservation_commit"].startswith("2026-08-12")
    assert "NOT the runner" in t["operator_record"]["source"]
    assert "perf_counter" in t["runner_log"]["source"]
    assert t["caps"]["bound"] is False


def test_the_runner_observability_defects_are_named():
    rec = _incident()
    defects = " ".join(rec["observability_defects"]).lower()
    assert "no g3 stress-point coordinates" in defects
    assert "publishes nothing on a fail-closed abort" in defects
    assert "no record is left" in defects
    required = " ".join(rec["required_before_any_rerun"]).lower()
    assert "write-once incident artifact" in required
    assert "new freeze and new g1/g2 scalars" in required


def test_the_inference_is_labelled_as_inference():
    inf = _incident()["inference_not_record"]
    assert "NOT a recorded fact" in inf["status"]
    assert "1245.8" in inf["basis"]


# ----------------------------------------------------- the derivation

def test_the_undecided_algebra_is_pinned_against_the_frozen_code():
    """The document's central claim is an identity, not a measurement:
    the outside probe's second leg sits a hard-coded 1e-6 from the edge
    it tests, so its undecided condition is `err >= 1e-6` and involves
    no geometry. Pin the three code facts it rests on, so a later edit
    to either file cannot silently falsify the write-up."""

    doc = (_PREREG / "p14_o4_incident.md").read_text(encoding="utf-8")
    plain = " ".join(doc.replace("**", "").split())
    assert "err₂ ≥ 1e-6" in plain
    assert "does not involve the geometry at all" in plain

    s1 = (_REPO / "experiments" / "positive_control"
          / "s1_schwarzschild_cost.py").read_text(encoding="utf-8")
    assert "if abs(dt - t_min) <= err:" in s1
    assert "return None" in s1.split("if abs(dt - t_min) <= err:")[1][:40]

    runner = (_REPO / "experiments" / "oracle"
              / "o4_volume_audit.py").read_text(encoding="utf-8")
    g3 = runner.split("def run_g3")[1]
    assert "(hi_edge + off, False)" in g3
    assert "0.5 * (lo_edge + hi_edge), True" in g3
    assert '"g3_stress_offset": 1e-6' in runner


def test_the_incident_document_refuses_to_grade_s1():
    doc = (_PREREG / "p14_o4_incident.md").read_text(encoding="utf-8")
    plain = " ".join(doc.replace("**", "").split())
    assert "not evidence about S1's correctness" in plain
    assert "Whether S1 is defective remains open" in plain
    assert "a replay is a reproduction, not an observation" in plain
