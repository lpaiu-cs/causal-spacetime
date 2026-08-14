"""The published census artifact: what it is, and what it is not.

The artifact is a REPLAY of a retired stream over a frozen stress set.
These tests hold it to that: it must cover the frozen set, separate
every (probe, leg, outcome), carry no key a reader could mistake for a
result, and leave the ledger and the reservation exactly as it found
them. They also pin the numbers the narrative quotes, so the two
cannot drift.

No test here asserts that anything passed. O4 has no verdict.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_g3_redesign as g3  # noqa: E402
import o4_replay_diagnostic as rd  # noqa: E402
import o4_volume_audit as o4  # noqa: E402
import probe_seed_ledger as ledger  # noqa: E402

_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o4_replay_diagnostic.json"
_NARRATIVE = (_REPO / "docs" / "prereg"
              / "p14_o4_replay_diagnostic_results.md")

_CODE_BASIS = "1257c4b"

#: Key substrings that would mean a scientific quantity had leaked in.
_FORBIDDEN = ("mean", "var", "half_width", "v_s1", "identified",
              "leak", "cp_upper", "concordant", "discordant", "power",
              "band", "verdict", "gate", "status", "estimate")


def _artifact() -> dict:
    return json.loads(_ARTIFACT.read_text(encoding="utf-8"))


def _prose() -> str:
    text = _NARRATIVE.read_text(encoding="utf-8")
    lines = [ln.lstrip("> ") for ln in text.replace("**", "").split("\n")]
    return " ".join(" ".join(lines).split())


def _walk_keys(node, path=""):
    if isinstance(node, dict):
        for key, value in node.items():
            yield f"{path}.{key}", key
            yield from _walk_keys(value, f"{path}.{key}")
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _walk_keys(value, f"{path}[{i}]")


# ------------------------------------------------------- what it is

def test_the_artifact_is_a_replay_over_the_whole_frozen_set():
    art = _artifact()
    assert art["run_kind"] == "replay"
    assert art["replayed_stream"]["seed"] == 40_000_281
    found = art["findings"]
    assert found["covers_frozen_stress_set"] is True
    assert found["clusters_probed"] == o4.FROZEN["g3_clusters"] == (
        100_000)


def test_the_narrative_names_the_code_it_ran_on():
    assert _CODE_BASIS in _prose()


def test_the_readout_never_disagreed_with_the_predicate():
    """A non-zero count would mean the census misread the calls it was
    reporting on."""

    found = _artifact()["findings"]
    assert found["rederivation_disagreements"] == 0
    assert found["census"]["dpsi_leg_disagreements"] == 0


# --------------------------------------------------- what it is not

def test_no_key_could_be_read_as_a_result():
    offenders = [(where, key) for where, key in _walk_keys(_artifact())
                 if any(bad in key.lower() for bad in _FORBIDDEN)]
    assert offenders == []


def test_the_narrative_awards_no_gate_and_no_verdict():
    plain = _prose()
    assert "게이트도, 과학적 verdict 도 부여하지 않는다" in plain
    assert "무판정(`verdict = null`)" in plain
    assert "독립 확인도 아니다" in plain
    for forbidden in ("CONCORDANT", "DISCORDANT", "통과했다"):
        assert forbidden not in plain


def test_the_census_did_not_move_the_ledger():
    """A replay is not an observation: nothing is spent by running it."""

    assert ledger.FRESH_PROBE_SCALARS == {}
    assert 40_000_301 not in ledger.spent_scalars()
    assert ledger.replay_scalar("o4_aborted_g1") == 40_000_281
    assert "40,000,301" in _prose() and "c4da162" in _prose()


def test_the_frozen_constants_were_not_retuned_to_the_result():
    plain = _prose()
    assert g3.ETA == 1e-12
    assert not hasattr(g3, "MAX_NUDGES")   # removed, not retuned
    assert "이 결과에 맞춰 아무 규칙도 바꾸지 않았다" in plain
    assert "η 와 `MAX_NUDGES` 는 여기 없다" in plain


# ------------------------------------------- the counts, pinned once

def _outcome_row(probe: str, leg: str) -> dict:
    rows = _artifact()["findings"]["census"][
        "outcomes_by_probe_and_leg"]
    return next(r for r in rows
                if r["probe"] == probe and r["leg"] == leg)


def test_only_one_of_the_twelve_cells_is_undecided():
    """The single realized cause. Every other cell is unanimous."""

    for probe, leg in ((rd.MIDPOINT, rd.LEG_PX),
                       (rd.MIDPOINT, rd.LEG_XQ),
                       (rd.OUTSIDE, rd.LEG_PX)):
        row = _outcome_row(probe, leg)
        assert row["undecided"] == 0
        assert row["true"] == 100_000

    row = _outcome_row(rd.OUTSIDE, rd.LEG_XQ)
    assert row["undecided"] == 626
    assert row["false"] == 99_374
    assert row["true"] == 0


def test_the_midpoint_path_never_fired_and_nothing_mismatched():
    found = _artifact()["findings"]
    counts = {c["cause"]: c["clusters"] for c in found["counts"]}
    assert counts[rd.CAUSE_MIDPOINT] == 0
    assert counts[rd.CAUSE_OUTSIDE] == 626
    assert counts[rd.CAUSE_MISMATCH] == 0
    assert found["clusters_with_no_cause"] == 99_374
    for row in found["census"]["boolean_mismatch_by_probe"]:
        assert row["mismatching_probes"] == 0


def test_the_first_undecided_is_cluster_604():
    first = _artifact()["findings"]["first_undecided"]
    assert first["cluster_index"] == 604
    assert first["cause"] == rd.CAUSE_OUTSIDE
    site = first["site"]
    assert site["probe"] == rd.OUTSIDE and site["leg"] == rd.LEG_XQ
    assert site["decision_margin"] < 0.0
    assert site["err"] > o4.FROZEN["g3_stress_offset"]
    assert "cluster 604" in _prose()


def test_nominal_and_realized_agreed_here_but_are_still_counted_apart():
    row = _artifact()["findings"]["census"][
        "outside_xq_err2_at_least_frozen_offset"]
    assert row["nominal"] == row["realized_undecided"] == 626
    assert row["disagreements"] == 0
    assert row["first_disagreement"] is None
    assert "원리적으로 가능한 것과 실제로 일어난 것은 다르며" in _prose()


def test_err_is_not_bounded_by_tol():
    """The failure's shape, measured: clearing `tol` by two decades
    buys nothing against `err`, whose tail runs five decades past it."""

    dists = _artifact()["findings"]["census"]["distributions"]
    assert o4.FROZEN["tol"] == 1e-8
    assert dists["err1"]["min"] >= o4.FROZEN["tol"]
    assert dists["err2"]["max"] > 1e-3
    assert dists["err2"]["max"] / o4.FROZEN["tol"] > 1e4
    assert dists["err2"]["max"] > o4.FROZEN["g3_stress_offset"] * 1000


def test_every_cluster_was_eligible_at_every_eta_on_this_set():
    rows = _artifact()["findings"]["census"]["eligibility_by_eta"]
    assert [r["eta"] for r in rows] == list(g3.ETA_GRID)
    for row in rows:
        assert row["eligible_clusters"] == 100_000
        assert row["lower_probe_t_x_negative"] == 0
        assert row["eligible_and_lower_probe_t_x_negative"] == 0
    assert "설계 입력일 뿐이다" in _prose()


def test_the_angle_term_was_zero_here_and_that_is_not_a_licence():
    """Zero by luck of the draw is exactly the kind of margin this
    stage keeps finding; the redesign makes it zero by construction."""

    dists = _artifact()["findings"]["census"]["distributions"]
    for name in ("abs_t_min_minus_T_in", "abs_t_min_minus_T_out"):
        assert dists[name]["min"] == dists[name]["max"] == 0.0
    plain = _prose()
    assert "그 항은 없다" in plain
    assert "뽑기 운이 아니라 구성이어야" in plain
