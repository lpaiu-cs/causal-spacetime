"""The G3 redesign design document: what it fixes, and that the code
agrees with it.

The document is the thing being frozen here -- there is no G3a/G3b
implementation yet and no campaign to run. So these tests do two jobs:
they pin the decisions the document is not allowed to lose, and they
tie every frozen number in it to the constant the census actually uses,
so a later edit cannot leave the two disagreeing.

Nothing here asserts a measured value. O4 remains verdict-free.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_g3_redesign as g3  # noqa: E402
import o4_replay_diagnostic as rd  # noqa: E402
import o4_sizing as sz  # noqa: E402
import o4_volume_audit as o4  # noqa: E402

_DOC = _REPO / "docs" / "prereg" / "p14_o4_g3_redesign.md"


def _prose() -> str:
    text = _DOC.read_text(encoding="utf-8")
    lines = [ln.lstrip("> ") for ln in text.replace("**", "").split("\n")]
    return " ".join(" ".join(lines).split())


# ------------------------------------------------ what it must not lose

def test_the_document_exists_and_publishes_no_result():
    plain = _prose()
    assert "결과 없음" in plain
    assert "verdict = null" in plain
    for forbidden in ("CONCORDANT", "DISCORDANT", "V_S1 ="):
        assert forbidden not in plain


def test_g3a_treats_a_legitimate_none_as_a_pass():
    """The contract fix. The frozen G3 aborted on any `None`; a
    tri-state predicate must be allowed to be undecided where it
    should be."""

    plain = _prose()
    assert "정당한 undecided" in plain
    assert "`None` 이 없다" in plain
    assert "`None` 이 나와야 할 때만" in plain


def test_g3a_row_b_is_separated_from_the_negative_dt_row():
    """Otherwise row B silently becomes a duplicate of row D: at a
    point where `t_min - err - eta < 0` the predicate short-circuits
    and row B passes without ever exercising the lower error-band
    path it exists to test."""

    plain = _prose()
    assert "행 B 와 행 D 는 겹치면 안 된다" in plain
    assert "construction -unavailable" in plain or (
        "construction-unavailable" in plain)


def test_g3a_runs_before_any_fresh_seed_is_touched():
    """The direct lesson of the abort: the frozen runner first touched
    the wrapper contract after spending twelve hours of sample."""

    plain = _prose()
    assert "fresh 시드를 건드리기 전에" in plain
    assert "통과한 경우에만 G1 / G2 를 실행한다" in plain


def test_g3b_states_the_three_probes_and_their_proofs():
    plain = _prose()
    for probe in ("determined-inside", "determined-outside-above",
                  "determined-outside-below"):
        assert probe in plain
    assert "(iii) 은 아래쪽 경계를 처음으로 검사한다" in plain


def test_the_lower_probe_separates_its_three_outcomes():
    plain = _prose()
    for outcome in ("lower-boundary-solver-probe",
                    "negative-dt-short-circuit",
                    "construction-unavailable"):
        assert outcome in plain
    assert "G3b 의 통과 증거로 세지 않는다" in plain


def test_availability_is_reported_now_and_gated_later():
    plain = _prose()
    assert "census 에는 gate 를 붙이지 않는다" in plain
    assert "최소 적격률" in plain
    assert "INVALID" in plain and "INCONCLUSIVE" in plain
    assert ("이전 stress set 의 availability 는 설계 입력일 뿐이며, 새 "
            "캠페인의 통과 증거로 재사용하지 않는다") in plain


def test_the_self_selection_of_eligible_points_is_admitted():
    plain = _prose()
    assert "자기선택된 부분집합" in plain
    assert "availability rate 를 함께 게시한다" in plain


def test_eta_is_a_realized_margin_not_an_addend():
    plain = _prose()
    assert "실현된 최소 판정 여유" in plain
    assert "realized_margin = |dt − t_min| − err" in plain
    assert "휴리스틱에 인증을 맡기지" in plain
    assert "nextafter" in plain
    assert "`err₂ ≥ 10⁻⁶` 의 빈도로 정하지 않는다" in plain


def test_the_angle_recovery_finding_is_recorded():
    plain = _prose()
    assert "ulp 단위로 유계가 아니고" in plain
    assert "술어 좌표에서 만든다" in plain


def test_using_err_to_place_probes_is_distinguished_from_inference():
    plain = _prose()
    assert "probe 를 배치" in plain
    assert "V_S1 의 구간을 만드는 것" in plain or (
        "구간을 만드는 것" in plain)


# ----------------------------------- the code agrees with the document

def test_eta_matches_the_constant_the_census_uses():
    assert g3.ETA == 1e-12
    assert rd.ETA is g3.ETA          # re-exported, never redefined
    assert "η = 10⁻¹² (frozen)" in _prose()


def test_the_outward_search_has_no_cap_and_says_why():
    """8,042 is evidence that an arbitrary cap is the defect, not that
    the case should be discarded. The cap is gone from the constant,
    the document and the verdict."""

    plain = _prose()
    assert not hasattr(g3, "MAX_NUDGES")
    assert "바깥 탐색에는 cap 이 없다" in plain
    assert "임의 cap 이 결함" in plain
    assert "거리는 진단값이며 비용이 아니다" in plain
    source = (_REPO / "experiments" / "oracle"
              / "o4_g3_redesign.py").read_text(encoding="utf-8")
    assert "There is deliberately NO cap" in source


def test_lower_probe_coverage_is_specified_as_a_joint():
    plain = _prose()
    assert "주변합만으로는 답이 나오지 않는다" in plain
    assert "eligible_and_lower_probe_t_x_negative" in plain
    assert "적격 cluster 에서만" in plain


def test_the_post_census_list_excludes_eta_and_has_no_cap_to_list():
    plain = _prose()
    assert "η 는 여기 없다" in plain
    assert "탐색 cap 도 여기 없다 — 아예 존재하지 않기 때문이다" in plain


def test_the_eight_thousand_case_is_a_SUCCESS_not_an_exclusion():
    """Run the loop: the margin IS reached, at 8,042 steps. The probe
    is perfectly constructible, and the first response -- capping at 64
    and recording it unavailable -- discarded a good probe."""

    t_min = 0.0026643057799644846
    err = 0.0026640756792607376
    dt = t_min - err - g3.ETA
    steps = 0
    while abs(dt - t_min) - err < g3.ETA:
        dt = math.nextafter(dt, -math.inf)
        steps += 1
        if steps > 10 ** 6:
            break
    assert steps == 8_042
    assert "8,042 걸음" in _prose()
    source = (_REPO / "experiments" / "oracle"
              / "o4_g3_redesign.py").read_text(encoding="utf-8")
    assert "8,042 steps" in source
    assert "151" not in source


def test_the_redesign_arithmetic_has_one_home():
    """A constant with two definitions is how a document and its code
    drift apart."""

    source = (_REPO / "experiments" / "oracle"
              / "o4_replay_diagnostic.py").read_text(encoding="utf-8")
    assert "ETA = g3.ETA" in source
    assert "ETA = 1e-12" not in source


def test_probe_times_place_each_probe_where_the_proofs_say():
    lo, hi, e1, e2 = 1.0, 7.0, 1e-6, 2e-6
    times = g3.probe_times(lo, hi, e1, e2)
    assert times["outside_above"] == hi + e2 + g3.ETA
    assert times["outside_below"] == lo - e1 - g3.ETA
    assert times["inside"] == 0.5 * ((lo + e1 + g3.ETA)
                                     + (hi - e2 - g3.ETA))
    # the inside probe's two margins sum to W_robust + 2*eta, so at the
    # midpoint each is at least eta exactly when the point is eligible
    margin_a = times["inside"] - lo - e1
    margin_b = (hi - times["inside"]) - e2
    assert margin_a + margin_b == pytest.approx(
        g3.w_robust(hi - lo, e1, e2) + 2 * g3.ETA, rel=1e-12)
    assert min(margin_a, margin_b) >= g3.ETA


def test_eligibility_is_strict_in_the_shared_module():
    assert g3.is_eligible(1.0, 0.4, 0.4, 0.1) is False   # exactly 0
    assert g3.is_eligible(1.0, 0.4, 0.4, 0.09) is True


def test_the_eta_grid_in_the_document_is_the_grid_in_the_code():
    text = _DOC.read_text(encoding="utf-8")
    listed = ("η ∈ {0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, "
              "1e-9, 1e-8, 1e-7, 1e-6}")
    assert listed in " ".join(text.split())
    assert rd.ETA_GRID == (0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11,
                           1e-10, 1e-9, 1e-8, 1e-7, 1e-6)


def test_the_quantile_declaration_matches_the_code():
    plain = _prose()
    assert "{0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}" in plain
    assert 'method="linear"' in plain
    assert rd.QUANTILE_PROBS == (0.001, 0.01, 0.1, 0.5, 0.9, 0.99,
                                 0.999)
    assert rd.QUANTILE_METHOD == "linear"


def test_the_bin_convention_matches_the_code():
    plain = _prose()
    assert "왼쪽 닫힘·오른쪽 열림" in plain
    got = rd.histogram([0.0], (0.0, 1.0))
    assert "left-closed right-open" in got["bin_convention"]
    assert rd.ERR_EDGES[0] == 1e-16
    assert rd.ERR_EDGES[-1] == 1e-2
    assert rd.LEN_EDGES[0] == 0.0
    assert rd.LEN_EDGES[-1] == sz.L_MAX_UB


def test_the_strictness_of_the_eligibility_comparison_matches():
    assert "W_robust > 0" in _prose()
    rows = {r["eta"]: r for r in rd.eligibility(
        {"lo": 1.0, "err1": 0.0, "err2": 0.0, "L": 0.0,
         "L_minus_errs": 0.0})}
    assert rows[0.0]["w_robust"] == 0.0
    assert rows[0.0]["eligible"] is False


def test_the_frozen_offset_the_redesign_replaces_is_still_1e_6():
    """The redesign's adaptive offset only means something against the
    fixed one it replaces."""

    assert o4.FROZEN["g3_stress_offset"] == 1e-6
    assert "고정 `10⁻⁶`" in _prose()


def test_the_document_quotes_the_measured_angle_error():
    """The table that forced the design decision. Recompute it rather
    than trusting the transcription."""

    plain = _prose()
    assert "4.14e-13" in plain
    got = abs(math.acos(max(-1.0, min(1.0, math.cos(1e-5)))) - 1e-5)
    assert f"{got:.2e}" == "4.14e-13"


def test_the_search_is_exponential_but_lands_where_the_walk_would():
    """No cap, and no linear walk either: G3b places three probes at
    100,000 clusters, where walking thousands of ulps would be a new
    way to stall. The margin is monotone in float order, so bracket and
    bisect -- and land on the SAME first-satisfying value."""

    import o4b_g3a as g3a

    t_min, err = 0.0026643057799644846, 0.0026640756792607376
    walked, steps = t_min - err - g3.ETA, 0
    while abs(walked - t_min) - err < g3.ETA:
        walked = math.nextafter(walked, -math.inf)
        steps += 1

    placed = g3a.place_row(t_min, err, g3.ETA, above=False)
    assert placed["reached"] is True          # a SUCCESS, not discarded
    assert placed["dt"] == walked
    assert placed["ulp_distance"] == steps == 8_042
    assert "bisection" in placed["search"]


def test_the_lower_probe_stops_at_zero_not_at_a_budget():
    """Termination is on representability. `t_min = 0` leaves nowhere
    below to go, and that is the only reason row B gives up there."""

    import o4b_g3a as g3a

    placed = g3a.place_row(0.0, 0.0, g3.ETA, above=False)
    assert placed["reached"] is False
    assert "short-circuit" in placed["why"]


def test_the_inside_probe_may_not_use_the_monotone_search():
    """Widening one leg narrows the other, so the success region is an
    interval and an exponential jump can step over it. The constraint
    has to be recorded before G3b is built, not discovered in it."""

    plain = _prose()
    assert "성공영역이 반직선이 아니라 구간" in plain
    assert "두 leg 의 표현 가능 범위 교집합을 직접 구하고" in plain
    assert "호출해서는 안 된다" in plain

    import o4b_g3a as g3a
    assert "must NOT call this function" in g3a.place_row.__doc__
    assert "MONOTONICITY IS A PRECONDITION" in g3a.place_row.__doc__


def test_the_only_accepted_omissions_are_the_seven_named_ones():
    """Frozen as a list of labels, and the run has to agree with it.

    Seven of the eleven specs are found by bisecting the solver's
    equal-perihelion band, and equal radii have no such band."""

    import o4b_g3a as g3a

    assert g3.EXPECTED_UNREACHABLE_LABELS == {
        f"equal-radius/{s}" for s in _BAND_SPECS}

    table = g3a.run_g3a()
    assert {u["case"] for u in table["unreachable"]} == (
        g3.EXPECTED_UNREACHABLE_LABELS)
    assert table["unexpected_unreachable"] == []
    assert table["unexpectedly_reachable"] == []
    assert table["case_build_failures"] == []
    assert all(u["expected"] for u in table["unreachable"])


#: The band-dependent specs that an equal-radius pair cannot host.
#: `one-turn-beyond` is deliberately NOT here -- it needs the band's
#: upper edge, which is known in closed form, and the angle above it
#: is reachable.
_BAND_SPECS = ("no-turn-mid", "equal-perihelion-inside",
               "equal-perihelion-lower-edge-below",
               "equal-perihelion-lower-edge-above",
               "equal-perihelion-upper-edge-below",
               "equal-perihelion-upper-edge-above")

#: Every spec whose construction goes through the band search, which
#: is what a broken search takes down at a geometry that has one.
_ALL_BAND_SPECS = _BAND_SPECS + ("one-turn-beyond",)


def _flaky_band(monkeypatch, radii: tuple[float, float]):
    """Make the band search fail with a plain `ValueError` at one
    geometry, the way a bracket or a bisection would."""

    real = g3.equal_perihelion_band

    def flaky(r1, r2, m, tol, steps=200):
        if (r1, r2) == radii:
            raise ValueError("bisection never landed in the band")
        return real(r1, r2, m, tol, steps)

    monkeypatch.setattr(g3, "equal_perihelion_band", flaky)


def test_a_search_that_fails_is_not_an_omission(monkeypatch):
    """The R2 defect. `equal_perihelion_band()` fails for two
    unrelated reasons, and swallowing both as `unreachable` dropped
    cases silently -- while the other geometries still covered all
    four families and no row had run to fail, so G3a would PASS with a
    hole in its table and fresh seeds would then be spent."""

    import o4b_g3a as g3a

    _flaky_band(monkeypatch, (13.0, 17.0))
    preflight = g3a.run_preflight()

    table = preflight["table"]
    lost = {f["case"] for f in table["case_build_failures"]}
    assert lost == {f"increasing-interior/{s}" for s in _ALL_BAND_SPECS}
    # the hole would not have shown up in any other condition
    assert table["covers_every_family"]
    assert not table["failures"]
    assert preflight["conditions"]["every_case_built"] is False
    assert preflight["passed"] is False
    assert "every_case_built" in preflight["failed_conditions"]


def test_a_search_that_fails_AT_AN_EXPECTED_NAME_is_still_a_failure(
        monkeypatch):
    """The R4 defect, and the reason the list of names is not enough.

    A frozen list separates the two reasons only by WHERE the failure
    happened. A bracket that went wrong on `equal-radius` lands on
    seven names already on the list, so the label check waves it
    through as the expected absence -- the branch goes unverified and
    G3a passes. Only the raiser knows which reason it had, so the
    geometry raises `ExpectedUnreachable` and everything else is a
    build failure regardless of where it happened."""

    import o4b_g3a as g3a

    def broken(r1, r2, m, tol):
        raise ValueError("the solver's small-angle structure moved")

    # equal radii do not reach the band SEARCH -- their band is known
    # in closed form -- so the analogous failure is their structure
    # check, which every band-dependent spec there goes through
    monkeypatch.setattr(g3, "verify_equal_radius_structure", broken)
    preflight = g3a.run_preflight()

    table = preflight["table"]
    failed = {f["case"] for f in table["case_build_failures"]}
    assert failed == {f"equal-radius/{s}" for s in _ALL_BAND_SPECS}
    listed = [f for f in table["case_build_failures"]
              if f["on_the_expected_list"]]
    assert len(listed) == len(_BAND_SPECS)
    # the label check alone is satisfied -- nothing is missing from
    # the list, and nothing on the list got built
    assert preflight["conditions"]["only_expected_unreachable"] is True
    assert preflight["conditions"]["every_case_built"] is False
    assert preflight["passed"] is False


def test_the_geometry_raises_its_own_exception_type():
    """`ExpectedUnreachable` is what "this branch does not exist here"
    is spelled as; a plain `ValueError` never means that."""

    import s1_schwarzschild_cost as s1

    with pytest.raises(g3.ExpectedUnreachable, match="zero arc"):
        g3.equal_perihelion_band(15.0, 15.0, s1.M, s1.DEFAULT_TOL)
    assert issubclass(g3.ExpectedUnreachable, Exception)
    assert not issubclass(g3.ExpectedUnreachable, ValueError)


def test_equal_radii_DO_have_a_band_it_is_simply_out_of_reach():
    """A correction to what this module said twice.

    "Equal radii are always one-turn" is false in `dpsi`: the solver
    reports `equal-perihelion` for every `dpsi <= tol`. What is true
    is that nothing lies BELOW the band, so no-turn never occurs and
    there is no transition to bracket -- and that the band sits below
    every angle the wrapper can hand over, so no case in theta space
    could land in it either."""

    import math

    import s1_schwarzschild_cost as s1

    assert g3.family_at(15.0, 15.0, 1e-9, s1.M,
                        s1.DEFAULT_TOL) == "equal-perihelion"
    assert not any(
        g3.family_at(15.0, 15.0, x, s1.M, s1.DEFAULT_TOL) == "no-turn"
        for x in (1e-13, 1e-11, 1e-9, 5e-9, 9e-9, 1.1e-8, 1e-7, 1e-3))

    smallest = math.acos(math.nextafter(1.0, -math.inf))
    assert smallest == 1.4901161193847656e-08
    assert smallest > s1.DEFAULT_TOL          # the band is unreachable
    reached = {g3.family_at(15.0, 15.0, g3.wrapper_dpsi(t), s1.M,
                            s1.DEFAULT_TOL)
               for t in (0.0, 1e-9, 1e-8, 1e-7, 1e-4, 1e-2, 0.5, 1.0)}
    assert reached == {"radial", "one-turn"}


def test_the_equal_radius_structure_is_checked_claim_by_claim(
        monkeypatch):
    """Four separate assertions, not one weak one (review R7).

    "The low probe is not no-turn" passes on `radial` too, so a solver
    whose small-angle behaviour had moved would go by as the expected
    structure. Each claim is therefore falsified on its own."""

    import s1_schwarzschild_cost as s1

    ok = g3.verify_equal_radius_structure(15.0, 15.0, s1.M,
                                          s1.DEFAULT_TOL)
    assert ok["families_seen"] == ["equal-perihelion", "one-turn",
                                   "radial"]
    assert ok["band"] == (0.0, s1.DEFAULT_TOL)
    assert ok["smallest_recoverable_dpsi"] > s1.DEFAULT_TOL

    real = g3.family_at

    def forced(want, at=None):
        def fam(r1, r2, x, m, tol):
            if at is None or x == at:
                return want
            return real(r1, r2, x, m, tol)
        return fam

    # the band is not there at all
    monkeypatch.setattr(g3, "family_at", forced("radial"))
    with pytest.raises(ValueError, match="should cover small positive"):
        g3.verify_equal_radius_structure(15.0, 15.0, s1.M,
                                         s1.DEFAULT_TOL)
    # above the band it is not one-turn
    monkeypatch.setattr(g3, "family_at", forced("radial", at=3.0))
    with pytest.raises(ValueError, match="should be one-turn"):
        g3.verify_equal_radius_structure(15.0, 15.0, s1.M,
                                         s1.DEFAULT_TOL)
    # a no-turn region appeared
    monkeypatch.setattr(g3, "family_at", forced("no-turn", at=1e-6))
    with pytest.raises(ValueError, match="no-turn should not occur"):
        g3.verify_equal_radius_structure(15.0, 15.0, s1.M,
                                         s1.DEFAULT_TOL)
    # the band is no longer out of the wrapper's reach
    monkeypatch.setattr(g3, "family_at", real)
    with pytest.raises(ValueError, match="no longer out of the wrapper"):
        g3.verify_equal_radius_structure(15.0, 15.0, s1.M, 1.0)


def test_one_turn_beyond_IS_constructible_at_equal_radii():
    """The R7 defect: a case the preflight table was leaving out.

    It asks for an angle above the band's upper edge. `1.5 * tol` is
    1.5e-08, the wrapper recovers 1.4901e-08 from it, and that is
    above `tol`, so the solver reports `one-turn` -- exactly what the
    spec's name says. The earlier version reached the band SEARCH
    first and gave up, which is a fact about the resolver, not about
    representability."""

    import o4b_g3a as g3a
    import s1_schwarzschild_cost as s1

    spec = next(s for s in g3.G3A_ANGLE_SPECS
                if s[0] == "one-turn-beyond")
    theta = g3.resolve_theta(spec, 15.0, 15.0, s1.M, s1.DEFAULT_TOL,
                             sz.PSI_MAX)
    assert theta == 1.5 * s1.DEFAULT_TOL
    dpsi = g3.wrapper_dpsi(theta)
    assert dpsi == 1.4901161193847656e-08 > s1.DEFAULT_TOL
    assert g3.family_at(15.0, 15.0, dpsi, s1.M,
                        s1.DEFAULT_TOL) == "one-turn"

    assert ("equal-radius", "one-turn-beyond") not in (
        g3.EXPECTED_UNREACHABLE)
    table = g3a.run_g3a()
    built = {c["case"]: c["family"] for c in table["detail"]}
    assert built["equal-radius/one-turn-beyond"] == "one-turn"
    assert table["cases"] == 71
    assert len(table["unreachable"]) == 6


def test_the_duplicate_recovered_angles_are_recorded_not_counted():
    """The case is valid and it runs, but at equal radii it recovers
    the SAME dpsi as `radial-reachable-edge-above`, so it adds no
    solver state. Recorded so a reader cannot count the row as
    independent coverage."""

    import o4b_g3a as g3a

    table = g3a.run_g3a()
    dupes = table["duplicate_recovered_angles"]
    incidental = [d for d in dupes if not d["by_construction"]]
    assert len(incidental) == 1
    assert incidental[0]["cases"] == ["equal-radius/one-turn-beyond",
                                      "equal-radius/"
                                      "radial-reachable-edge-above"]
    assert incidental[0]["dpsi"] == 1.4901161193847656e-08
    assert incidental[0]["family"] == "one-turn"

    # the other seven are the reachability-edge pairs, where agreeing
    # with theta = 0 is what the case is there to demonstrate
    assert len(dupes) - len(incidental) == 7
    assert all(set(d["cases"][0].split("/")) & {"radial"}
               for d in dupes if d["by_construction"])

    assert table["cases"] == 71
    assert table["distinct_solver_states"] == 63
    assert "must not be counted as independent" in (
        table["why_duplicates_are_recorded"])


def test_a_named_omission_that_turns_out_buildable_also_fails(
        monkeypatch):
    """The list is an equality, not a lower bound: if a case the freeze
    recorded as impossible does resolve, the frozen table no longer
    describes the solver and the preflight has to say so."""

    import o4b_g3a as g3a

    monkeypatch.setattr(
        g3, "EXPECTED_UNREACHABLE_LABELS",
        g3.EXPECTED_UNREACHABLE_LABELS | {"anchor-pair/no-turn-mid"})
    preflight = g3a.run_preflight()

    assert preflight["table"]["unexpectedly_reachable"] == [
        "anchor-pair/no-turn-mid"]
    assert preflight["conditions"]["only_expected_unreachable"] is False
    assert preflight["passed"] is False


def test_distance_and_work_are_separate_fields():
    """`ulp_distance` is how far the nominal placement sat from a
    satisfying one; it is not a count of anything performed."""

    import o4b_g3a as g3a
    placed = g3a.place_row(0.0026643057799644846,
                           0.0026640756792607376, g3.ETA, above=False)
    assert placed["ulp_distance"] == 8_042
    assert placed["search_comparisons"] < 40
    assert "nudges" not in placed
    assert "거리가 실행 횟수로 오독된다" in _prose()
