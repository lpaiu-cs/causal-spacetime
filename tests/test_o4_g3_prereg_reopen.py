"""The prereg re-opening: the eight things it freezes, and the code
that must agree with each.

The document is the artifact. These tests keep its decisions from being
lost in a later edit, and tie every number in it to the constant that
will actually be used -- including the scan cap, which is read from the
campaign's own sizing rather than transcribed.

Nothing here asserts a measured value. O4 has no verdict.
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
import o4_sizing as sz  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402

_DOC = _REPO / "docs" / "prereg" / "p14_o4_g3_prereg_reopen.md"
_TOL = s1.DEFAULT_TOL


def _prose() -> str:
    text = _DOC.read_text(encoding="utf-8")
    lines = [ln.lstrip("> ") for ln in text.replace("**", "").split("\n")]
    return " ".join(" ".join(lines).split())


# ------------------------------------------- the eight frozen items

def test_it_freezes_rules_and_publishes_no_result():
    plain = _prose()
    assert "결과 없음" in plain
    assert "무판정(`verdict = null`)" in plain
    for forbidden in ("CONCORDANT", "DISCORDANT", "통과했다"):
        assert forbidden not in plain


def test_the_availability_structure_is_double_and_says_why():
    plain = _prose()
    assert g3.N_AVAIL == 100_000
    assert g3.K_G3B == 100_000
    assert "고정분모 availability 추정량이 아니다" in plain
    assert "분모가 결과에 의해 정해진다" in plain
    assert "고정 prefix 가 availability 를 보고" in plain


def test_the_scan_cap_is_read_from_the_frozen_sizing():
    """Not a transcribed integer: the cap must move if the new freeze
    re-derives its own G1 size."""

    assert g3.scan_cap() == sz.N_G1
    plain = _prose()
    assert "새 임의 숫자를 들이지 않는다" in plain
    assert "o4_g3_redesign.scan_cap()" in plain


def test_the_quoted_cap_discrepancy_is_flagged_not_silently_resolved():
    """The ruling cited 26,012,722, which is the pre-rounding threshold
    in a comment, not the frozen size. Resolving that quietly is
    exactly what this program does not do."""

    plain = _prose()
    assert "26,012,722` 은 동결값이 아니다" in plain
    assert "26,200,000" in plain
    assert sz.N_G1 == 26_200_000
    sizing = (_REPO / "experiments" / "oracle"
              / "o4_sizing.py").read_text(encoding="utf-8")
    # it appears in a comment only -- never as a bound
    assert "26,012,722" in sizing
    assert "26_012_722" not in sizing


def test_a_shortfall_is_inconclusive_not_invalid():
    plain = _prose()
    assert "`INCONCLUSIVE`" in plain
    assert "기기 위반을 뜻하는 `INVALID` 로도, mismatch 로도 처리하지 않는다" in (
        plain)


def test_fully_testable_has_all_three_conditions():
    plain = _prose()
    assert "W_robust = L − err₁ − err₂ − 2η > 0" in plain
    assert "하단 probe 가 `t_x ≥ 0`" in plain
    assert "세 probe 모두" in plain and "실현 여유 `≥ η` 달성" in plain


def test_mismatch_zero_is_the_gate_and_cp_is_only_reported():
    plain = _prose()
    assert "mismatch 0" in plain
    assert "characterization 으로 보고할 뿐 별도 accuracy gate 로 쓰지 않는다" in (
        plain)
    assert "coverage · sidedness · power" in plain


def test_the_cp_bound_in_the_document_is_the_one_the_code_computes():
    got = sz.g3_upper(g3.K_G3B, g3.ALPHA_G3B)
    assert got == pytest.approx(2.995687401941005e-05, rel=1e-12)
    assert repr(got) in _prose() or "2.995687401941005e-05" in _prose()


def test_eta_and_max_nudges_are_untouched_by_the_reopening():
    assert g3.ETA == 1e-12
    assert g3.MAX_NUDGES == 64
    plain = _prose()
    assert "census 도 이 재개방도 바꾸지 않는다" in plain
    # the value was not re-selected; only its justification changed
    assert "충분성 논증은 구현 중 철회" in plain
    assert "값이 재선택된 것이 아니라" in plain


def test_every_abort_path_leaves_a_write_once_incident():
    plain = _prose()
    assert "모든 abort 경로" in plain
    for path in ("`INVALID`(G3a 실패)", "`INCONCLUSIVE`(스캔 소진)",
                 "fail-closed 솔버 경로", "상한 도달"):
        assert path in plain


def test_the_census_hundred_percent_is_not_reused_as_evidence():
    plain = _prose()
    assert "비용 예상과 설계 입력" in plain
    assert "새 캠페인의 통과 증거로 재사용되지 않는다" in plain


# ------------------------------------------------ the G3a case table

def test_the_case_table_is_frozen_in_code_and_named_in_the_document():
    plain = _prose()
    assert "G3A_GEOMETRIES" in plain and "G3A_ANGLE_SPECS" in plain
    assert len(g3.G3A_GEOMETRIES) == 7
    assert len(g3.G3A_ANGLE_SPECS) == 11
    assert g3.G3A_ROWS == ("A", "B", "C")
    assert set(g3.FAMILIES) == {"radial", "no-turn",
                                "equal-perihelion", "one-turn"}


def test_the_wrapper_cannot_reach_the_radial_threshold():
    """The finding that moved the table into theta space: the wrapper
    recovers `dpsi` as `acos(clamp(cos theta))`, whose smallest nonzero
    value is 1.49e-08 -- so the solver's own 1e-12 threshold sits in a
    gap no coordinate can produce."""

    smallest = math.acos(math.nextafter(1.0, -math.inf))
    assert smallest == pytest.approx(1.4901161193847656e-08, rel=1e-12)
    assert g3.RADIAL_THRESHOLD < smallest
    reachable = {g3.wrapper_dpsi(t) for t in
                 (0.0, 1e-12, 1e-10, 1e-9, 1e-8, 2e-8, 1e-7)}
    assert not any(0.0 < d < smallest for d in reachable)
    assert "도달 불가 구간 안에 있다" in _prose()


def test_the_radial_edge_the_table_straddles_is_the_reachable_one():
    below = g3.resolve_theta(("x", "radial-edge", "below"), 13.0, 17.0,
                             s1.M, _TOL, sz.PSI_MAX)
    above = g3.resolve_theta(("x", "radial-edge", "above"), 13.0, 17.0,
                             s1.M, _TOL, sz.PSI_MAX)
    assert math.nextafter(below, math.inf) == above   # adjacent thetas
    assert g3.wrapper_dpsi(below) == 0.0
    assert g3.wrapper_dpsi(above) > 0.0
    assert g3.family_at(13.0, 17.0, g3.wrapper_dpsi(below), s1.M,
                        _TOL) == "radial"
    assert g3.family_at(13.0, 17.0, g3.wrapper_dpsi(above), s1.M,
                        _TOL) != "radial"
    assert "radial-reachable-edge-*" in _prose()


def test_the_radial_branch_reports_exactly_zero_error():
    """Which is why row C works there at all: `|dt - t_min| <= 0` only
    at exact equality."""

    _, err = s1.flight_time(13.0, 17.0, 0.0, s1.M, _TOL)
    assert err == 0.0


def test_the_band_is_located_by_asking_the_solver_not_re_deriving():
    lower, upper = g3.equal_perihelion_band(13.0, 17.0, s1.M, _TOL)
    assert lower < upper
    for x in (lower, upper, 0.5 * (lower + upper)):
        assert g3.family_at(13.0, 17.0, x, s1.M,
                            _TOL) == "equal-perihelion"
    assert g3.family_at(13.0, 17.0, math.nextafter(lower, -math.inf),
                        s1.M, _TOL) != "equal-perihelion"
    assert g3.family_at(13.0, 17.0, math.nextafter(upper, math.inf),
                        s1.M, _TOL) != "equal-perihelion"
    assert "솔버가 보고하는 `family` 라벨을 이분법으로 좁혀" in _prose()


def test_a_geometry_with_no_band_raises_instead_of_pretending():
    """Equal radii have a zero arc and are one-turn everywhere. That is
    a fact about the geometry, so it surfaces."""

    assert g3.family_at(15.0, 15.0, 0.3, s1.M, _TOL) == "one-turn"
    with pytest.raises(ValueError, match="no no-turn/one-turn"):
        g3.equal_perihelion_band(15.0, 15.0, s1.M, _TOL)


def test_every_case_is_a_theta_the_wrapper_can_be_handed():
    """The R1 defect: a case frozen at an unreachable `dpsi` would
    pre-compute on one branch while the wrapper took another."""

    import numpy as np
    for spec in g3.G3A_ANGLE_SPECS:
        theta = g3.resolve_theta(spec, 13.0, 17.0, s1.M, _TOL,
                                 sz.PSI_MAX)
        p = np.array([0.0, 13.0, 0.0, 0.0])
        x = np.array([1.0, 17.0, theta, 0.0])
        cosang = (math.sin(p[2]) * math.sin(x[2])
                  * math.cos(p[3] - x[3])
                  + math.cos(p[2]) * math.cos(x[2]))
        assert g3.wrapper_dpsi(theta) == math.acos(
            max(-1.0, min(1.0, cosang)))


def test_the_table_covers_every_family_on_a_representative_geometry():
    seen = set()
    for spec in g3.G3A_ANGLE_SPECS:
        theta = g3.resolve_theta(spec, 13.0, 17.0, s1.M, _TOL,
                                 sz.PSI_MAX)
        seen.add(g3.family_at(13.0, 17.0, g3.wrapper_dpsi(theta),
                              s1.M, _TOL))
    assert seen == set(g3.FAMILIES)


def test_the_band_edges_are_straddled_by_adjacent_thetas():
    for which in ("lower", "upper"):
        below = g3.resolve_theta(("x", "band-edge", f"{which}-below"),
                                 13.0, 17.0, s1.M, _TOL, sz.PSI_MAX)
        above = g3.resolve_theta(("x", "band-edge", f"{which}-above"),
                                 13.0, 17.0, s1.M, _TOL, sz.PSI_MAX)
        assert math.nextafter(below, math.inf) == above
        fams = [g3.family_at(13.0, 17.0, g3.wrapper_dpsi(t), s1.M,
                             _TOL) for t in (below, above)]
        assert ("equal-perihelion" in fams) and fams[0] != fams[1]
    assert "인접한 binary64 `θ` 쌍" in _prose()


def test_every_frozen_geometry_resolves_and_covers_what_it_can():
    """Each geometry either reaches all four families, or says which
    it cannot reach and why. Equal radii have a zero arc, so no-turn
    and equal-perihelion do not exist there -- that is geometry, not a
    gap in the table."""

    radii = {"R_LO": sz.R_LO, "R_HI": sz.R_HI}
    for name, r1, r2 in g3.G3A_GEOMETRIES:
        if r1 is None:
            r1, r2 = sz.R_IN, sz.R_OUT
        r1 = radii.get(r1, r1)
        r2 = radii.get(r2, r2)
        families, bandless = set(), False
        for spec in g3.G3A_ANGLE_SPECS:
            try:
                theta = g3.resolve_theta(spec, r1, r2, s1.M, _TOL,
                                         sz.PSI_MAX)
            except ValueError:
                bandless = True
                continue
            families.add(g3.family_at(r1, r2, g3.wrapper_dpsi(theta),
                                      s1.M, _TOL))
        if name == "equal-radius":
            assert bandless is True
            assert families == {"radial", "one-turn"}
        else:
            assert bandless is False, name
            assert families == set(g3.FAMILIES), name


def test_row_b_must_not_collapse_into_the_negative_dt_row():
    plain = _prose()
    assert "B 는 D 와 겹치면 안 된다" in plain
    assert "construction-unavailable" in plain


def test_the_negative_dt_row_must_prove_the_solver_was_not_called():
    """A returned `False` does not show the solver was skipped; only
    making the call fail does."""

    import numpy as np

    def explode(*args, **kwargs):
        raise AssertionError("flight_time must not be called")

    original = s1.flight_time
    s1.flight_time = explode
    try:
        p = np.array([5.0, 13.0, 0.1, 0.0])
        q = np.array([1.0, 17.0, 0.2, 0.0])       # dt = -4 < 0
        assert s1.causal_relation(p, q, s1.M, _TOL) is False
    finally:
        s1.flight_time = original
    assert "반환값만 보는 검사는" in _prose()


def test_bit_identity_needs_the_recovery_not_just_determinism():
    """Determinism alone does not establish that the wrapper used OUR
    arguments -- the R1 defect in one sentence."""

    plain = _prose()
    assert "결정론만 확인하는 검사는 이 조건을 확인하지 못한다" in plain
    assert "wrapper_dpsi(θ) = acos(clamp(cos θ))` 에서" in plain

    theta = 0.3
    args = (13.0, 17.0, g3.wrapper_dpsi(theta), s1.M, _TOL)
    assert s1.flight_time(*args) == s1.flight_time(*args)
    # and the recovery is not the identity on theta, which is the point
    assert g3.wrapper_dpsi(theta) != theta


def test_the_error_jump_at_the_band_edge_is_recorded_as_measured():
    """The structural reason a fixed margin cannot work: `err` moves
    four decades across one ulp of `dpsi`."""

    lower, upper = g3.equal_perihelion_band(13.0, 17.0, s1.M, _TOL)
    _, inside = s1.flight_time(13.0, 17.0, 0.5 * (lower + upper),
                               s1.M, _TOL)
    _, below = s1.flight_time(13.0, 17.0,
                              math.nextafter(lower, -math.inf),
                              s1.M, _TOL)
    _, above = s1.flight_time(13.0, 17.0,
                              math.nextafter(upper, math.inf),
                              s1.M, _TOL)
    assert below > 1e-3 and above > 1e-3
    assert below / max(inside, 1e-12) > 1e3
    plain = _prose()
    assert "한 ulp 사이에 4자릿수 점프" in plain
    assert "가설이지 확인된 사실이 아니고" in plain


def test_g3a_failure_stops_before_the_seeds():
    plain = _prose()
    assert "fresh 시드 접촉 전 `INVALID` 이며, G3b 는 실행하지 않는다" in plain
    assert "incident 형식" in plain
