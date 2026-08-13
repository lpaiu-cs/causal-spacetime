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

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

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
    assert rd.ETA == 1e-12
    assert "η = 10⁻¹² (frozen)" in _prose()


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
        {"lo": 1.0, "err1": 0.0, "L_minus_errs": 0.0})}
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
