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


def test_max_nudges_is_frozen_now_not_after_the_census():
    """The cap decides whether a cluster yields a valid probe or an
    unavailable one, so choosing it later would let the census pick
    which clusters count -- the same trap eta was kept out of."""

    plain = _prose()
    assert g3.MAX_NUDGES == 64
    assert "`MAX_NUDGES = 64` (frozen)" in plain
    assert "census 가 어떤 cluster 를 셀지 고르게 된다" in plain
    assert "사전 동결 계산 예산" in plain


def test_lower_probe_coverage_is_specified_as_a_joint():
    plain = _prose()
    assert "주변합만으로는 답이 나오지 않는다" in plain
    assert "eligible_and_lower_probe_t_x_negative" in plain
    assert "적격 cluster 에서만" in plain


def test_the_post_census_list_excludes_both_frozen_constants():
    plain = _prose()
    assert "η 와 `MAX_NUDGES` 는 여기 없다" in plain


def test_the_nudge_bound_is_a_budget_and_says_it_is_not_a_guarantee():
    """The original derivation assumed a step gains `ulp(1.0)`. It is
    worse than that: the margin is formed at the magnitude of `err`, so
    a one-ulp move of `dt` moves the margin by less again. The
    retraction has to survive in both the constant and the document."""

    plain = _prose()
    assert "철회한다" in plain
    assert "유한한 예산으로 모든 부족분이 닫힌다는 보장은 없다" in plain
    assert "사전 동결 계산 예산" in plain
    # the phrase survives ONLY inside the retraction that quotes it
    assert ("약 17 걸음이면 충분하고 64 는 3.7 배 여유\" 라고 "
            "적었다. 이 논증은 철회한다") in plain

    source = (_REPO / "experiments" / "oracle"
              / "o4_g3_redesign.py").read_text(encoding="utf-8")
    assert "RETRACTED" in source
    assert "No finite budget closes every shortfall" in source


def test_the_quoted_step_count_is_recomputed_not_transcribed():
    """8,042 is asserted in three places. Run the loop and check it --
    an earlier draft quoted 151 from a linear estimate that ignored
    where the margin is formed, and it was wrong by 53x."""

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
    assert steps > g3.MAX_NUDGES
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
