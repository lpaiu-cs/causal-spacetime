"""G3b's inside probe places a time in an INTERVAL.

The two outside probes move one way and the margin grows without
bound, so the monotone search is right for them. The inside probe
carries two legs on one time, and raising it widens one margin while
narrowing the other. These tests hold that difference: that the
success region really is bounded on both sides, that a search assuming
a ray fails on it, and that the verdict comes from recomputing both
conditions rather than from the search that located the ends.
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
import o4b_g3a as g3a  # noqa: E402
import o4b_g3b as g3b  # noqa: E402

ETA = g3.ETA

#: A cluster with room to spare: the window is 30 and the two legs
#: need 22, so the satisfying interval is wide.
ROOMY = dict(t_p=100.0, t_q=130.0,
             t_min_1=10.0, err_1=1e-9,
             t_min_2=12.0, err_2=2e-9, eta=ETA)


def test_a_roomy_cluster_places_a_time_that_satisfies_both_legs():
    placed = g3b.place_inside(**ROOMY)
    assert placed["reached"] is True
    assert placed["realized_margin_1"] >= ETA
    assert placed["realized_margin_2"] >= ETA
    assert ROOMY["t_p"] < placed["t_x"] < ROOMY["t_q"]


def test_the_margins_are_recomputed_from_the_formed_differences():
    """Not carried over from the endpoints that located the interval.

    The differences round, so the numbers that decide have to be the
    ones the predicate would form at the time actually chosen."""

    placed = g3b.place_inside(**ROOMY)
    t_x = placed["t_x"]
    assert placed["dt_1"] == t_x - ROOMY["t_p"]
    assert placed["dt_2"] == ROOMY["t_q"] - t_x
    assert placed["realized_margin_1"] == g3b.leg_margin(
        placed["dt_1"], ROOMY["t_min_1"], ROOMY["err_1"])
    assert placed["realized_margin_2"] == g3b.leg_margin(
        placed["dt_2"], ROOMY["t_min_2"], ROOMY["err_2"])
    assert placed["realized_margin"] == min(placed["realized_margin_1"],
                                            placed["realized_margin_2"])


def test_the_success_region_is_bounded_on_BOTH_sides():
    """The property the outside probes do not have. One representable
    step past either end and a leg fails -- which is why an
    unbounded-above search is the wrong instrument."""

    placed = g3b.place_inside(**ROOMY)
    lo, hi = placed["lo"], placed["hi"]

    def both(t_x: float) -> bool:
        return (g3b.leg_margin(t_x - ROOMY["t_p"], ROOMY["t_min_1"],
                               ROOMY["err_1"]) >= ETA
                and g3b.leg_margin(ROOMY["t_q"] - t_x,
                                   ROOMY["t_min_2"],
                                   ROOMY["err_2"]) >= ETA)

    assert both(lo) and both(hi)
    assert not both(math.nextafter(lo, -math.inf))
    assert not both(math.nextafter(hi, math.inf))


def test_a_probe_placed_at_an_endpoint_would_sit_one_step_from_failing():
    """So the chosen time is the interval's bit-space midpoint. An
    endpoint is by construction the first value that passes."""

    placed = g3b.place_inside(**ROOMY)
    assert placed["lo"] < placed["t_x"] < placed["hi"]
    assert "not an endpoint" in g3b.place_inside.__doc__.replace(
        "\n    ", " ")


def test_an_exponential_bracket_steps_clean_over_a_narrow_interval():
    """The concrete reason the inside probe may not use the monotone
    search. The bracket doubles until the predicate holds; if the
    satisfying set is a few steps wide and sits between two of those
    doublings, every probe misses and the search reports that nothing
    exists."""

    seed = 100.0
    base = g3a._bits(seed)
    lo, hi = base + 100, base + 103          # four representable times

    def both(t_x: float) -> bool:
        return lo <= g3a._bits(t_x) <= hi

    found = g3b._bisect_endpoint(seed, both, upward=True)
    assert found["reached"] is False         # 64 too low, 128 too high
    # the interval was there the whole time
    assert both(g3a._from_bits(lo + 1))


def test_a_cluster_whose_legs_do_not_overlap_is_reported_not_papered():
    """Emptiness is a fact about the cluster. It is what `W_robust > 0`
    is supposed to exclude in advance, so reaching here empty means the
    eligibility predicate and the construction disagree."""

    tight = dict(ROOMY, t_q=121.0)           # window 21, legs need 22
    placed = g3b.place_inside(**tight)
    assert placed["reached"] is False
    assert "do not overlap" in placed["why"]
    assert placed["range"]["non_empty"] is False


def test_the_inside_probe_does_not_call_the_monotone_search():
    """Stated in `place_row`'s docstring and enforced here: the
    docstring is a claim about the code, and this is the check."""

    calls = {"n": 0}

    def forbidden(*args, **kwargs):
        calls["n"] += 1
        raise AssertionError(
            "the inside probe called place_row, whose exponential "
            "bracket is unsound on an interval")

    original = g3a.place_row
    g3a.place_row = forbidden
    try:
        assert g3b.place_inside(**ROOMY)["reached"] is True
    finally:
        g3a.place_row = original
    assert calls["n"] == 0


def test_the_two_endpoints_are_bisected_separately_and_cheaply():
    """Each leg alone IS monotone in `t_x`, so each end is a bisection.
    It is only their conjunction that is not, which is why there are
    two searches and an intersection rather than one search."""

    placed = g3b.place_inside(**ROOMY)
    assert placed["search_comparisons"] < 200
    assert placed["ulp_width"] > 0


def test_the_free_variable_is_t_x_and_the_module_says_why():
    """`place_row` searches in `dt`. Here each leg's `dt` is a
    difference formed from `t_x`, and each subtraction rounds, so a
    satisfying `dt` need not be reachable by any representable `t_x`."""

    assert "THE FREE VARIABLE IS `t_x`, NOT `dt`" in g3b.__doc__
    assert "wrong coordinate" in g3b.__doc__


def test_the_document_freezes_the_execution_order_with_g3b_first():
    doc = (_REPO / "docs" / "prereg"
           / "p14_o4_g3_prereg_reopen.md").read_text(encoding="utf-8")
    assert "G3a → G3b(고정 prefix) → 나머지 G1 → G2" in doc
    assert "G3b 는 계측 전제조건이다" in doc


def test_the_document_freezes_the_prefix_being_unbiased():
    doc = (_REPO / "docs" / "prereg"
           / "p14_o4_g3_prereg_reopen.md").read_text(encoding="utf-8")
    assert ("prefix 에서 스캔된 모든 점은 적격 여부와 무관하게 G1 "
            "누적통계에\n> 들어간다" in doc)
    assert "`Z = 0` 인 점" in doc
    assert "availability\n   보고의 분모" in doc


def test_the_document_states_what_survives_each_failure():
    doc = (_REPO / "docs" / "prereg"
           / "p14_o4_g3_prereg_reopen.md").read_text(encoding="utf-8")
    assert "### 5.1 실패 지점별로 무엇이 남는가" in doc
    # the claim itself is what "no seed was spent" now rests on: the
    # ref is the authority, so not constructing the generator was
    # necessary and never sufficient
    assert "예약도 청구되지 않는다" in doc
    assert "reservation_claimed: false" in doc
    assert "seeds_spent: false" in doc
    assert "그때까지 누적된 G1 prefix 통계" in doc


@pytest.mark.parametrize("eta", [0.0, 1e-15, g3.ETA, 1e-9])
def test_a_larger_eta_can_only_narrow_the_interval(eta):
    """Both legs demand more, from opposite directions. The interval
    shrinks from both ends; it never moves."""

    placed = g3b.place_inside(**dict(ROOMY, eta=eta))
    assert placed["reached"] is True
    assert placed["realized_margin"] >= eta
