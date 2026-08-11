"""PR-O1 contract tests for the certified flight-time enclosure
(docs/theory/schwarzschild_volume_oracle_certification.md L0-L5).

The load-bearing checks are CONTAINMENT checks against references
that are themselves directed-rounding enclosures at higher precision
(200-bit), never bare Python floats: a float reference can sit
outside a correct 96-bit certified interval by more than the
interval's width, so float equality proves nothing either way. The
S1 heuristic comparison is a consistency diagnostic, not a
certification ground."""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import gmpy2
import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import certified_flight_time as cft  # noqa: E402
import certified_interval as ci  # noqa: E402
from certified_flight_time import (  # noqa: E402
    FROZEN_ANCHORS,
    containment_certificate,
    flight_time_certified,
)
from certified_interval import CertificationError, Iv, iv_sum  # noqa: E402
from s1_schwarzschild_cost import flight_time as s1_flight_time  # noqa: E402

_DN200 = gmpy2.context(precision=200, round=gmpy2.RoundDown)
_UP200 = gmpy2.context(precision=200, round=gmpy2.RoundUp)


def _ref_chord(r1: float, r2: float, psi: float) -> tuple:
    """200-bit directed enclosure of the flat chord
    sqrt(r1^2 + r2^2 - 2 r1 r2 cos psi)."""

    lo = _DN200.sqrt(
        _DN200.add(
            _DN200.add(_DN200.mul(gmpy2.mpfr(r1), gmpy2.mpfr(r1)),
                       _DN200.mul(gmpy2.mpfr(r2), gmpy2.mpfr(r2))),
            -_UP200.mul(_UP200.mul(gmpy2.mpfr(2 * r1), gmpy2.mpfr(r2)),
                        _UP200.cos(gmpy2.mpfr(psi)))))
    hi = _UP200.sqrt(
        _UP200.add(
            _UP200.add(_UP200.mul(gmpy2.mpfr(r1), gmpy2.mpfr(r1)),
                       _UP200.mul(gmpy2.mpfr(r2), gmpy2.mpfr(r2))),
            -_DN200.mul(_DN200.mul(gmpy2.mpfr(2 * r1), gmpy2.mpfr(r2)),
                        _DN200.cos(gmpy2.mpfr(psi)))))
    return lo, hi


def _ref_tortoise_diff(r1: float, r2: float, m: float) -> tuple:
    """200-bit directed enclosure of rho(r2) - rho(r1), r2 > r1."""

    def rho(ctx, r):
        # every step is monotone increasing in its operand, so a
        # single rounding direction yields a one-sided bound
        arg = ctx.sub(ctx.div(gmpy2.mpfr(r), gmpy2.mpfr(2 * m)),
                      gmpy2.mpfr(1))
        return ctx.add(gmpy2.mpfr(r),
                       ctx.mul(gmpy2.mpfr(2 * m), ctx.log(arg)))

    lo = _DN200.sub(rho(_DN200, r2), rho(_UP200, r1))
    hi = _UP200.sub(rho(_UP200, r2), rho(_DN200, r1))
    return lo, hi


# ---------------------------------------------------------------------
# the interval core
# ---------------------------------------------------------------------


def test_interval_core_encloses_exact_rationals():
    """Directed arithmetic must enclose exact rational results;
    gmpy2 compares mpfr against mpq exactly, so these are proofs,
    not tolerance checks."""

    third = Iv(1) / Iv(3)
    assert third.lo < gmpy2.mpq(1, 3) < third.hi
    tenth_sum = iv_sum([Iv(0.1)] * 10)
    exact = gmpy2.mpq(Fraction(0.1) * 10)
    assert tenth_sum.lo <= exact <= tenth_sum.hi
    for a in (-2.5, -0.3, 0.7, 3.0):
        for b in (-1.1, -0.2, 0.4, 2.0):
            p = Iv(a) * Iv(b)
            e = gmpy2.mpq(Fraction(a) * Fraction(b))
            assert p.lo <= e <= p.hi, (a, b)
    sq = Iv(-2.0, 3.0).sq()
    assert sq.lo <= 0 and sq.hi >= 9


def test_interval_core_is_fail_closed():
    with pytest.raises(CertificationError):
        Iv(1) / Iv(-1.0, 1.0)
    with pytest.raises(CertificationError):
        Iv(-1.0, 2.0).sqrt()
    with pytest.raises(CertificationError):
        Iv(0.0, 2.0).log()
    with pytest.raises(CertificationError):
        Iv(-0.5, 0.5).cos()
    with pytest.raises(CertificationError):
        Iv(2.0, 4.0).sin()
    with pytest.raises(CertificationError):
        Iv(0.0, 1.0).intersect(Iv(2.0, 3.0))
    with pytest.raises(CertificationError):
        Iv(0.0, 1.0).widen(Iv(-1.0, 1.0))
    with pytest.raises(CertificationError):
        Iv(2.0, 1.0)


def test_interval_transcendentals_are_directed():
    pi = Iv.pi()
    assert pi.lo < pi.hi
    assert float(pi.lo) == pytest.approx(3.14159265358979, abs=1e-12)
    c = Iv(0.5).cos()
    assert c.lo < c.hi
    s = Iv(0.5).sin()
    assert s.lo < s.hi
    lg = Iv(5.0).log()
    assert lg.lo < lg.hi


# ---------------------------------------------------------------------
# flat control: the closed-form chord (M = 0)
# ---------------------------------------------------------------------

_FLAT_CASES = [
    (10.0, 20.0, 0.3),      # no-turn, extreme radial aspect
    (12.0, 12.0, 0.5),      # equal radii: always one-turn
    (10.0, 10.0, 2.0),      # patch-corner angle
    (15.0, 10.0, 1.0),      # one-turn, reversed order
    (12.0, 18.0, 0.818),    # the frozen-anchor box scale
    (12.0, 18.0, 1e-6),     # near-radial
    (12.0, 18.0, 0.8410686705679303),  # ~arccos(2/3): A_eq corridor
]


@pytest.mark.parametrize("r1,r2,psi", _FLAT_CASES)
def test_flat_control_contains_the_chord(r1, r2, psi):
    ft = flight_time_certified(r1, r2, psi, m=0.0)
    lo, hi = _ref_chord(r1, r2, psi)
    assert ft.t.lo <= lo and hi <= ft.t.hi, (
        float(ft.t.lo), float(lo), float(hi), float(ft.t.hi))


def test_flat_radial_is_exact_distance():
    ft = flight_time_certified(10.0, 20.0, 0.0, m=0.0)
    assert ft.family == "radial"
    assert ft.t.contains(10)


# ---------------------------------------------------------------------
# radial control: the exact tortoise difference (M = 1)
# ---------------------------------------------------------------------


def test_radial_control_contains_the_tortoise_difference():
    ft = flight_time_certified(10.0, 20.0, 0.0, m=1.0)
    assert ft.family == "radial"
    lo, hi = _ref_tortoise_diff(10.0, 20.0, 1.0)
    assert ft.t.lo <= lo and hi <= ft.t.hi


# ---------------------------------------------------------------------
# Schwarzschild consistency and convergence
# ---------------------------------------------------------------------

_SCHW_CASES = [
    (12.0, 12.0, 0.5),
    (10.0, 10.0, 2.0),
    (15.0, 10.0, 1.0),
    (12.0, 18.0, 0.818),
    (11.5, 18.5, 0.1),
    (12.0, 18.0, 1e-6),
]


@pytest.mark.parametrize("r1,r2,psi", _SCHW_CASES)
def test_s1_heuristic_lands_inside_the_certified_interval(r1, r2, psi):
    """Consistency diagnostic (not a certification ground): the S1
    solver's value, within its own claimed error, must intersect the
    certified enclosure; disjointness would mean one of the two is
    wrong."""

    ft = flight_time_certified(r1, r2, psi, m=1.0)
    t1, e1 = s1_flight_time(r1, r2, psi)
    assert float(ft.t.lo) - e1 - 1e-9 <= t1 <= float(ft.t.hi) + e1 + 1e-9


def test_refinement_narrows_and_nests():
    """Doubling the quadrature must narrow the enclosure, and all
    enclosures of the same true value must mutually intersect
    (L6d's nesting rule at the flight-time level)."""

    r1, r2, psi = 12.0, 18.0, 0.5
    ladder = [flight_time_certified(r1, r2, psi, m=1.0, n_sub=n).t
              for n in (64, 128, 256)]
    for a, b in zip(ladder[:-1], ladder[1:], strict=True):
        assert b.width() < a.width()
        a.intersect(b)  # raises CertificationError if disjoint


def test_certified_widths_are_useful_on_the_anchor_box_scale():
    """Quality floor, not a soundness claim: at the default
    quadrature the relative width on frozen-box-scale calls must be
    below 1% (the L6 budget is allocated per-cell in PR-O2; this
    guards against silent quality regressions)."""

    for r1, r2, psi in [(12.0, 18.0, 0.818), (11.5, 18.5, 0.1),
                        (12.0, 12.0, 0.5)]:
        ft = flight_time_certified(r1, r2, psi, m=1.0)
        assert ft.t.width() / ft.t.mid() < 0.01, (r1, r2, psi)


# ---------------------------------------------------------------------
# fail-closed behavior
# ---------------------------------------------------------------------


def test_patch_preconditions_are_enforced_before_interval_work():
    with pytest.raises(ValueError):
        flight_time_certified(5.0, 18.0, 0.5)
    with pytest.raises(ValueError):
        flight_time_certified(12.0, 25.0, 0.5)
    with pytest.raises(ValueError):
        flight_time_certified(12.0, 18.0, 2.5)
    with pytest.raises(ValueError):
        flight_time_certified(12.0, 18.0, -0.1)
    with pytest.raises(ValueError):
        flight_time_certified(12.0, 18.0, 0.5, m=-1.0)


def test_q_positivity_guard_refuses_a_bad_turning_point():
    m = Iv(1.0)
    bad_u_t = Iv(1.0 / 2.9)  # inside the photon sphere: Q(u_t) < 0
    with pytest.raises(CertificationError):
        cft._arc_enclosure(Iv(0.05), bad_u_t, Iv(6.0), m, 16, False)


# ---------------------------------------------------------------------
# L4: the containment certificate
# ---------------------------------------------------------------------


def test_frozen_anchor_certificate_matches_the_review_numbers():
    """The review-supplied independent values for (12, 18, 8.5):
    margins 1.5600/3.3326/2.9111/1.8923, box r in
    [11.353642, 18.694960], psi_max 0.817913, perihelion floor
    10.417M. Certified intervals must sit on them."""

    cert = containment_certificate(*FROZEN_ANCHORS)
    expected = {"nonempty": 1.559993, "inner_shell": 3.332581,
                "outer_shell": 2.911139, "polar_cap": 1.892305}
    for key, val in expected.items():
        margin = cert["margins"][key]
        assert margin.certainly_gt(Iv(0))
        assert abs(margin.mid() - val) < 1e-5, key
    r_lo, r_hi = cert["r_box"]
    assert abs(r_lo.mid() - 11.353642) < 1e-5
    assert abs(r_hi.mid() - 18.694960) < 1e-5
    assert abs(cert["psi_max"].mid() - 0.817913) < 1e-5
    assert abs(cert["perihelion_floor"].mid() - 10.417379) < 1e-4
    # every solver call inside the box respects the patch cap, and
    # the perihelion floor keeps orbits inside the S1 shell
    assert cert["psi_max"].hi < cft.PSI_MAX
    assert cert["perihelion_floor"].certainly_gt(Iv(cft.R_MIN))


def test_containment_gate_refuses_oversized_diamonds():
    with pytest.raises(CertificationError):
        containment_certificate(12.0, 18.0, 12.0)
    with pytest.raises(CertificationError):
        containment_certificate(12.0, 18.0, 5.0)  # empty diamond


def test_frozen_anchor_calls_stay_certifiable_across_the_box():
    """Representative T1/T2 calls the PR-O2 integrator will make:
    anchor to box points at the box extremes must certify without
    family failures and with sub-percent width."""

    cert = containment_certificate(*FROZEN_ANCHORS)
    psi_hi = float(cert["psi_max"].lo)
    for r_x, psi in [(11.36, psi_hi), (18.69, psi_hi),
                     (15.0, 0.5 * psi_hi), (12.0, psi_hi)]:
        for r_anchor in FROZEN_ANCHORS[:2]:
            ft = flight_time_certified(r_anchor, r_x, psi, m=1.0)
            assert ft.t.certainly_gt(Iv(0))
            assert ft.t.width() / ft.t.mid() < 0.01, (r_anchor, r_x)


def test_interval_module_precision_is_the_documented_contract():
    assert ci.PRECISION == 96
