"""Certified minimal null flight time on the S1 Schwarzschild patch.

This is the PR-O1 contract of the oracle certification
(docs/theory/schwarzschild_volume_oracle_certification.md): every
returned time is a directed-rounding MPFR interval PROVEN to contain
the true optical distance, by L0-L5 of that document. The S1
heuristic solver (s1_schwarzschild_cost.py) is untouched and remains
non-certified; nothing here imports it.

Structure of a certified call, mirroring the S1 branches:

- radial pairs: exact tortoise difference (closed form, interval).
- otherwise the equal-perihelion arc A_eq classifies the family;
  the matched Clairaut constant b is bisected with certified
  three-valued angle comparisons (an undecidable step STOPS the
  refinement -- sound, merely wider); arc integrals are interval
  Riemann enclosures over the factored regular integrand 2/sqrt(Q)
  (L2b guarantees Q > 0; every subinterval re-certifies it); the
  residual angle mismatch converts to time error through the L1
  Lipschitz rule |dT/dpsi| <= min(w(r1), w(r2)) -- soundness never
  depends on the bisection having converged.

Fail-closed: patch-precondition violations raise ValueError before
any interval work; certification failures raise CertificationError.
"""

from __future__ import annotations

from dataclasses import dataclass

from certified_interval import (
    CertificationError,
    Iv,
    iv_max,
    iv_min,
    iv_sum,
)

# Frozen S1 patch (the domain L2/L4 certify)
M_DEFAULT = 1.0
R_MIN, R_MAX = 10.0, 20.0
PSI_MAX = 2.0

# Review-approved frozen anchor configuration (certification doc L4):
# radially aligned anchors at 12M and 18M, coordinate-time separation
# 8.5M. PR-O2/O3 freeze and execute on this triple.
FROZEN_ANCHORS = (12.0, 18.0, 8.5)

_TWO = Iv(2)


def f_iv(r: Iv, m: Iv) -> Iv:
    return Iv(1) - _TWO * m / r


def tortoise_iv(r: Iv, m: Iv) -> Iv:
    """rho = r + 2M ln(r/2M - 1); the m = 0 limit is rho = r."""

    if m.lo == 0 and m.hi == 0:
        return r
    return r + _TWO * m * (r / (_TWO * m) - Iv(1)).log()


def w_iv(r: Iv, m: Iv) -> Iv:
    """Optical radius w = r / sqrt(1 - 2M/r)."""

    return r / f_iv(r, m).sqrt()


def w_glob_iv(m: Iv) -> Iv:
    """Global exterior minimum of w: 3 sqrt(3) M at the photon
    sphere (L1)."""

    return Iv(27).sqrt() * m


def angle_cost_iv(m: Iv) -> Iv:
    """Certified constant c with d_opt >= c * dpsi between shell
    points: 3 sqrt(3) M for M > 0 (w's global exterior minimum, L4).
    For M = 0 there is NO photon-sphere floor (w = r has infimum 0),
    so the flat constant comes from the chord instead:
    chord >= 2 sqrt(r1 r2) sin(psi/2) >= (2 R_MIN / pi) psi on the
    shell (sin x >= 2x/pi on [0, pi/2])."""

    if m.lo == 0 and m.hi == 0:
        return Iv(2) * Iv(R_MIN) / Iv.pi()
    return w_glob_iv(m)


def big_r_iv(u: Iv, b: Iv, m: Iv) -> Iv:
    return Iv(1) / b.sq() - u.sq() * (Iv(1) - _TWO * m * u)


def q_iv(u: Iv, u_t: Iv, m: Iv) -> Iv:
    """Q(u) = (u + u_t) - 2M (u^2 + u u_t + u_t^2), the regular
    cofactor of R(u) = (u_t - u) Q(u) (L2b)."""

    return (u + u_t) - _TWO * m * (u.sq() + u * u_t + u_t.sq())


@dataclass(frozen=True)
class CertifiedFlightTime:
    """A certified enclosure of the optical distance T.

    `b` is a certified enclosure of the TRUE Clairaut constant
    b* = dT/dpsi at the requested configuration, built from theory
    (L1/L2a: one-turn b* in [w(r_floor), b_eq]; no-turn and the
    equal-perihelion corridor only certify b* in [0, b_eq]) -- NEVER
    from the bisection search bracket, which narrows by comparisons
    whose monotonicity is not certified and so encloses nothing.
    L6a consumers must hull their gradient over this interval."""

    t: Iv
    family: str
    b: Iv | None
    ang: Iv | None
    n_sub: int


def _grid(lo: float, hi: float, n: int) -> list[float]:
    step = (hi - lo) / n
    pts = [lo + k * step for k in range(n)]
    pts.append(hi)
    return pts


def _even_poly(coeffs: list[Iv], s: Iv) -> Iv:
    """c0 + c1 x + c2 x^2 + ... with x = s^2 (even polynomial in s),
    interval Horner."""

    x = s.sq()
    acc = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        acc = acc * x + c
    return acc


def _arc_polys(u_t: Iv, b: Iv, m: Iv,
               ) -> tuple[list[Iv], list[Iv]]:
    """Even-polynomial forms of the factored arc integrands in
    s = sqrt(u_t - u) (exact substitution; at an exact root u_t the
    leading coefficient satisfies d0 = 1/b -- a sanity diagnostic,
    not an enforced identity, since u_t arrives as an interval):

      Q~(s) = Q(u_t - s^2) = c0 + c2 s^2 + c4 s^4,
        c0 = 2 u_t (1 - 3 m u_t), c2 = 6 m u_t - 1, c4 = -2m
      D~(s) = b u^2 (1 - 2 m u)|_{u = u_t - s^2}
            = d0 + d2 s^2 + d4 s^4 + d6 s^6,
        d0 = b u_t^2 (1 - 2 m u_t) = 1/b,  d2 = -2 b u_t (1 - 3 m u_t),
        d4 = b (1 - 6 m u_t),              d6 = 2 m b

    Coefficient lists are in powers of x = s^2."""

    one = Iv(1)
    q_c = [_TWO * u_t * (one - Iv(3) * m * u_t),
           Iv(6) * m * u_t - one,
           -_TWO * m]
    d_c = [b * u_t.sq() * (one - _TWO * m * u_t),
           -_TWO * b * u_t * (one - Iv(3) * m * u_t),
           b * (one - Iv(6) * m * u_t),
           _TWO * m * b]
    return q_c, d_c


def _even_d1(coeffs: list[Iv], s: Iv) -> Iv:
    """d/ds of an even polynomial: sum 2k c_k s^(2k-1)."""

    x = s.sq()
    acc = Iv(0)
    for k in range(len(coeffs) - 1, 0, -1):
        acc = acc * x + Iv(2 * k) * coeffs[k]
    return acc * s


def _even_d2(coeffs: list[Iv], s: Iv) -> Iv:
    """d^2/ds^2 of an even polynomial: sum 2k(2k-1) c_k s^(2k-2)."""

    x = s.sq()
    acc = Iv(0)
    for k in range(len(coeffs) - 1, 0, -1):
        acc = acc * x + Iv(2 * k * (2 * k - 1)) * coeffs[k]
    return acc


def _midpoint_term(f_mid: Iv, f_dd: Iv, s0: float, s1: float) -> Iv:
    """Certified composite-midpoint step: the integral over [s0, s1]
    lies in h f(mid) +- h^3 |f''|_sup / 24."""

    h = Iv(s1) - Iv(s0)
    err = h * h.sq() * Iv(f_dd.abs().hi) / Iv(24)
    return (h * f_mid).widen(err)


def _arc_enclosure(u_from: Iv, u_t: Iv, b: Iv, m: Iv,
                   n_sub: int, want_time: bool) -> tuple[Iv, Iv]:
    """(angle, time) enclosure of one arc from u_from to the turning
    point u_t, over s = sqrt(u_t - u): certified composite midpoint
    rule (second order) on [0, s_hi.lo] with interval bounds on the
    exact even-polynomial second derivatives, plus a first-order
    [0, +] tail covering the uncertain top of the range (u_t is
    itself an interval). Every step re-certifies Q > 0 (L2b)."""

    du = u_t - u_from
    if du.hi < 0:
        raise CertificationError(f"arc with u_from above u_t: {du}")
    s_hi = du.sqrt0()
    q_c, d_c = _arc_polys(u_t, b, m)
    ang_terms, tim_terms = [], []

    def eval_fs(s: Iv) -> tuple[Iv, Iv]:
        """(g, tau) = (2 Q~^-1/2, 2 Q~^-1/2 / D~) at s."""

        q = _even_poly(q_c, s)
        g = _TWO / q.sqrt()
        if not want_time:
            return g, g
        return g, g / _even_poly(d_c, s)

    def eval_dd(s_box: Iv) -> tuple[Iv, Iv]:
        """Interval bounds of g'' and tau'' over s_box, via
        A = Q~^-1/2, B = D~^-1:
          A'  = -1/2 Q^-3/2 Q', A'' = 3/4 Q^-5/2 Q'^2 - 1/2 Q^-3/2 Q''
          B'  = -D^-2 D',       B'' = 2 D^-3 D'^2 - D^-2 D''
          g'' = 2 A'',          tau'' = 2 (A'' B + 2 A' B' + A B'')."""

        q = _even_poly(q_c, s_box)
        qp = _even_d1(q_c, s_box)
        qpp = _even_d2(q_c, s_box)
        sq_ = q.sqrt()
        q32 = q * sq_
        q52 = q.sq() * sq_
        a2 = Iv(3) / Iv(4) * qp.sq() / q52 - qpp / (_TWO * q32)
        g_dd = _TWO * a2
        if not want_time:
            return g_dd, g_dd
        d = _even_poly(d_c, s_box)
        dp = _even_d1(d_c, s_box)
        dpp = _even_d2(d_c, s_box)
        a0 = Iv(1) / sq_
        a1 = -qp / (_TWO * q32)
        b0 = Iv(1) / d
        b1 = -dp / d.sq()
        b2 = _TWO * dp.sq() / (d * d.sq()) - dpp / d.sq()
        t_dd = _TWO * (a2 * b0 + _TWO * a1 * b1 + a0 * b2)
        return g_dd, t_dd

    # directed float conversions: the body top rounds DOWN and the
    # tail top rounds UP, so [0, body_hi] u [body_hi, tail_hi] COVERS
    # the true [0, s_hi*]
    body_hi = s_hi.lo_float()
    tail_hi = s_hi.hi_float()
    if body_hi > 0.0:
        pts = _grid(0.0, body_hi, n_sub)
        for s0, s1 in zip(pts[:-1], pts[1:], strict=True):
            g_dd, t_dd = eval_dd(Iv(s0, s1))
            g_mid, t_mid = eval_fs(Iv(0.5 * (s0 + s1)))
            ang_terms.append(_midpoint_term(g_mid, g_dd, s0, s1))
            if want_time:
                tim_terms.append(_midpoint_term(t_mid, t_dd, s0, s1))
    if tail_hi > body_hi:
        # first-order tail: the true upper limit lies inside, so the
        # contribution is between 0 and the full subinterval bound
        s_box = Iv(body_hi, tail_hi)
        q = _even_poly(q_c, s_box)
        g = _TWO / q.sqrt()
        ds = Iv(tail_hi) - Iv(body_hi)
        ang_terms.append((ds * g).hull(Iv(0)))
        if want_time:
            t = ds * g / _even_poly(d_c, s_box)
            tim_terms.append(t.hull(Iv(0)))
    ang = iv_sum(ang_terms) if ang_terms else Iv(0)
    tim = iv_sum(tim_terms) if tim_terms else Iv(0)
    return ang, tim


def _poly(coeffs: list[Iv], x: Iv) -> Iv:
    acc = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        acc = acc * x + c
    return acc


def _poly_d(coeffs: list[Iv]) -> list[Iv]:
    return [Iv(k) * c for k, c in enumerate(coeffs)][1:]


def _direct_enclosure(u_lo: Iv, u_hi: Iv, b: Iv, m: Iv,
                      n_sub: int, want_time: bool) -> tuple[Iv, Iv]:
    """(angle, time) enclosure of a monotone sub-critical arc by a
    certified composite midpoint rule in u (no turning point exists:
    L2d certified): integrands R^-1/2 and R^-1/2 / D with
    R = 1/b^2 - u^2 + 2 m u^3 and D = b (u^2 - 2 m u^3), second
    derivatives bounded by interval evaluation of the exact
    polynomial derivatives."""

    r_c = [Iv(1) / b.sq(), Iv(0), Iv(-1), _TWO * m]
    d_c = [Iv(0), Iv(0), b, -_TWO * m * b]
    r_c1, d_c1 = _poly_d(r_c), _poly_d(d_c)
    r_c2, d_c2 = _poly_d(r_c1), _poly_d(d_c1)
    ang_terms, tim_terms = [], []
    # endpoint-uncertainty slivers get their OWN [0, +] cells so the
    # two-sided body runs exactly over the PROVEN-inside range
    # [u_lo.hi, u_hi.lo] -- excluding whole grid cells at the
    # boundary would cost an O(1/n) lower-bound deficit that buries
    # the midpoint rule's O(1/n^2)
    a0, a1 = u_lo.lo_float(), u_lo.hi_float()
    b0_f, b1_f = u_hi.lo_float(), u_hi.hi_float()
    if a1 >= b0_f:
        # endpoint uncertainties overlap (near-degenerate range):
        # one [0, +] cell over the whole covering range is sound
        u_box = Iv(a0, b1_f)
        g = Iv(1) / _poly(r_c, u_box).sqrt()
        ds = Iv(b1_f) - Iv(a0)
        ang = (ds * g).hull(Iv(0))
        if not want_time:
            return ang, Iv(0)
        return ang, (ds * g / _poly(d_c, u_box)).hull(Iv(0))
    for s_lo, s_hi in ((a0, a1), (b0_f, b1_f)):
        if s_hi > s_lo:
            u_box = Iv(s_lo, s_hi)
            g = Iv(1) / _poly(r_c, u_box).sqrt()
            ds = Iv(s_hi) - Iv(s_lo)
            ang_terms.append((ds * g).hull(Iv(0)))
            if want_time:
                t = ds * g / _poly(d_c, u_box)
                tim_terms.append(t.hull(Iv(0)))
    pts = _grid(a1, b0_f, n_sub)
    for x0, x1 in zip(pts[:-1], pts[1:], strict=True):
        u_box = Iv(x0, x1)
        r_box = _poly(r_c, u_box)
        sq_ = r_box.sqrt()
        rp, rpp = _poly(r_c1, u_box), _poly(r_c2, u_box)
        r32, r52 = r_box * sq_, r_box.sq() * sq_
        a2 = Iv(3) / Iv(4) * rp.sq() / r52 - rpp / (_TWO * r32)
        u_mid = Iv(0.5 * (x0 + x1))
        g_mid = Iv(1) / _poly(r_c, u_mid).sqrt()
        if want_time:
            d_box = _poly(d_c, u_box)
            dp, dpp = _poly(d_c1, u_box), _poly(d_c2, u_box)
            a0 = Iv(1) / sq_
            a1 = -rp / (_TWO * r32)
            b0 = Iv(1) / d_box
            b1 = -dp / d_box.sq()
            b2 = (_TWO * dp.sq() / (d_box * d_box.sq())
                  - dpp / d_box.sq())
            t_dd = a2 * b0 + _TWO * a1 * b1 + a0 * b2
            t_mid = g_mid / _poly(d_c, u_mid)
            tim_terms.append(_midpoint_term(t_mid, t_dd, x0, x1))
        ang_terms.append(_midpoint_term(g_mid, a2, x0, x1))
    return (iv_sum(ang_terms),
            iv_sum(tim_terms) if want_time else Iv(0))


def _perihelion_bracket(b: Iv, m: Iv, iters: int = 80) -> Iv | None:
    """Certified bracket [u_lo, u_hi] for the smallest positive root
    of R(u), or None when R is certified positive up to the photon
    sphere (sub-critical b: no turning point on the patch).
    Undecidable criticality raises (L2a keeps oracle calls far from
    the critical impact parameter)."""

    if m.lo == 0 and m.hi == 0:
        return Iv(1) / b  # flat: exact turning point u_t = 1/b
    u_cap = Iv(1) / (Iv(3) * m)
    r_cap = big_r_iv(u_cap, b, m)
    if r_cap.certainly_gt(Iv(0)):
        return None
    if not r_cap.certainly_lt(Iv(0)):
        raise CertificationError(
            f"impact parameter {b} not certifiably super/sub-critical")
    # outward-rounded float cap, with the bracket invariant
    # re-certified AT the float point (nearest-rounding float() could
    # land outside the interval the sign was certified on)
    hi = u_cap.hi_float()
    if not big_r_iv(Iv(hi), b, m).certainly_lt(Iv(0)):
        raise CertificationError(
            f"root-bracket top {hi} not certifiably past the root")
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        r_mid = big_r_iv(Iv(mid), b, m)
        if r_mid.certainly_gt(Iv(0)):
            lo = mid
        elif r_mid.certainly_lt(Iv(0)):
            hi = mid
        else:
            break  # sound: the root stays inside [lo, hi]
    return Iv(lo, hi)


def _swept_enclosure(u_a: Iv, u_b: Iv, b: Iv, m: Iv, one_turn: bool,
                     n_sub: int, want_time: bool) -> tuple[Iv, Iv]:
    """(angle, time) enclosure of the connecting arc at Clairaut
    constant b: sum of to-the-turn arcs (one-turn), difference of
    them (no-turn with a real turning point), or direct integration
    (sub-critical, no turning point) -- S1's `_swept`, certified."""

    u_t = _perihelion_bracket(b, m)
    if one_turn:
        if u_t is None:
            raise CertificationError(
                "one-turn arc requested at sub-critical impact "
                "parameter")
        a_a, t_a = _arc_enclosure(u_a, u_t, b, m, n_sub, want_time)
        a_b, t_b = _arc_enclosure(u_b, u_t, b, m, n_sub, want_time)
        return a_a + a_b, t_a + t_b
    # no-turn: two regular representations, picked by where the
    # turning point sits. With u_t well above the range (certified
    # u_t > 2 u_in) the DIRECT form is well-conditioned -- R(u_in)
    # >= (u_t - u_in) Q > 0 with a fat margin -- and the
    # arc-difference form would subtract two near-equal large-domain
    # integrals (the flat near-radial case has u_t = 1/b -> huge:
    # catastrophically wide). With u_t close to the range (b near
    # b_eq) the direct integrand is near-singular at u_in and the
    # s-substituted arc DIFFERENCE is the regular representation.
    u_in = iv_max(u_a, u_b)
    if u_t is None or u_t.certainly_gt(u_in * _TWO):
        return _direct_enclosure(iv_min(u_a, u_b), u_in,
                                 b, m, n_sub, want_time)
    a_a, t_a = _arc_enclosure(u_a, u_t, b, m, n_sub, want_time)
    a_b, t_b = _arc_enclosure(u_b, u_t, b, m, n_sub, want_time)
    return (a_a - a_b).abs(), (t_a - t_b).abs()


def flight_time_certified(r1: float, r2: float, dpsi: float,
                          m: float = M_DEFAULT, n_sub: int = 128,
                          bisect_iters: int = 200,
                          ) -> CertifiedFlightTime:
    """Certified enclosure of T(r1, r2, dpsi) on the S1 patch.

    Preconditions (ValueError, before any interval work): both radii
    in [R_MIN, R_MAX], 0 <= dpsi <= PSI_MAX, m >= 0. Soundness of the
    returned interval does NOT depend on the bisection converging:
    whatever arc the search lands on is a geodesic, hence THE
    minimizer for its own swept angle (L5), and the residual angle
    mismatch is converted to certified time error by the L1 rule."""

    if not (R_MIN <= r1 <= R_MAX and R_MIN <= r2 <= R_MAX):
        raise ValueError(f"radius outside patch shell: {r1}, {r2}")
    if not 0.0 <= dpsi <= PSI_MAX:
        raise ValueError(f"dpsi outside patch cap: {dpsi}")
    if not m >= 0.0:
        raise ValueError(f"negative mass: {m}")
    m_i = Iv(m)
    r1_i, r2_i = Iv(r1), Iv(r2)

    if dpsi == 0.0:
        t = (tortoise_iv(r2_i, m_i) - tortoise_iv(r1_i, m_i)).abs()
        return CertifiedFlightTime(t, "radial", None, None, n_sub)

    u1, u2 = Iv(1) / r1_i, Iv(1) / r2_i
    u_in, u_out = iv_max(u1, u2), iv_min(u1, u2)
    dpsi_i = Iv(dpsi)
    lip = iv_min(w_iv(r1_i, m_i), w_iv(r2_i, m_i))
    # monotone box (L1): T(0) <= T(psi) <= T(0) + Lip psi, since
    # dT/dpsi = b in [0, Lip]; every family result is intersected
    # with it (both enclosures are certified, so an empty
    # intersection would be a genuine certification failure)
    t_rad = (tortoise_iv(r2_i, m_i) - tortoise_iv(r1_i, m_i)).abs()
    mono_box = Iv(t_rad.lo, (t_rad + lip * dpsi_i).hi)

    # equal-perihelion arc: turning point exactly at the inner
    # endpoint, the analytic family boundary (S1's b_eq)
    b_eq = Iv(1) / (u_in.sq() * (Iv(1) - _TWO * m_i * u_in)).sqrt()
    a_eq, t_eq = _arc_enclosure(u_out, u_in, b_eq, m_i, n_sub, True)

    if not (dpsi_i.certainly_gt(a_eq) or dpsi_i.certainly_lt(a_eq)):
        # decidability corridor: the request angle may sit on either
        # family side of the true A_eq, so the certified b* hull is
        # the union of both families' ranges, [0, b_eq] -- returning
        # the razor-thin b_eq here would hand L6a a gradient the true
        # b* can sit outside of (review R1)
        t = t_eq.widen(lip * (dpsi_i - a_eq).abs())
        return CertifiedFlightTime(t.intersect(mono_box),
                                   "equal-perihelion",
                                   Iv(0).hull(b_eq), a_eq, n_sub)

    one_turn = dpsi_i.certainly_gt(a_eq)
    if one_turn:
        # L2a: the matched Clairaut constant is at least w at the
        # certified perihelion floor r_min cos(dpsi/2)
        r_floor = iv_min(r1_i, r2_i) * (dpsi_i / _TWO).cos()
        if not r_floor.certainly_gt(Iv(3) * m_i):
            raise CertificationError(
                f"perihelion floor {r_floor} not above photon sphere")
        # certified b* hull (L1/L2a): b* = w(r_p) with
        # r_floor <= r_p <= r_inner, so w(r_floor) <= b* <= b_eq
        b_cert = Iv(w_iv(r_floor, m_i).lo, b_eq.hi)
        b_lo = b_cert.lo_float()
        b_hi = b_cert.hi_float()
    else:
        # no-turn: ang(b) is strictly increasing (L2e: fixed limits,
        # pointwise-decreasing R), so a bracket whose endpoint angles
        # certifiably straddle dpsi ENCLOSES the unique b*. The upper
        # straddle is the branch condition itself (dpsi < a_eq); the
        # lower is certified by evaluation, with the family hull
        # [0, b_eq] as the sound fallback.
        b_cert = Iv(0).hull(b_eq)
        b_lo = b_eq.lo_float() * 1e-8
        b_hi = b_eq.hi_float()

    lo, hi = b_lo, b_hi
    n_scan = max(16, n_sub // 4)
    n_scan_cap = 8 * n_sub
    straddle = False
    if not one_turn:
        try:
            ang0, _ = _swept_enclosure(u1, u2, Iv(lo), m_i, False,
                                       n_scan, want_time=False)
            straddle = ang0.certainly_lt(dpsi_i)
        except CertificationError:
            straddle = False
    for _ in range(bisect_iters):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi or (hi - lo) <= 1e-13 * mid:
            break
        side = None
        while True:  # decidability retry: double the quadrature
            try:
                ang, _ = _swept_enclosure(u1, u2, Iv(mid), m_i,
                                          one_turn, n_scan,
                                          want_time=False)
            except CertificationError:
                break  # cannot classify this step: sound to stop
            if ang.certainly_gt(dpsi_i):
                side = "gt"
            elif ang.certainly_lt(dpsi_i):
                side = "lt"
            elif n_scan < n_scan_cap:
                n_scan *= 2
                continue
            break
        if side is None:
            break  # undecidable at the quadrature cap: stay sound
        # one-turn: angle grows as b falls; no-turn: grows with b
        if (side == "gt") == one_turn:
            lo = mid
        else:
            hi = mid
    if not one_turn and straddle:
        # L2e: every accepted lo carried a certified ang < dpsi and
        # every hi a certified ang > dpsi (b_eq's via the branch
        # condition dpsi < a_eq), so [lo, hi] encloses b*
        b_cert = Iv(lo, hi).intersect(b_cert)
    b_hat = 0.5 * (lo + hi)
    ang, tim = _swept_enclosure(u1, u2, Iv(b_hat), m_i, one_turn,
                                max(n_sub, n_scan), want_time=True)
    t = tim.widen(lip * (dpsi_i - ang).abs())
    return CertifiedFlightTime(
        t.intersect(mono_box), "one-turn" if one_turn else "no-turn",
        b_cert, ang, n_sub)


# ---------------------------------------------------------------------
# L4: anchor-diamond containment certificate
# ---------------------------------------------------------------------


def containment_certificate(r_in: float, r_out: float, dt: float,
                            m: float = M_DEFAULT,
                            shell: tuple[float, float] = (R_MIN, R_MAX),
                            cap_half_angle: float = 1.0) -> dict:
    """Certified L4 margin table for radially aligned anchors at
    r_in < r_out with coordinate-time separation dt. Raises
    CertificationError unless EVERY margin is certainly positive --
    the oracle's entry gate. Returns the certified margins, the
    containment box, and the perihelion floor for the box."""

    m_i, dt_i = Iv(m), Iv(dt)
    rho_in = tortoise_iv(Iv(r_in), m_i)
    rho_out = tortoise_iv(Iv(r_out), m_i)
    rho_smin = tortoise_iv(Iv(shell[0]), m_i)
    rho_smax = tortoise_iv(Iv(shell[1]), m_i)
    wg = angle_cost_iv(m_i)

    margins = {
        "nonempty": dt_i - (rho_out - rho_in),
        "inner_shell": (rho_in - rho_smin) + (rho_out - rho_smin)
        - dt_i,
        "outer_shell": (rho_smax - rho_in) + (rho_smax - rho_out)
        - dt_i,
        "polar_cap": _TWO * wg * Iv(cap_half_angle) - dt_i,
    }
    bad = [k for k, v in margins.items()
           if not v.certainly_gt(Iv(0))]
    if bad:
        raise CertificationError(
            f"containment margins not certified positive: {bad}")

    # containment box from the same 1-Lipschitz functionals
    rho_lo = (rho_in + rho_out - dt_i) / _TWO
    rho_hi = (rho_in + rho_out + dt_i) / _TWO
    r_lo = _invert_tortoise(rho_lo, m_i, shell)
    r_hi = _invert_tortoise(rho_hi, m_i, shell)
    psi_max = dt_i / (_TWO * wg)
    return {"margins": margins, "r_box": (r_lo, r_hi),
            "psi_max": psi_max,
            "perihelion_floor": r_lo * (psi_max / _TWO).cos()}


def _invert_tortoise(rho_target: Iv, m: Iv,
                     shell: tuple[float, float]) -> Iv:
    """Certified bracket for r with rho(r) = rho_target, by bisection
    on the strictly increasing tortoise map (drho/dr = 1/f > 0)."""

    lo, hi = shell
    if not (tortoise_iv(Iv(lo), m).certainly_lt(rho_target)
            and tortoise_iv(Iv(hi), m).certainly_gt(rho_target)):
        raise CertificationError(
            f"tortoise inversion target {rho_target} not bracketed "
            f"by shell {shell}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        rho_mid = tortoise_iv(Iv(mid), m)
        if rho_mid.certainly_lt(rho_target):
            lo = mid
        elif rho_mid.certainly_gt(rho_target):
            hi = mid
        else:
            break
    return Iv(lo, hi)
