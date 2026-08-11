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


def big_r_iv(u: Iv, b: Iv, m: Iv) -> Iv:
    return Iv(1) / b.sq() - u.sq() * (Iv(1) - _TWO * m * u)


def q_iv(u: Iv, u_t: Iv, m: Iv) -> Iv:
    """Q(u) = (u + u_t) - 2M (u^2 + u u_t + u_t^2), the regular
    cofactor of R(u) = (u_t - u) Q(u) (L2b)."""

    return (u + u_t) - _TWO * m * (u.sq() + u * u_t + u_t.sq())


@dataclass(frozen=True)
class CertifiedFlightTime:
    """A certified enclosure of the optical distance T, with the
    matched family and diagnostics."""

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


def _arc_enclosure(u_from: Iv, u_t: Iv, b: Iv, m: Iv,
                   n_sub: int, want_time: bool) -> tuple[Iv, Iv]:
    """(angle, time) enclosure of one arc from u_from to the turning
    point u_t, over s = sqrt(u_t - u): interval Riemann sum of
    2/sqrt(Q) on [0, s_hi.lo] plus a [0, +] tail covering the
    uncertain top of the range (u_t is itself an interval)."""

    du = u_t - u_from
    if du.hi < 0:
        raise CertificationError(f"arc with u_from above u_t: {du}")
    s_hi = du.sqrt0()
    ang_terms, tim_terms = [], []

    def push(s0: float, s1: float, tail: bool) -> None:
        s_box = Iv(s0, s1)
        u_box = u_t - s_box.sq()
        q = q_iv(u_box, u_t, m)
        g = _TWO / q.sqrt()
        ds = Iv(s1) - Iv(s0)
        a = ds * g
        ang_terms.append(a.hull(Iv(0)) if tail else a)
        if want_time:
            t = a / (b * u_box.sq() * (Iv(1) - _TWO * m * u_box))
            tim_terms.append(t.hull(Iv(0)) if tail else t)

    body_hi = float(s_hi.lo)
    if body_hi > 0.0:
        pts = _grid(0.0, body_hi, n_sub)
        for s0, s1 in zip(pts[:-1], pts[1:], strict=True):
            push(s0, s1, tail=False)
    if float(s_hi.hi) > body_hi:
        # the true upper limit lies in [s_hi.lo, s_hi.hi]; the tail
        # contribution is between 0 and the full subinterval bound
        push(body_hi, float(s_hi.hi), tail=True)
    ang = iv_sum(ang_terms) if ang_terms else Iv(0)
    tim = iv_sum(tim_terms) if tim_terms else Iv(0)
    return ang, tim


def _direct_enclosure(u_lo: Iv, u_hi: Iv, b: Iv, m: Iv,
                      n_sub: int, want_time: bool) -> tuple[Iv, Iv]:
    """(angle, time) enclosure of a monotone sub-critical arc by
    direct interval Riemann integration of 1/sqrt(R) in u (no
    turning point exists: L2d certified)."""

    ang_terms, tim_terms = [], []
    pts = _grid(float(u_lo.lo), float(u_hi.hi), n_sub)
    for x0, x1 in zip(pts[:-1], pts[1:], strict=True):
        u_box = Iv(x0, x1)
        r_val = big_r_iv(u_box, b, m)
        g = Iv(1) / r_val.sqrt()
        ds = Iv(x1) - Iv(x0)
        a = ds * g
        ang_terms.append(a)
        if want_time:
            tim_terms.append(
                a / (b * u_box.sq() * (Iv(1) - _TWO * m * u_box)))
    # the grid range [u_lo.lo, u_hi.hi] COVERS the true integration
    # range [u_lo*, u_hi*] (positive integrand: upper bound), while
    # the lower bound may only count cells PROVEN inside it -- cells
    # with x0 >= u_lo.hi and x1 <= u_hi.lo
    ang_full = iv_sum(ang_terms)
    tim_full = iv_sum(tim_terms) if want_time else Iv(0)
    inside = [i for i, (x0, x1) in
              enumerate(zip(pts[:-1], pts[1:], strict=True))
              if x0 >= u_lo.hi and x1 <= u_hi.lo]
    if len(inside) == len(ang_terms):
        return ang_full, tim_full
    core_a = iv_sum(ang_terms[i] for i in inside)
    ang = Iv(min(core_a.lo, ang_full.lo), ang_full.hi)
    if not want_time:
        return ang, Iv(0)
    core_t = iv_sum(tim_terms[i] for i in inside)
    return ang, Iv(min(core_t.lo, tim_full.lo), tim_full.hi)


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
    lo, hi = 0.0, float(u_cap.lo)
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
    if u_t is not None:
        a_a, t_a = _arc_enclosure(u_a, u_t, b, m, n_sub, want_time)
        a_b, t_b = _arc_enclosure(u_b, u_t, b, m, n_sub, want_time)
        if one_turn:
            return a_a + a_b, t_a + t_b
        return (a_a - a_b).abs(), (t_a - t_b).abs()
    if one_turn:
        raise CertificationError(
            "one-turn arc requested at sub-critical impact parameter")
    return _direct_enclosure(iv_min(u_a, u_b), iv_max(u_a, u_b),
                             b, m, n_sub, want_time)


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
        t = t_eq.widen(lip * (dpsi_i - a_eq).abs())
        return CertifiedFlightTime(t.intersect(mono_box),
                                   "equal-perihelion", b_eq,
                                   a_eq, n_sub)

    one_turn = dpsi_i.certainly_gt(a_eq)
    if one_turn:
        # L2a: the matched Clairaut constant is at least w at the
        # certified perihelion floor r_min cos(dpsi/2)
        r_floor = iv_min(r1_i, r2_i) * (dpsi_i / _TWO).cos()
        if not r_floor.certainly_gt(Iv(3) * m_i):
            raise CertificationError(
                f"perihelion floor {r_floor} not above photon sphere")
        b_lo = float(w_iv(r_floor, m_i).lo)
        b_hi = float(b_eq.hi)
    else:
        b_lo = float(b_eq.lo) * 1e-8
        b_hi = float(b_eq.hi)

    lo, hi = b_lo, b_hi
    n_scan = max(16, n_sub // 4)
    n_scan_cap = 8 * n_sub
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
    b_hat = 0.5 * (lo + hi)
    ang, tim = _swept_enclosure(u1, u2, Iv(b_hat), m_i, one_turn,
                                max(n_sub, n_scan), want_time=True)
    t = tim.widen(lip * (dpsi_i - ang).abs())
    return CertifiedFlightTime(
        t.intersect(mono_box), "one-turn" if one_turn else "no-turn",
        Iv(lo, hi), ang, n_sub)


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
    wg = w_glob_iv(m_i) if m > 0.0 else None

    margins = {
        "nonempty": dt_i - (rho_out - rho_in),
        "inner_shell": (rho_in - rho_smin) + (rho_out - rho_smin)
        - dt_i,
        "outer_shell": (rho_smax - rho_in) + (rho_smax - rho_out)
        - dt_i,
    }
    if wg is not None:
        margins["polar_cap"] = (_TWO * wg * Iv(cap_half_angle)
                                - dt_i)
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
    out = {"margins": margins, "r_box": (r_lo, r_hi)}
    if wg is not None:
        psi_max = dt_i / (_TWO * wg)
        out["psi_max"] = psi_max
        out["perihelion_floor"] = r_lo * (psi_max / _TWO).cos()
    return out


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
