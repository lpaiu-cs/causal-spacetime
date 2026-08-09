"""S1: the cost of causality on a Schwarzschild patch -- a PRICE, not
a verdict (p14_weyl_curvature.md §8.1; scope frozen by the P14
preregistration review: the number closes or holds only THIS
solver-domain-budget path, and even an affordable price leaves the
Schwarzschild volume/campaign path separately unresolved).

**The predicate.** `p prec q` between events in the exterior patch:
`q` is causally after `p` iff `dt >= T_min(x1, x2)`, the minimal
coordinate flight time of a DIRECT null geodesic between the spatial
points. Spherical symmetry reduces the flight to the plane through
both points and the center; the orbit is parameterized by the impact
parameter `b` with, in `u = 1/r`,

    (du/dphi)^2 = R(u) = 1/b^2 - u^2 (1 - 2 M u)
    dt/du       = 1 / (b u^2 (1 - 2 M u) sqrt(R(u)))

The boundary-value solve brackets `b` so the swept angle matches the
endpoint separation `dpsi`: the no-turn arc (monotone `r`) covers
angles up to the equal-perihelion arc `A_eq` (perihelion at the inner
endpoint), the one-turn arc continues beyond it; on the frozen patch
(shell `r/M in [10, 20]`, `dpsi <= 2` by a polar cap) the direct
solution is unique and the perihelion stays above `~5 M`, far from
the photon sphere -- winding branches never enter. Radial pairs use
the exact tortoise-coordinate difference.

**Tolerance discipline** mirrors P1: quadrature carries a Richardson
error estimate, the solve returns `T_min` WITH an error bound, and a
pair whose `|dt - T_min|` lands inside the bound is UNDECIDED --
counted, never silently classified.

**The measurement.** Wall-clock per predicate call on
coordinate-uniform pairs in the patch (a COST benchmark sampling, not
a physics sprinkling -- no volume claim attaches), the scaling of
that cost with tolerance, and the projection onto P11-P13-like
sample sizes (`n` in the low thousands, `O(n^2)` relations per
sample). Wall-clock is host-dependent; the artifact records the
host and the per-pair microseconds it was measured on.

Run:  python experiments/positive_control/s1_schwarzschild_cost.py
"""

from __future__ import annotations

import json
import math
import platform
import sys
import time
from functools import cache
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------
# Frozen S1 domain (the budget path this price covers)
# ---------------------------------------------------------------------

M = 1.0                      # geometric units G = c = 1
R_MIN, R_MAX = 10.0, 20.0    # exterior shell, photon sphere at 3 M
CAP_HALF_ANGLE = 1.0         # polar cap -> pairwise dpsi <= 2.0
T_EXTENT = 40.0              # coordinate-time extent of the patch
PSI_MAX = 2.0 * CAP_HALF_ANGLE

DEFAULT_TOL = 1e-8
TOL_LADDER = (1e-6, 1e-8, 1e-10)
BENCH_SEED = 40_000_101      # cost benchmark only; above every ledger
BENCH_N = 300                # events per benchmark sample
PROJECTION_N = (1_000, 2_000, 4_000)

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_s1_cost.json")


def _f(r: float, m: float) -> float:
    return 1.0 - 2.0 * m / r


def tortoise(r: float, m: float) -> float:
    """r* = r + 2M ln(r/2M - 1), exact radial flight times."""

    if m == 0.0:
        return r
    return r + 2.0 * m * math.log(r / (2.0 * m) - 1.0)


def _big_r(u: float, b: float, m: float) -> float:
    return 1.0 / (b * b) - u * u * (1.0 - 2.0 * m * u)


def _perihelion_u(b: float, m: float) -> float:
    """Smallest positive root of R(u) = 0 -- the direct orbit's
    turning point -- by bisection from u = 0 (where R > 0). The flat
    limit is exact: R(u) = 1/b^2 - u^2 turns at u = 1/b."""

    if m == 0.0:
        return 1.0 / b
    lo, hi = 0.0, 1.0 / (3.0 * m)
    if _big_r(hi, b, m) > 0.0:
        raise ValueError(f"no turning point: b={b} under critical")
    # bisect to adjacent floats and STOP -- this runs at every step
    # of the impact-parameter solve, and iterations past float
    # convergence would inflate the very wall-clock this measurement
    # exists to report (PR #50 review)
    while True:
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            return mid
        if _big_r(mid, b, m) > 0.0:
            lo = mid
        else:
            hi = mid


@cache
def _nodes(n: int) -> tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(n)


def _gauss(fn, a: float, c: float, nodes: int) -> float:
    """Gauss-Legendre over a VECTORIZED integrand (fn maps arrays to
    arrays) -- the per-call price is dominated by these quadratures,
    so they run as numpy dot products, not Python loops."""

    x, w = _nodes(nodes)
    mid, half = 0.5 * (a + c), 0.5 * (c - a)
    return float(half * np.dot(w, fn(mid + half * x)))


def _refine(fn, a: float, c: float) -> tuple[float, float]:
    """Gauss-Legendre with node doubling until the Richardson-style
    difference stabilizes; returns (value, error bound)."""

    prev = _gauss(fn, a, c, 16)
    for nodes in (32, 64, 128, 256):
        cur = _gauss(fn, a, c, nodes)
        err = abs(cur - prev)
        if err <= 1e-13 * max(1.0, abs(cur)):
            return cur, err
        prev = cur
    return cur, err


def _arc_integrals(u_from: float, u_turn: float, b: float, m: float,
                   want_time: bool) -> tuple[float, float, float]:
    """(angle, time, err) along one arc from `u_from` up to the
    turning point, via s = sqrt(u_turn - u) AND the exact cubic
    factorization R(u) = (u_turn - u) Q(u) with

        Q(u) = -(2m u^2 + (2m u_turn - 1)(u + u_turn))

    so the integrand is 2/sqrt(Q) -- no cancellation anywhere on the
    arc. Evaluating R directly near its root loses everything to
    float cancellation (the equal-radius perihelion probe produced
    1e+133 garbage angles before this form). The time integral is
    skipped during the impact-parameter bisection
    (`want_time=False`) -- only the matched arc pays for it."""

    s_hi = math.sqrt(max(u_turn - u_from, 0.0))

    def q_of(s: np.ndarray) -> np.ndarray:
        u = u_turn - s * s
        return -(2.0 * m * u * u
                 + (2.0 * m * u_turn - 1.0) * (u + u_turn))

    def ang(s: np.ndarray) -> np.ndarray:
        return 2.0 / np.sqrt(q_of(s))

    a1, e1 = _refine(ang, 0.0, s_hi)
    if not want_time:
        return a1, 0.0, e1

    def tim(s: np.ndarray) -> np.ndarray:
        u = u_turn - s * s
        return ang(s) / (b * u * u * (1.0 - 2.0 * m * u))

    t1, e2 = _refine(tim, 0.0, s_hi)
    return a1, t1, e1 + e2


def _swept(u_a: float, u_b: float, b: float, m: float,
           one_turn: bool,
           want_time: bool = True) -> tuple[float, float, float]:
    """(angle, time, err) of the connecting arc at impact parameter
    `b`: sum of to-the-turn arcs for the one-turn path, difference
    for the monotone path when a turning point exists, and direct
    smooth integration when it does not (sub-critical `b` has
    `R > 0` everywhere -- the monotone arc needs no perihelion, and
    on the frozen patch the hull stays far from R's near-zero)."""

    if one_turn:
        u_t = _perihelion_u(b, m)
    else:
        try:
            u_t = _perihelion_u(b, m)
        except ValueError:
            u_t = None
    if u_t is not None:
        a_a, t_a, e_a = _arc_integrals(u_a, u_t, b, m, want_time)
        a_b, t_b, e_b = _arc_integrals(u_b, u_t, b, m, want_time)
        if one_turn:
            return a_a + a_b, t_a + t_b, e_a + e_b
        return abs(a_a - a_b), abs(t_a - t_b), e_a + e_b

    def ang(u: np.ndarray) -> np.ndarray:
        return 1.0 / np.sqrt(_big_r(u, b, m))

    lo_u, hi_u = min(u_a, u_b), max(u_a, u_b)
    a1, e1 = _refine(ang, lo_u, hi_u)
    if not want_time:
        return a1, 0.0, e1

    def tim(u: np.ndarray) -> np.ndarray:
        return ang(u) / (b * u * u * (1.0 - 2.0 * m * u))

    t1, e2 = _refine(tim, lo_u, hi_u)
    return a1, t1, e1 + e2


def flight_time(r1: float, r2: float, dpsi: float,
                m: float = M, tol: float = DEFAULT_TOL,
                details: dict | None = None) -> tuple[float, float]:
    """(T_min, err): minimal null coordinate flight time between
    exterior spatial points separated by center angle `dpsi`.

    Radial pairs are the exact tortoise difference. Otherwise the
    impact parameter is bracketed on the no-turn family (angles up
    to the equal-perihelion arc A_eq) or the one-turn family (angles
    beyond A_eq), then bisected until the angle matches to `tol`.
    Pass a dict as `details` to receive the matched family, impact
    parameter, and perihelion radius (patch-safety diagnostics)."""

    if dpsi < 1e-12:
        if details is not None:
            details.update(family="radial", b=0.0, r_perihelion=None)
        return abs(tortoise(r2, m) - tortoise(r1, m)), 0.0
    u_a, u_b = 1.0 / r1, 1.0 / r2
    u_in, u_out = max(u_a, u_b), min(u_a, u_b)
    # b_eq: perihelion exactly at the inner endpoint. A_eq is the
    # single outer arc with u_turn = u_in passed ANALYTICALLY --
    # root-finding a perturbed b here sits exactly on the degenerate
    # point and the factorized integrand is what keeps it finite
    # (zero arc for equal radii, so equal-radius pairs are always
    # one-turn).
    b_eq = 1.0 / math.sqrt(u_in * u_in * (1.0 - 2.0 * m * u_in))
    a_eq, t_eq, e_eq = _arc_integrals(u_out, u_in, b_eq, m, True)
    if abs(dpsi - a_eq) <= tol:
        if details is not None:
            details.update(family="equal-perihelion", b=b_eq,
                           r_perihelion=1.0 / u_in)
        # dT/dpsi = b exactly, so the skipped angle mismatch is a
        # TIME error of b_eq * |dpsi - a_eq| -- on this patch b_eq
        # is >> 1 and adding bare `tol` under-reported the bound by
        # over an order, letting dt values inside the true band be
        # silently classified (PR #50 review)
        return t_eq, e_eq + abs(dpsi - a_eq) * b_eq
    one_turn = dpsi > a_eq
    if one_turn:
        lo, hi = 3.0 * math.sqrt(3.0) * m * (1.0 + 1e-9), b_eq
        # angle grows as b falls toward the photon sphere
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            ang, _, _ = _swept(u_a, u_b, mid, m, True,
                               want_time=False)
            if ang > dpsi:
                lo = mid
            else:
                hi = mid
            if hi - lo <= tol * mid:
                break
        b_star = 0.5 * (lo + hi)
        ang, t, err = _swept(u_a, u_b, b_star, m, True)
    else:
        lo, hi = 1e-9, b_eq
        # angle grows with b on the monotone family
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            ang, _, _ = _swept(u_a, u_b, mid, m, False,
                               want_time=False)
            if ang < dpsi:
                lo = mid
            else:
                hi = mid
            if hi - lo <= tol * max(mid, 1e-12):
                break
        b_star = 0.5 * (lo + hi)
        ang, t, err = _swept(u_a, u_b, b_star, m, False)
    if details is not None:
        try:
            r_p = 1.0 / _perihelion_u(b_star, m)
        except ValueError:
            r_p = None  # sub-critical monotone arc: no perihelion
        details.update(family="one-turn" if one_turn else "no-turn",
                       b=b_star, r_perihelion=r_p)
    # angle mismatch converts to a time uncertainty via dT/dpsi = b
    return t, err + abs(ang - dpsi) * b_star + tol


def causal_relation(p: np.ndarray, q: np.ndarray,
                    m: float = M, tol: float = DEFAULT_TOL):
    """None if undecided (|dt - T_min| inside the error bound), else
    the boolean `p prec q` for events (t, r, theta, phi)."""

    dt = q[0] - p[0]
    if dt < 0.0:
        return False
    cosang = (math.sin(p[2]) * math.sin(q[2]) * math.cos(p[3] - q[3])
              + math.cos(p[2]) * math.cos(q[2]))
    dpsi = math.acos(max(-1.0, min(1.0, cosang)))
    t_min, err = flight_time(p[1], q[1], dpsi, m, tol)
    if abs(dt - t_min) <= err:
        return None
    return bool(dt > t_min)


# ---------------------------------------------------------------------
# The benchmark: coordinate-uniform events on the patch
# ---------------------------------------------------------------------


def sample_events(n: int, rng: np.random.Generator) -> np.ndarray:
    """Coordinate-uniform events in the frozen patch (cost sampling
    only -- no volume-measure claim): t uniform, r^3 uniform in the
    shell, direction uniform on the polar cap."""

    t = rng.uniform(0.0, T_EXTENT, n)
    r = (rng.uniform(R_MIN ** 3, R_MAX ** 3, n)) ** (1.0 / 3.0)
    cos_lo = math.cos(CAP_HALF_ANGLE)
    theta = np.arccos(rng.uniform(cos_lo, 1.0, n))
    phi = rng.uniform(0.0, 2.0 * math.pi, n)
    return np.column_stack([t, r, theta, phi])


def bench(n: int = BENCH_N, seed: int = BENCH_SEED,
          tol: float = DEFAULT_TOL,
          k_pairs: int = 1_500) -> dict:
    """Predicate cost on `k_pairs` uniformly sampled ordered pairs
    from one `n`-event sample -- the PER-PAIR price is the quantity;
    an all-pairs census at probe sizes is precisely the cost this
    measurement exists to project, not to pay."""

    rng = np.random.default_rng(seed)
    ev = sample_events(n, rng)
    order = np.argsort(ev[:, 0], kind="stable")
    ev = ev[order]
    ii = rng.integers(0, n - 1, k_pairs)
    jj = ii + 1 + rng.integers(0, n - 1 - ii)
    related = undecided = 0
    start = time.perf_counter()
    for i, j in zip(ii, jj, strict=True):
        rel = causal_relation(ev[i], ev[j], M, tol)
        if rel is None:
            undecided += 1
        elif rel:
            related += 1
    seconds = time.perf_counter() - start
    return {"n": n, "pairs_sampled": k_pairs, "tol": tol,
            "seconds": seconds,
            "us_per_pair": seconds / k_pairs * 1e6,
            "related": related, "undecided": undecided}


def _plane_wave_reference() -> dict:
    """Same-host per-pair price of the PLANE-WAVE predicate at the
    probe chain's operating point -- the apples-to-apples baseline
    the Schwarzschild price is quoted against. Cost sampling only
    (seed recorded, never a statistics stream)."""

    from p14_plane_wave import Slab, arms, sprinkle
    from p14_probe_p1 import relation_census
    slab = Slab(du=1.0, dv=1.0, dx=2.0, dy=6.0)
    rho = 300 / slab.coordinate_volume
    curved, _ = arms(slab, 1.0)
    rng = np.random.default_rng(BENCH_SEED + 1)
    pts = sprinkle(curved, rho, rng)
    n = len(pts)
    start = time.perf_counter()
    relation_census(curved, pts)
    seconds = time.perf_counter() - start
    pairs = n * (n - 1) // 2
    return {"point": "aniso-a1.0", "n": n, "pairs": pairs,
            "seed": BENCH_SEED + 1,
            "us_per_pair": seconds / pairs * 1e6}


def main() -> None:
    print(f"S1 cost bench: shell [{R_MIN}, {R_MAX}] M, "
          f"cap {CAP_HALF_ANGLE} rad, n = {BENCH_N}")
    ladder = []
    for tol in TOL_LADDER:
        row = bench(tol=tol)
        ladder.append(row)
        print(f"  tol {tol:g}: {row['us_per_pair']:.1f} us/pair "
              f"({row['related']} related, "
              f"{row['undecided']} undecided)", flush=True)
    base = next(r for r in ladder if r["tol"] == DEFAULT_TOL)
    ref = _plane_wave_reference()
    print(f"  plane-wave reference: {ref['us_per_pair']:.2f} us/pair "
          f"-> price ratio {base['us_per_pair'] / ref['us_per_pair']:.0f}x",
          flush=True)
    projection = [
        {"n": n, "pairs": n * (n - 1) // 2,
         "seconds_per_sample":
             n * (n - 1) / 2 * base["us_per_pair"] / 1e6,
         "hours_per_4800_samples":
             4_800 * n * (n - 1) / 2 * base["us_per_pair"]
             / 1e6 / 3600.0}
        for n in PROJECTION_N]
    art = {
        "domain": {"m": M, "r_shell": [R_MIN, R_MAX],
                   "cap_half_angle": CAP_HALF_ANGLE,
                   "t_extent": T_EXTENT, "psi_max": PSI_MAX,
                   "bench_seed": BENCH_SEED, "bench_n": BENCH_N},
        "scope": "이 가격은 해당 solver·patch·N 경로만 연다/보류한다; "
                 "감당 가능해도 Schwarzschild 부피·캠페인 경로는 별도 "
                 "미해결이다 (p14_weyl_curvature.md §8.1)",
        "host": {"machine": platform.machine(),
                 "python": platform.python_version()},
        "ladder": ladder,
        "default_tol": DEFAULT_TOL,
        "projection_at_default_tol": projection,
        "plane_wave_reference": ref,
        "price_ratio_at_default_tol":
            base["us_per_pair"] / ref["us_per_pair"],
    }
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(art, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"artifact: {_ARTIFACT}")


if __name__ == "__main__":
    sys.exit(main())
