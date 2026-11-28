"""Referee check (e), part 1: the finite-component smoothing lemma, exactly.

p(tau) = Vol{ t in [0,1]^k : sigma_e (t_j - t_i) >= tau_e for every row e }
for independent UNIFORM times (the h-scaling is p_h(tau) = p_1(tau/h), so a
Lipschitz constant L on the unit cube is L/h^2 on a cube of side h).

Everything is exact rational arithmetic: volumes by the Lasserre recursion and
the gradient as exact facet volumes,  d p / d tau_e = -(projected facet volume).
No smoothing, no finite-difference gradient, no square-root bound is assumed.

Referee's explicit constants for uniform density on the unit cube:
    |d^2 p / d tau_e d tau_e'| <= 1  (two different rows, sharing a vertex or not)
    |d^2 p / d tau_e^2|        <= 2  (slice of a polytope none of whose facets
                                      is parallel to the slicing hyperplane)
    => Hessian operator norm <= m + 1 =: L   (Gershgorin)

Tests, for k = 2, 3, 4 and many sign patterns (including the cyclic triangle,
rows sharing a vertex, and the complete graph on 4 vertices whose 6 rows have
rank 3):
  E1  entrywise difference quotients of the exact gradient: mixed <= 1,
      repeated <= 2, at step sizes 1/8 ... 1/4096 (no blow-up at small steps).
  E2  |grad p(tau') - grad p(tau)| <= (m+1) |tau' - tau|.
  E3  |grad p|^2 <= 2 (m+1) p, including points at distance 2^-13 from the
      zero set of p; the largest observed |grad p|^2 / (2p) is reported.
  E4  (sqrt p(tau') - sqrt p(tau))^2 <= ((m+1)/2) |tau'-tau|^2, including pairs
      straddling the zero set (patterns of zero reference probability).
  E5  cross-checks of the volume code: the eight triangle patterns sum to 1,
      the cyclic pattern equals (9/2) eps^2 - 4 eps^3, and the author's
      recorded Hellinger number 13.44677 at eps = 0.001 is reproduced.
  E6  sandwich -> Hellinger: for thresholds moved by +-eps in the tightening /
      relaxing direction, sum_S (sqrt p_S(relaxed) - sqrt p_S(tightened))^2
      <= 2^m ((m+1)/2) m (2 eps)^2, and sum_S p_S(tightened) <= 1 <= sum_S p_S(relaxed).
  E7  SymPy, k = 2: p'' is piecewise +-1 (so L = 1 = 1/h^2 after scaling) and
      sqrt(p) has slope exactly 1/sqrt(2) at the zero set: the constant
      sqrt(L/2) cannot be improved.

Deterministic (fixed seed).  Run: python check_e_smoothing_lemma.py
"""

from __future__ import annotations

import itertools
import random
from fractions import Fraction as Fr

import sympy as sp
from audit_common import Reporter, pattern_grad, pattern_prob

rep = Reporter("check (e1): smoothing lemma, exact rational polytope volumes")
rnd = random.Random(20261124)

CONFIGS = [
    ("k2 single row", 2, [(0, 1)]),
    ("k3 chain (rows share vertex 1)", 3, [(0, 1), (1, 2)]),
    ("k3 star (rows share vertex 0)", 3, [(0, 1), (0, 2)]),
    ("k3 vee (rows share vertex 2)", 3, [(0, 2), (1, 2)]),
    ("k3 triangle (3 dependent rows)", 3, [(0, 1), (1, 2), (0, 2)]),
    ("k4 path", 4, [(0, 1), (1, 2), (2, 3)]),
    ("k4 cycle (4 dependent rows)", 4, [(0, 1), (1, 2), (2, 3), (0, 3)]),
    (
        "k4 complete (6 rows, rank 3)",
        4,
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    ),
    ("k5 cycle (5 dependent rows)", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]),
    (
        "k5 complete (10 rows, rank 4)",
        5,
        [(i, j) for i in range(5) for j in range(i + 1, 5)],
    ),
]


def sign_patterns(m):
    allp = list(itertools.product((1, -1), repeat=m))
    if m <= 3:
        return allp
    if m == 4:
        return [allp[0], allp[5], allp[10], allp[15], allp[1], allp[6]]
    if m != 6:
        # k = 5: the all-true pattern plus two seeded random ones
        return [tuple([1] * m)] + [
            tuple(rnd.choice((1, -1)) for _ in range(m)) for _ in range(2)
        ]
    # complete graph on 4 vertices: a transitive pattern, cyclic ones, mixed ones
    return [
        (1, 1, 1, 1, 1, 1),
        (1, -1, 1, 1, 1, 1),  # 0<1, 1<2, but not 0<2: cyclic triangle inside
        (1, 1, -1, 1, 1, 1),
        (-1, 1, -1, 1, -1, 1),
    ]


def rand_tau(m, den=64, span=80):
    return [Fr(rnd.randint(-span, span), den) for _ in range(m)]


def norm2(vec):
    return sum(x * x for x in vec)


def sqrt_lipschitz_ok(p1, p2, c_d2):
    """Exact test of (sqrt p1 - sqrt p2)^2 <= c_d2 for rationals."""
    lhs = p1 + p2 - c_d2
    if lhs <= 0:
        return True
    return lhs * lhs <= 4 * p1 * p2


STEPS = [Fr(1, 8), Fr(1, 64), Fr(1, 512), Fr(1, 4096)]

for name, k, pairs in CONFIGS:
    m = len(pairs)
    L = m + 1
    npts = {2: 10, 3: 6, 4: 2, 5: 1}[k]
    near_exp = {2: (4, 7, 10, 13), 3: (4, 7, 10, 13), 4: (5, 11), 5: (9,)}[k]
    corner_exp = {2: (4, 8, 12), 3: (4, 8, 12), 4: (4, 10), 5: (8,)}[k]
    steps = STEPS if k <= 4 else [Fr(1, 64), Fr(1, 4096)]
    max_mixed = Fr(0)
    max_rep = Fr(0)
    max_gl = Fr(0)
    max_ratio = Fr(0)
    e1 = e2 = e3 = e4 = True
    n_zero_cross = 0
    n_eval = 0
    pats = sign_patterns(m)
    for sg in pats:
        rows = [(i, j, s) for (i, j), s in zip(pairs, sg, strict=False)]
        points = []
        # generic points
        for _ in range(npts):
            points.append(rand_tau(m))
        # the degenerate corner used in the author's control: every signed
        # threshold equal to -eps (for the cyclic triangle p = 9/2 eps^2 - 4 eps^3)
        for r in corner_exp:
            points.append([Fr(-1, 2**r)] * m)
        # points close to the zero set: walk along the diagonal from a point
        # with p > 0 until p = 0 (dyadic bisection), then step back by 2^-r
        tries = 0
        base = None
        while base is None and tries < 40:
            cand = rand_tau(m, span=48)
            tries += 1
            if pattern_prob(k, rows, cand) > 0:
                base = cand
        if base is not None:
            lo, hi = Fr(0), Fr(4)
            for _ in range(15):
                mid = (lo + hi) / 2
                if pattern_prob(k, rows, [t + mid for t in base]) > 0:
                    lo = mid
                else:
                    hi = mid
            for r in near_exp:
                s_in = hi - Fr(1, 2**r)
                if s_in > 0:
                    points.append([t + s_in for t in base])
            # a pair straddling the zero set
            inside = [t + hi - Fr(1, 2**9) for t in base]
            outside = [t + hi + Fr(1, 2**9) for t in base]
            p_in = pattern_prob(k, rows, inside)
            p_out = pattern_prob(k, rows, outside)
            d2 = norm2([a - b for a, b in zip(inside, outside, strict=False)])
            if p_out == 0 and p_in > 0:
                n_zero_cross += 1
            e4 &= sqrt_lipschitz_ok(p_in, p_out, Fr(L, 2) * d2)
        for tau in points:
            p0 = pattern_prob(k, rows, tau)
            g0 = pattern_grad(k, rows, tau)
            n_eval += 1
            g2 = norm2(g0)
            if g2 > 2 * L * p0:
                e3 = False
            if p0 > 0:
                max_ratio = max(max_ratio, g2 / (2 * p0))
            elif g2 != 0:
                e3 = False
            for eta in steps:
                # coordinate direction -> entrywise Hessian difference quotients
                ep = rnd.randrange(m)
                tau1 = list(tau)
                tau1[ep] += eta
                g1 = pattern_grad(k, rows, tau1)
                for e in range(m):
                    quo = abs(g1[e] - g0[e]) / eta
                    if e == ep:
                        max_rep = max(max_rep, quo)
                        if quo > 2:
                            e1 = False
                    else:
                        max_mixed = max(max_mixed, quo)
                        if quo > 1:
                            e1 = False
                # random direction -> operator norm and sqrt Lipschitz
                dvec = [Fr(rnd.choice((-1, 0, 1, 1))) * eta for _ in range(m)]
                if norm2(dvec) == 0:
                    dvec[0] = eta
                tau2 = [a + b for a, b in zip(tau, dvec, strict=False)]
                g2v = pattern_grad(k, rows, tau2)
                d2 = norm2(dvec)
                dg2 = norm2([a - b for a, b in zip(g2v, g0, strict=False)])
                if dg2 > L * L * d2:
                    e2 = False
                max_gl = max(max_gl, dg2 / d2)
                p2 = pattern_prob(k, rows, tau2)
                e4 &= sqrt_lipschitz_ok(p0, p2, Fr(L, 2) * d2)
    tag = f"{name}, m={m}, L=m+1={L}, {len(pats)} sign patterns, {n_eval} base points"
    rep.check(
        "E1 Hessian entries bounded (mixed<=1, repeated<=2) " + tag,
        e1,
        f"max mixed quotient {float(max_mixed):.4f}, max repeated quotient "
        f"{float(max_rep):.4f}",
    )
    rep.check(
        "E2 gradient Lipschitz with constant m+1 " + name,
        e2,
        "max |dgrad|/|dtau| = %.4f" % float(max_gl) ** 0.5,
    )
    rep.check(
        "E3 |grad p|^2 <= 2 L p " + name,
        e3,
        f"max |grad p|^2/(2p) = {float(max_ratio):.4f} (lower bound on the true L)",
    )
    rep.check(
        "E4 sqrt(p) Lipschitz with constant sqrt(L/2), incl. zero-probability side "
        + name,
        e4,
        f"zero-set crossings tested: {n_zero_cross}",
    )

# ---------------------------------------------------------------- E5
tri = [(0, 1), (1, 2), (2, 0)]  # the author's cyclic orientation: t2-t1, t3-t2, t1-t3


def tri_patterns(eps):
    """Author's control: True means edge.t >= -eps, False means edge.t <= -eps."""
    out = []
    for bits in itertools.product((False, True), repeat=3):
        rows = []
        tau = []
        for (i, j), bt in zip(tri, bits, strict=False):
            if bt:
                rows.append((i, j, 1))
                tau.append(-eps)
            else:
                rows.append((i, j, -1))
                tau.append(eps)
        out.append(pattern_prob(3, rows, tau))
    return out


ok5 = True
detail = []
base = tri_patterns(Fr(0))
for eps in (Fr(0), Fr(1, 1000), Fr(1, 100), Fr(1, 10), Fr(2, 5)):
    pr = tri_patterns(eps)
    cyc = pr[-1]
    expected = Fr(9, 2) * eps**2 - 4 * eps**3
    if sum(pr) != 1 or cyc != expected:
        ok5 = False
    detail.append(f"eps={eps} cyc={cyc}")
rep.check(
    "E5a triangle: patterns sum to 1 and cyclic mass = (9/2)eps^2 - 4 eps^3 exactly",
    ok5,
    "; ".join(detail),
)
pr = tri_patterns(Fr(1, 1000))
hell = sum(
    (float(a) ** 0.5 - float(b) ** 0.5) ** 2 for a, b in zip(pr, base, strict=False)
)
rep.check(
    "E5b author's recorded H^2/eps^2 = 13.44677 at eps = 0.001 reproduced",
    abs(hell * 1e6 - 13.44677) < 1e-4,
    "referee value %.5f" % (hell * 1e6),
)

# ---------------------------------------------------------------- E6
for name, k, pairs in CONFIGS[:5]:
    m = len(pairs)
    L = m + 1
    ok6 = True
    worst = 0.0
    for trial in range(6):
        theta = rand_tau(m, span=40) if trial else [Fr(0)] * m
        for eps in (Fr(1, 16), Fr(1, 256), Fr(1, 4096)):
            tight_sum = Fr(0)
            relax_sum = Fr(0)
            base_sum = Fr(0)
            h2 = 0.0
            for sg in itertools.product((1, -1), repeat=m):
                rows = [(i, j, s) for (i, j), s in zip(pairs, sg, strict=False)]
                tau = [s * t for s, t in zip(sg, theta, strict=False)]
                pt = pattern_prob(k, rows, [t + eps for t in tau])
                pr_ = pattern_prob(k, rows, [t - eps for t in tau])
                pb = pattern_prob(k, rows, tau)
                if not (pt <= pb <= pr_):
                    ok6 = False
                tight_sum += pt
                relax_sum += pr_
                base_sum += pb
                h2 += (float(pr_) ** 0.5 - float(pt) ** 0.5) ** 2
            bound = (2**m) * (L / 2.0) * m * float(2 * eps) ** 2
            worst = max(worst, h2 / bound)
            if (
                base_sum != 1
                or tight_sum > 1
                or relax_sum < 1
                or h2 > bound * (1 + 1e-12)
            ):
                ok6 = False
    rep.check(
        "E6 sandwich -> Hellinger bound " + name, ok6, f"max H^2/bound = {worst:.4f}"
    )

# ---------------------------------------------------------------- E7
tau = sp.symbols("tau", real=True)
t1, t2 = sp.symbols("t1 t2", real=True)
# p(tau) = P(t2 - t1 >= tau), uniform on the unit square
p_mid_neg = 1 - (1 + tau) ** 2 / 2  # -1 <= tau <= 0
p_mid_pos = (1 - tau) ** 2 / 2  # 0 <= tau <= 1
chk = sp.integrate(sp.integrate(1, (t2, t1 + tau, 1)), (t1, 0, 1 - tau))
ok7 = sp.simplify(chk - p_mid_pos) == 0
ok7 &= sp.diff(p_mid_neg, tau, 2) == -1 and sp.diff(p_mid_pos, tau, 2) == 1
ok7 &= sp.simplify(p_mid_neg.subs(tau, 0) - p_mid_pos.subs(tau, 0)) == 0
ok7 &= (
    sp.simplify(
        sp.diff(p_mid_neg, tau).subs(tau, 0) - sp.diff(p_mid_pos, tau).subs(tau, 0)
    )
    == 0
)
slope = sp.simplify(sp.diff(sp.sqrt(p_mid_pos), tau).subs(tau, sp.Rational(9, 10)))
ok7 &= sp.simplify(slope + 1 / sp.sqrt(2)) == 0
exact_k2 = pattern_prob(2, [(0, 1, 1)], [Fr(3, 10)]) == Fr(49, 200)
rep.check(
    "E7 k=2 closed form: C^{1,1}, p'' = -1 / +1, sqrt(p) slope = 1/sqrt(2) = "
    "sqrt(L/2) with L=1",
    bool(ok7) and exact_k2,
    f"sqrt-slope {sp.sstr(slope)}",
)

rep.finish()
