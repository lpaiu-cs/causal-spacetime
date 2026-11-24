"""Referee check (f), (g): assembly facts, the algebra (7) -> (8), the regime.

Part 1 (probability facts used in Section 5, with the convention
H^2(P,Q) = sum (sqrt p - sqrt q)^2, i.e. NOT divided by 2):
  H1  TV <= H                                   (random discrete laws)
  H2  H^2(product) <= sum of H^2                (random product laws)
  H3  H^2(mixture) <= mixture of H^2            (averaging over spatial coords)
  H4  TV(joint of (W,O)) = E_W TV(O|W) when W has the same law in both models
  H5  E sqrt(U) <= sqrt(E U)                    (the Cauchy-Schwarz step)

Part 2 (SymPy, exact exponents).  With n = log N, d = log(1/delta) > 0 and
log h = (2/(2K+1)) (-d - K n):
  X1  N^{K+1} h^K = N delta / sqrt(h) = (N^{3K+1} delta^{2K})^{1/(2K+1)}
  X2  N h = (N delta^2)^{1/(2K+1)}
  X3  N^2 h^4 = (N h)(N h^3) and N h^3 <= (N h)^3
  X4  h is the unique balance point of the two h-dependent main terms
Part 3 (regime, N = lambda^-alpha, delta = lambda^2, lambda -> 0): every term
of (7) has a positive lambda-exponent and h/delta -> infinity, h -> 0,
N h -> 0 exactly when alpha < 4K/(3K+1); threshold sequence 8/7, 6/5, 16/13;
it increases to 4/3 and never reaches it; K=1 returns alpha < 1 (the old
N lambda -> 0 bound); at alpha = 4/3 the main exponent is -4/(3(2K+1)) < 0 for
every K (no information).  A numeric minimisation over h confirms that no
other h does better for the bound (7) as written.

Deterministic.  Run: python check_f_assembly_exponents.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from audit_common import Reporter

rep = Reporter("check (f,g): assembly inequalities, exponent algebra, regime")
rng = np.random.default_rng(99)


def hell2(p, q):
    return float(((np.sqrt(p) - np.sqrt(q)) ** 2).sum())


def tv(p, q):
    return 0.5 * float(np.abs(p - q).sum())


def rand_law(k, sparse=False):
    x = rng.gamma(0.5, size=k)
    if sparse:
        x[rng.uniform(size=k) < 0.4] = 0.0
        if x.sum() == 0:
            x[0] = 1.0
    return x / x.sum()


# ---------------------------------------------------------------- H1..H5
ok1 = ok2 = ok3 = ok4 = True
worst1 = worst2 = 0.0
for _ in range(4000):
    k = int(rng.integers(2, 9))
    p, q = rand_law(k, True), rand_law(k, True)
    h = np.sqrt(hell2(p, q))
    t = tv(p, q)
    worst1 = max(worst1, t / h if h > 0 else 0.0)
    ok1 &= t <= h * (1 + 1e-12)
    # products of 3 independent factors
    ps = [rand_law(int(rng.integers(2, 5)), True) for _ in range(3)]
    qs = [rand_law(len(x), True) for x in ps]
    P = ps[0][:, None, None] * ps[1][None, :, None] * ps[2][None, None, :]
    Q = qs[0][:, None, None] * qs[1][None, :, None] * qs[2][None, None, :]
    lhs = hell2(P.ravel(), Q.ravel())
    rhs = sum(hell2(a, b) for a, b in zip(ps, qs, strict=False))
    worst2 = max(worst2, lhs / rhs if rhs > 0 else 0.0)
    ok2 &= lhs <= rhs * (1 + 1e-12) + 1e-15
    # mixtures with common weights
    wts = rand_law(4)
    comp_p = [rand_law(k, True) for _ in range(4)]
    comp_q = [rand_law(k, True) for _ in range(4)]
    mp_ = sum(w * c for w, c in zip(wts, comp_p, strict=False))
    mq_ = sum(w * c for w, c in zip(wts, comp_q, strict=False))
    ok3 &= (
        hell2(mp_, mq_)
        <= sum(w * hell2(a, b) for w, a, b in zip(wts, comp_p, comp_q, strict=False))
        * (1 + 1e-12)
        + 1e-15
    )
    # joint with common first marginal
    jp = np.stack([w * c for w, c in zip(wts, comp_p, strict=False)])
    jq = np.stack([w * c for w, c in zip(wts, comp_q, strict=False)])
    ok4 &= (
        abs(
            tv(jp.ravel(), jq.ravel())
            - sum(w * tv(a, b) for w, a, b in zip(wts, comp_p, comp_q, strict=False))
        )
        < 1e-12
    )
rep.check(
    "H1 TV <= H with H^2 = sum (sqrt p - sqrt q)^2", ok1, f"max TV/H = {worst1:.4f}"
)
rep.check(
    "H2 H^2 of independent products <= sum of H^2", ok2, f"max ratio {worst2:.4f}"
)
rep.check(
    "H3 H^2 jointly convex (mixing over spatial coordinates cannot increase it)", ok3
)
rep.check(
    "H4 TV of (marks, order) = E_marks TV(order | marks) for a common mark law", ok4
)
u_counts = rng.poisson(3.0, size=100000)
rep.check(
    "H5 E sqrt(U) <= sqrt(E U)",
    float(np.sqrt(u_counts).mean()) <= float(np.sqrt(u_counts.mean())),
)

# ---------------------------------------------------------------- X1..X4
K = sp.symbols("K", positive=True)
n, d = sp.symbols("n d", positive=True)  # n = log N, d = log(1/delta)
logh = sp.Rational(2) / (2 * K + 1) * (-d - K * n)
tree = (K + 1) * n + K * logh
hellt = n - d - logh / 2
target = ((3 * K + 1) * n - 2 * K * d) / (2 * K + 1)
rep.check(
    "X1 N^{K+1} h^K = N delta h^{-1/2} = (N^{3K+1} delta^{2K})^{1/(2K+1)}",
    sp.simplify(tree - target) == 0 and sp.simplify(hellt - target) == 0,
)
rep.check(
    "X2 N h = (N delta^2)^{1/(2K+1)}",
    sp.simplify(n + logh - (n - 2 * d) / (2 * K + 1)) == 0,
)
rep.check(
    "X3 log(N^2 h^4) - log(N h) = log(N h^3) = 3 log(N h) - 2 n",
    sp.simplify((2 * n + 4 * logh) - (n + logh) - (3 * (n + logh) - 2 * n)) == 0,
)
x = sp.symbols("x", real=True)  # x = log h free
bal = sp.solve(sp.Eq((K + 1) * n + K * x, n - d - x / 2), x)
rep.check(
    "X4 the stated h is the balance point of N^{K+1}h^K and N delta/sqrt(h)",
    len(bal) == 1 and sp.simplify(bal[0] - logh) == 0,
)

# ---------------------------------------------------------------- regime
alpha = sp.symbols("alpha", positive=True)
# N = lambda^-alpha, delta = lambda^2; write everything as lambda^(exponent)
subs = {
    n: -alpha,
    d: -2,
}  # units of log(lambda): log N = -alpha log(lam), log(1/delta) = -2 log(lam)
ex_main = sp.simplify(target.subs(subs))  # exponent of lambda in the main term
ex_Nh = sp.simplify((n + logh).subs(subs))
ex_N2h4 = sp.simplify((2 * n + 4 * logh).subs(subs))
ex_Ndelta = sp.simplify((n - d).subs(subs))
ex_h = sp.simplify(logh.subs(subs))
ex_h_over_delta = sp.simplify((logh + d).subs(subs))
thr = sp.solve(sp.Eq(ex_main, 0), alpha)
rep.check(
    "R1 main exponent (4K - (3K+1) alpha)/(2K+1): positive iff alpha < 4K/(3K+1)",
    sp.simplify(ex_main - (4 * K - (3 * K + 1) * alpha) / (2 * K + 1)) == 0
    and sp.simplify(thr[0] - 4 * K / (3 * K + 1)) == 0,
    f"exponent {sp.sstr(sp.factor(ex_main))}",
)
ok = True
rows = []
for Kv in range(1, 41):
    a_thr = sp.Rational(4 * Kv, 3 * Kv + 1)
    for av in (sp.Rational(1, 100), a_thr / 2, a_thr - sp.Rational(1, 10**6)):
        vals = [
            e.subs({K: Kv, alpha: av})
            for e in (ex_main, ex_Nh, ex_N2h4, ex_Ndelta, ex_h)
        ]
        if not all(v > 0 for v in vals):
            ok = False
        if (
            not ex_h_over_delta.subs({K: Kv, alpha: av}) < 0
        ):  # h/delta = lambda^(negative) -> infinity
            ok = False
    # just above the threshold the main term no longer vanishes
    if not ex_main.subs({K: Kv, alpha: a_thr + sp.Rational(1, 10**6)}) < 0:
        ok = False
    if not a_thr < sp.Rational(4, 3):
        ok = False
    if Kv <= 4:
        rows.append(f"K={Kv}:{a_thr}")
rep.check(
    "R2 for K=1..40 and alpha below 4K/(3K+1): all five terms of (7) vanish, h->0, "
    "h/delta->infinity",
    ok,
    " ".join(rows),
)
rep.check(
    "R3 thresholds 4K/(3K+1) increase strictly to 4/3 (limit not attained)",
    sp.limit(4 * K / (3 * K + 1), K, sp.oo) == sp.Rational(4, 3)
    and sp.simplify(sp.diff(4 * K / (3 * K + 1), K) - 4 / (3 * K + 1) ** 2) == 0,
)
rep.check(
    "R4 K=1 gives alpha < 1 only (N^4 lambda^4 -> 0, i.e. the old N lambda -> 0 bound)",
    sp.Rational(4 * 1, 3 * 1 + 1) == 1,
)
at43 = sp.simplify(ex_main.subs(alpha, sp.Rational(4, 3)))
rep.check(
    "R5 at alpha = 4/3 the main exponent is -4/(3(2K+1)) < 0 for every K: (8) says "
    "nothing",
    sp.simplify(at43 + sp.Rational(4, 3) / (2 * K + 1)) == 0,
    f"exponent {sp.sstr(at43)}",
)
# h/delta -> infinity needs alpha < 2 - 1/K; this is implied by alpha < 4K/(3K+1)
gap = sp.simplify((2 - 1 / K) - 4 * K / (3 * K + 1))
rep.check(
    "R6 the side condition h >> delta (alpha < 2 - 1/K) is implied: (2 - 1/K) - "
    "4K/(3K+1) = (2K+1)(K-1)/(K(3K+1)) >= 0",
    sp.simplify(gap - (2 * K + 1) * (K - 1) / (K * (3 * K + 1))) == 0,
)
# for a given alpha < 4/3 the smallest admissible K
ok = True
detail = []
for av in (
    sp.Rational(1),
    sp.Rational(11, 10),
    sp.Rational(6, 5),
    sp.Rational(5, 4),
    sp.Rational(13, 10),
    sp.Rational(133, 100),
):
    Kmin = sp.floor(av / (4 - 3 * av)) + 1
    if not (
        sp.Rational(4 * Kmin, 3 * Kmin + 1) > av
        and (Kmin == 1 or sp.Rational(4 * (Kmin - 1), 3 * (Kmin - 1) + 1) <= av)
    ):
        ok = False
    detail.append(f"alpha={av}->K={Kmin}")
rep.check(
    "R7 every fixed alpha < 4/3 is covered by K = floor(alpha/(4-3alpha)) + 1",
    ok,
    " ".join(detail),
)


# numeric: (7) at the stated h versus (8), and versus the best h for (7)
def bound7(Nv, dv, hs, Kv):
    return (
        Nv * dv
        + Nv * hs
        + Nv**2 * hs**4
        + Nv ** (Kv + 1) * hs**Kv
        + Nv * dv / np.sqrt(hs)
    )


ok = True
worst = 0.0
for Kv in (2, 3, 4, 6):
    for lam_e in (3, 5, 8):
        lamv = 10.0 ** (-lam_e)
        for frac in (0.5, 0.9):
            av = frac * 4 * Kv / (3 * Kv + 1)
            Nv = lamv ** (-av)
            dv = lamv**2
            hs = np.logspace(np.log10(dv) + 0.01, 0, 20001)
            fmin = bound7(Nv, dv, hs, Kv).min()
            h_star = (dv / Nv**Kv) ** (2.0 / (2 * Kv + 1))
            f_star = bound7(Nv, dv, h_star, Kv)
            rhs8 = (
                Nv * lamv**2
                + (Nv * lamv**4) ** (1.0 / (2 * Kv + 1))
                + (Nv ** (3 * Kv + 1) * lamv ** (4 * Kv)) ** (1.0 / (2 * Kv + 1))
            )
            worst = max(worst, f_star / fmin)
            # (8) with constant 3 bounds (7) at h_star (N^2 h^4 <= N h there)
            if not (f_star <= 3.0 * rhs8 * (1 + 1e-9) and h_star > dv and h_star < 1.0):
                ok = False
rep.check(
    "R8 numeric: (7) evaluated at the stated h is <= 3 x RHS(8), with delta < h < 1", ok
)
rep.info(
    "R8 note: the stated h is NOT the minimiser of (7) when N h dominates (small "
    "alpha); "
    f"max f(h*)/min_h f = {worst:.1f}.  This costs rate, not validity: (8) is still "
    "an upper bound."
)
# the vanishing condition of (7) optimised over h is exactly (1): no better
# exponent can be extracted from (7) at fixed K by another choice of h.
ok = True
detail = []
for Kv in (2, 3, 4, 8):
    thr_v = 4 * Kv / (3 * Kv + 1)
    below = []
    above = []
    for lam_e in (4, 8, 12):
        lamv = 10.0 ** (-lam_e)
        dv = lamv**2
        hs = np.logspace(np.log10(dv) + 0.01, 0, 40001)
        below.append(bound7(lamv ** (-0.97 * thr_v), dv, hs, Kv).min())
        above.append(bound7(lamv ** (-1.03 * thr_v), dv, hs, Kv).min())
    if not (below[0] > below[1] > below[2] and above[0] < above[1] < above[2]):
        ok = False
    detail.append(
        f"K={Kv} below:{below[0]:.2e}->{below[2]:.2e} "
        f"above:{above[0]:.2e}->{above[2]:.2e}"
    )
rep.check(
    "R9 min over h of (7) tends to 0 just below alpha = 4K/(3K+1) and grows just "
    "above it",
    ok,
    "; ".join(detail),
)

rep.finish()
