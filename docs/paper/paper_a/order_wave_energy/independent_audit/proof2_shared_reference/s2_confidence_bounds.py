# ruff: noqa: E501, UP031
"""Referee check of Section 2 of local_null_challenge.md (statistics only).

A. Bernstein radius algebra behind eq. (3).
B. Exact binomial tails: the Bernstein radius with the TRUE variance has two-sided
   failure <= 2 exp(-B).
C. Exact binomial computation of the data-driven radius r_D of eq. (3)
   (Bernstein with the empirical upper bound d_U, on the Hoeffding event):
   failure <= 4 exp(-B) = Hoeffding budget + Bernstein budget.
D. Order-two U-statistic with nested {0,1} kernels on a toy partial order:
   nonnegativity of the three increments, Var(U) <= sigma^2/k with k=floor(M/2),
   coverage of r_D, and a negative control with the naive k=choose(M,2).
E. Union-bound arithmetic: 2^B >= 32 N^4 gives simultaneous coverage >= 1-N^-2.
F. Deterministic propagation bound (4) on the domain R = B_b + hull(0,B_d).
G. Exponent bookkeeping for (6)-(8).

Deterministic (fixed seeds).  Prints PASS/FAIL lines; exits nonzero on any FAIL.
Run:  python s2_confidence_bounds.py      (system Python: numpy + sympy)
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np
import sympy as sp

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


# --------------------------------------------------------------------------- A
def part_a():
    k, b, s2 = sp.symbols("k B sigma2", positive=True)
    tstar = b / (3 * k) + sp.sqrt(b**2 / (9 * k**2) + 2 * s2 * b / k)
    # t* solves k t^2 = B (2 sigma^2 + 2 t / 3): the Bernstein exponent equals B.
    residual = sp.simplify(sp.expand(k * tstar**2 - b * (2 * s2 + 2 * tstar / 3)))
    check("A1 Bernstein exponent equals B at t*", residual == 0, f"residual={residual}")
    stated = sp.sqrt(2 * s2 * b / k) + 2 * b / (3 * k)
    # (stated - B/3k)^2 - (t* - B/3k)^2 must be a manifestly nonnegative expression.
    gap = sp.simplify(
        sp.expand((stated - b / (3 * k)) ** 2 - (tstar - b / (3 * k)) ** 2)
    )
    target = 2 * b * sp.sqrt(2 * s2 * b / k) / (3 * k)
    check(
        "A2 stated radius sqrt(2 s2 B/k)+2B/(3k) dominates the exact root",
        sp.simplify(gap - target) == 0,
        f"gap={gap}",
    )
    rh = sp.sqrt(b / (2 * k))
    check(
        "A3 Hoeffding radius r_H=sqrt(B/2k): 2 exp(-2 k r^2) = 2 exp(-B)",
        sp.simplify(2 * k * rh**2 - b) == 0,
    )


# --------------------------------------------------------------------------- B, C
def log_pmf(k, mu):
    j = np.arange(k + 1)
    lg = np.array(
        [math.lgamma(k + 1) - math.lgamma(i + 1) - math.lgamma(k - i + 1) for i in j]
    )
    with np.errstate(divide="ignore"):
        return lg + j * math.log(mu) + (k - j) * math.log1p(-mu)


def part_bc():
    worst_true, worst_data, worst_data_vs_single = 0.0, 0.0, 0.0
    cases = 0
    for k in (5, 20, 100, 500, 3000):
        for mu in (1e-4, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.8, 0.97):
            pmf = np.exp(log_pmf(k, mu))
            mean = np.arange(k + 1) / k
            for b in (0.5, 1, 2, 3, 5, 8, 12):
                cases += 1
                sigma2 = mu * (1 - mu)
                r_true = math.sqrt(2 * sigma2 * b / k) + 2 * b / (3 * k)
                fail_true = pmf[np.abs(mean - mu) > r_true].sum()
                worst_true = max(worst_true, fail_true / (2 * math.exp(-b)))
                r_h = math.sqrt(b / (2 * k))
                d_u = np.minimum(1.0, mean + r_h)
                r_d = np.minimum(r_h, np.sqrt(2 * d_u * b / k) + 2 * b / (3 * k))
                fail_data = pmf[np.abs(mean - mu) > r_d].sum()
                worst_data = max(worst_data, fail_data / (4 * math.exp(-b)))
                worst_data_vs_single = max(
                    worst_data_vs_single, fail_data / (2 * math.exp(-b))
                )
    check(
        "B  exact binomial: Bernstein(true variance) failure <= 2 exp(-B)",
        worst_true <= 1.0,
        f"{cases} cases, worst ratio to budget {worst_true:.3e}",
    )
    check(
        "C  exact binomial: data-driven r_D of eq.(3) failure <= 4 exp(-B)",
        worst_data <= 1.0,
        f"worst ratio to 4exp(-B) {worst_data:.3e}; to 2exp(-B) {worst_data_vs_single:.3e}",
    )


# --------------------------------------------------------------------------- D
def part_d():
    rng = np.random.default_rng(20261124)
    reps, m = 20000, 40
    k2 = m // 2
    x = np.array([0.5, 0.5])
    y = np.array([0.56, 0.58])  # X < Y in the product order on the unit square
    z = rng.uniform(0.0, 1.0, (reps, m, 2))
    less = (z[:, :, None, 0] < z[:, None, :, 0]) & (z[:, :, None, 1] < z[:, None, :, 1])
    comparable = less | np.swapaxes(less, 1, 2)
    below_x = np.all(z < x, axis=2)
    below_y = np.all(z < y, axis=2)
    above_x = np.all(z > x, axis=2)
    above_y = np.all(z > y, axis=2)
    pairs = math.comb(m, 2)

    def a_hat(mask):
        both = mask[:, :, None] & mask[:, None, :] & comparable
        return both.sum(axis=(1, 2)) / 2 / pairs

    p_x, p_y = below_x.mean(axis=1), below_y.mean(axis=1)
    f_x, f_y = above_x.mean(axis=1), above_y.mean(axis=1)
    a_x, a_y = a_hat(below_x), a_hat(below_y)
    inc = np.stack([p_y - p_x, f_x - f_y, a_y - a_x], axis=1)
    check(
        "D1 nestedness: the three empirical increments are >= 0 in every replicate",
        bool(np.all(inc >= 0)),
        f"min={inc.min():.3g} over {reps} replicates",
    )
    # exact population values on the unit square (two uniform points of a rectangle
    # are comparable with probability 1/2)
    px, py = x.prod(), y.prod()
    fx, fy = (1 - x).prod(), (1 - y).prod()
    truth = np.array([py - px, fx - fy, (py**2 - px**2) / 2])
    bias = np.abs(inc.mean(axis=0) - truth) / (inc.std(axis=0) / math.sqrt(reps))
    check(
        "D2 a = p^2 c and the increment means (toy order, c=1/2) are unbiased",
        bool(np.all(bias < 5.0)),
        f"z-scores {np.array2string(bias, precision=2)}",
    )
    sizes = np.array([m, m, k2])
    # kernel variance of a {0,1} kernel = mean (1-mean) <= mean
    var_ratio = inc.var(axis=0) / (truth * (1 - truth) / sizes)
    check(
        "D3 Var(estimate) <= kernel variance / k, k=M (p,f) and floor(M/2) (a)",
        bool(np.all(var_ratio <= 1.05)),
        f"variance ratios {np.array2string(var_ratio, precision=3)}",
    )
    ok_all, report = True, []
    for b in (1.0, 2.0, 3.0, 4.0):
        fails = []
        for i in range(3):
            kk = sizes[i]
            r_h = math.sqrt(b / (2 * kk))
            d_u = np.minimum(1.0, inc[:, i] + r_h)
            r_d = np.minimum(r_h, np.sqrt(2 * d_u * b / kk) + 2 * b / (3 * kk))
            fails.append(float(np.mean(np.abs(inc[:, i] - truth[i]) > r_d)))
        budget = 4 * math.exp(-b)
        ok_all &= max(fails) <= budget
        report.append(f"B={b:g}: max fail {max(fails):.4f} <= {budget:.4f}")
    check(
        "D4 U-statistic increments: r_D coverage with k=floor(M/2)",
        ok_all,
        "; ".join(report),
    )
    # Negative control: treating the choose(M,2) overlapping pairs as independent.
    # (B=10: nominal budget 1.8e-4, i.e. at most ~4 of 20000 replicates.)
    b = 10.0
    r_h = math.sqrt(b / (2 * pairs))
    d_u = np.minimum(1.0, a_x + r_h)
    r_d = np.minimum(r_h, np.sqrt(2 * d_u * b / pairs) + 2 * b / (3 * pairs))
    naive = float(np.mean(np.abs(a_x - px**2 / 2) > r_d))
    check(
        "D5 negative control: naive k=choose(M,2) violates its nominal budget",
        naive > 5 * 4 * math.exp(-b),
        f"B=10: observed failure {naive:.5f} vs nominal {4 * math.exp(-b):.5f}",
    )
    # Same B with the correct k=floor(M/2): no violation.
    r_h = math.sqrt(b / (2 * k2))
    d_u = np.minimum(1.0, a_x + r_h)
    r_d = np.minimum(r_h, np.sqrt(2 * d_u * b / k2) + 2 * b / (3 * k2))
    proper = float(np.mean(np.abs(a_x - px**2 / 2) > r_d))
    check(
        "D6 same B with k=floor(M/2) stays inside the budget",
        proper <= 4 * math.exp(-b),
        f"B=10: observed failure {proper:.5f} <= {4 * math.exp(-b):.5f}",
    )


# --------------------------------------------------------------------------- E
def part_e():
    worst_paper, worst_actual, ok = 0.0, 0.0, True
    for n in itertools.chain(range(4, 2000), (10**j for j in range(4, 13))):
        b = (32 * n**4).bit_length()  # smallest-ish integer with 2^B >= 32 N^4
        ok &= 2**b >= 32 * n**4
        events_paper = 3 * n + 3 * n + 3 * n * n + 3 * n * n
        fail_paper = events_paper * 2 * math.exp(-b)
        probes = n - n // 2
        events_actual = 6 * probes + 6 * probes * (probes - 1) // 2
        fail_actual = events_actual * 2 * math.exp(-b)
        worst_paper = max(worst_paper, fail_paper * n * n)
        worst_actual = max(worst_actual, fail_actual * n * n)
        ok &= fail_paper <= n**-2
    check(
        "E1 union bound: (6N+6N^2) 2exp(-B) <= N^-2 whenever 2^B >= 32 N^4",
        ok,
        f"worst N^2*failure: paper count {worst_paper:.3e}, actual probe/pair count {worst_actual:.3e}",
    )
    n = sp.symbols("N", positive=True)
    # with exp(-B) <= 2^-B <= 1/(32 N^4): (6N+6N^2)*2/(32N^4) <= N^-2  <=>  12N+12N^2 <= 32N^2
    expr = sp.simplify(32 * n**2 - 12 * n - 12 * n**2)
    check(
        "E2 symbolic slack 20N^2-12N > 0 for N >= 1",
        sp.simplify(expr - (20 * n**2 - 12 * n)) == 0,
    )


# --------------------------------------------------------------------------- F
def part_f():
    rng = np.random.default_rng(7)

    def g(t):  # all first/second partials have nonnegative coefficients on t>0
        p, f, a = t
        return p * p * f + a**3 + p * f * a + 0.5 * f * f * a

    def grad_sup(hi):
        p, f, a = hi
        return np.array(
            [2 * p * f + f * a, p * p + p * a + f * a, 3 * a * a + p * f + 0.5 * f * f]
        )

    def hess_sup(hi):
        p, f, a = hi
        return np.array(
            [[2 * f, 2 * p + a, f], [2 * p + a, a, p + f], [f, p + f, 6 * a]]
        )

    worst_coarse, worst_fine, trials = 0.0, 0.0, 20000
    for _ in range(trials):
        e_b = rng.uniform(0.0, 0.05, 3)
        e_d = rng.uniform(0.0, 0.02, 3)
        b_hat = rng.uniform(0.3, 0.7, 3)
        d_hat = rng.uniform(0.0, 0.1, 3) * np.array([1, -1, 1])
        b = b_hat + rng.uniform(-1, 1, 3) * e_b
        d = d_hat + rng.uniform(-1, 1, 3) * e_d
        # boxes: B_b = b_hat +- e_b, B_d = d_hat +- e_d, R = B_b + hull(0, B_d)
        hi = b_hat + e_b + np.maximum(0.0, d_hat + e_d)
        big_d = np.maximum(np.abs(d_hat - e_d), np.abs(d_hat + e_d))
        lhs = abs((g(b_hat + d_hat) - g(b_hat)) - (g(b + d) - g(b)))
        l_i, h_ij = grad_sup(hi), hess_sup(hi)
        coarse = l_i.max() * e_d.sum() + h_ij.max() * big_d.sum() * e_b.sum()
        fine = float(l_i @ e_d + e_b @ h_ij @ big_d)
        worst_coarse = max(worst_coarse, lhs / coarse)
        worst_fine = max(worst_fine, lhs / fine)
    check(
        "F1 propagation bound (4), coarse form L e_d + H D e_b",
        worst_coarse <= 1.0,
        f"{trials} random boxes, worst lhs/bound {worst_coarse:.3f}",
    )
    check(
        "F2 propagation bound (4), componentwise form",
        worst_fine <= 1.0,
        f"worst lhs/bound {worst_fine:.3f}",
    )
    # The domain matters: with derivative bounds taken on B_b alone the bound can fail.
    bad = 0
    for _ in range(trials):
        e_b = rng.uniform(0.0, 0.01, 3)
        e_d = rng.uniform(0.0, 0.002, 3)
        b_hat = rng.uniform(0.3, 0.7, 3)
        d_hat = rng.uniform(0.2, 0.3, 3)
        b = b_hat + rng.choice([-1.0, 1.0], 3) * e_b
        d = d_hat + rng.choice([-1.0, 1.0], 3) * e_d
        hi_small = b_hat + e_b
        big_d = d_hat + e_d
        lhs = abs((g(b_hat + d_hat) - g(b_hat)) - (g(b + d) - g(b)))
        fine = float(grad_sup(hi_small) @ e_d + e_b @ hess_sup(hi_small) @ big_d)
        bad += lhs > fine
    check(
        "F3 negative control: bounds on B_b alone (not R) are violated",
        bad > 0,
        f"{bad}/{trials} violations when the hull(0,B_d) part of R is dropped",
    )


# --------------------------------------------------------------------------- G
def part_g():
    lam, rho, lg = sp.symbols("lambda rho L", positive=True)
    n = rho * lg / lam**3  # N lambda^3 / log N = rho
    t1 = sp.simplify(sp.sqrt(lam * lg / n) / lam**2)
    t2 = sp.simplify((lg / n) / lam**2)
    t3 = sp.simplify(lam * sp.sqrt(lg / n) / lam**2)
    ok = (
        sp.simplify(t1 - 1 / sp.sqrt(rho)) == 0
        and sp.simplify(t2 - lam / rho) == 0
        and sp.simplify(t3 - sp.sqrt(lam / rho)) == 0
    )
    check(
        "G1 width (8) / lambda^2 = rho^-1/2 + lambda/rho + sqrt(lambda/rho), rho=N lambda^3/log N",
        ok,
        f"terms: {t1}, {t2}, {t3}",
    )
    nn, c1, c2 = sp.symbols("N c C", positive=True)
    cheb = sp.simplify(
        c2 * (nn**2 * lam**6 + nn**3 * lam**12) / (c1 * nn**2 * lam**6) ** 2
    )
    target = c2 / c1**2 * (1 / (nn**2 * lam**6) + 1 / nn)
    check(
        "G2 Chebyshev: Var W/(E W)^2 <= (C/c^2)(1/(N lambda^3)^2 + 1/N)",
        sp.simplify(cheb - target) == 0,
    )
    check(
        "G3 pair-mass exponent 1 + 2*(3/2) + 2 = 6", 1 + 2 * sp.Rational(3, 2) + 2 == 6
    )
    # action budget: |dz|^2/(2 tau lambda f) with |dz| ~ lambda^(3/2) is lambda^2;
    # equal-z pair: |dy| ~ lambda^2 gives lambda^4/lambda = lambda^3.
    check(
        "G4 action exponents: (3/2)*2-1 = 2 and 2*2-1 = 3",
        sp.Rational(3, 2) * 2 - 1 == 2 and 2 * 2 - 1 == 3,
    )
    # Separation sweep of local_null_results.md: h = lambda^x, signal d = lambda*min(h,lambda).
    best, best_x = None, None
    for x in np.linspace(0.0, 3.0, 301):
        d_exp = 1 + max(x, 1.0)
        witness = x + d_exp  # N >> 1/(h d)
        width = 2 * d_exp - x  # N >> h/d^2
        price = max(witness, width)
        if best is None or price < best - 1e-12:
            best, best_x = price, x
    check(
        "G5 separation sweep: max(1/(h d), h/d^2) is minimised at h=lambda with exponent 3",
        abs(best - 3.0) < 1e-9 and abs(best_x - 1.0) < 1e-9,
        f"min exponent {best:.3f} at h=lambda^{best_x:.2f}",
    )
    # log N bookkeeping: B ~ 4 log2 N; r_H^2 = B/(2k) = o(lambda^2) and d_U = O(lambda).
    worst = 0.0
    for expo in (6, 9, 12, 18, 24, 36):
        for rho_val in (1e2, 1e4):
            n_val = 10.0**expo
            lam_val = (rho_val * math.log(n_val) / n_val) ** (1 / 3)
            b_val = math.ceil(math.log2(32 * n_val**4))
            k_val = n_val / 4
            worst = max(worst, math.sqrt(b_val / (2 * k_val)) / lam_val)
    check(
        "G6 initial Hoeffding radius / lambda -> 0 under N lambda^3/log N >= 100",
        worst < 0.2,
        f"largest r_H/lambda on the grid {worst:.3e}",
    )


if __name__ == "__main__":
    part_a()
    part_bc()
    part_d()
    part_e()
    part_f()
    part_g()
    print(f"SUMMARY s2_confidence_bounds: {len(FAILS)} FAIL")
    sys.exit(1 if FAILS else 0)
