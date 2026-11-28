# ruff: noqa: E501, UP031
"""Referee arithmetic for Section 7 of gram_fourier_noise.md (block / majority /
grid scan) and for the last paragraph of Section 3 (domination on the scan
range).  Equations (5), (10), (11) are taken as hypotheses here.

L1  exact binomial tail: P(Bin(L, p) > L/2) <= exp(-L/8) for p <= 1/4.
L2  L = ceil(16 log N) gives exp(-L/8) <= N^-2.
L3  the grid has fewer than 4 sqrt(N) points; whole-scan size bound -> 0.
L4  on 1 <= omega <= n^(K/(2K+2)) the right side of (5) is at most K times
    the bracket of the threshold (12); remainder and 1/n are <= first term.
L4b the same domination up to sqrt(n) (referee remark on the scan cap).
L5  exponent bookkeeping: which requirement on alpha is binding.
L6  the same with the log factors and the block splitting kept, in logs.

Deterministic; prints PASS/FAIL; exits nonzero on any FAIL.
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction

import sympy as sp

FAILS = []


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL") + " " + name + (" :: " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def l1_majority():
    worst = 0.0
    for prob in (Fraction(1, 4), Fraction(1, 5), Fraction(1, 10)):
        for blocks in range(1, 401):
            tail = sum(
                math.comb(blocks, k) * prob**k * (1 - prob) ** (blocks - k)
                for k in range(blocks // 2 + 1, blocks + 1)
            )
            worst = max(worst, float(tail) / math.exp(-blocks / 8))
    check(
        "L1 exact P(Bin(L,p) > L/2) <= exp(-L/8) for p in {1/4,1/5,1/10}, L = 1..400",
        worst <= 1.0,
        f"largest tail/exp(-L/8) = {worst:.4f}",
    )


def l2_l3_union():
    worst_major, worst_grid, sizes = -math.inf, -math.inf, []
    for exponent in range(2, 61):
        log_n_total = exponent * math.log(10.0)
        blocks = math.ceil(16 * log_n_total)
        worst_major = max(
            worst_major, -blocks / 8 + 2 * log_n_total
        )  # log(exp(-L/8) N^2)
        for k in (1, 2, 3, 5, 8, 50):
            # n <= N, so the number of grid points is at most 4 n^(K/(2K+2)) < 4 sqrt(N)
            log_points = math.log(4.0) + k / (2 * k + 2) * log_n_total
            worst_grid = max(worst_grid, log_points - math.log(4.0) - 0.5 * log_n_total)
        sizes.append(math.log(4.0) + 0.5 * log_n_total - blocks / 8)
    check(
        "L2 L=ceil(16 log N) gives exp(-L/8) <= N^-2, N = 1e2..1e60",
        worst_major <= 1e-12,
        f"max log(exp(-L/8) N^2) = {worst_major:.3f}",
    )
    check(
        "L3 grid points < 4 sqrt(N) and whole-scan bound 4 sqrt(N) exp(-L/8) decreases to 0",
        worst_grid < 0
        and all(b < a for a, b in zip(sizes, sizes[1:], strict=False))
        and sizes[-1] < -100,
        f"log(whole-scan bound) from {sizes[0]:.1f} (N=1e2) to {sizes[-1]:.1f} (N=1e60)",
    )


def l4_domination():
    worst_ratio, worst_first, ok = 0.0, 0.0, True
    for k in range(1, 9):
        for exponent in range(1, 13):
            n = 10.0**exponent
            cap = n ** (k / (2 * k + 2))
            omegas = [1.0, cap] + [
                math.exp(math.log(cap) * j / 200) for j in range(1, 200)
            ]
            for w in omegas:
                first = 1 / (w * math.sqrt(n))
                series = sum(w ** (j - 2) * n ** (-j / 2) for j in range(1, k + 1))
                rem = w**k / n ** ((k + 1) / 2)
                bound5 = series + 1 / n + rem
                bracket12 = first + 1 / n + rem
                worst_ratio = max(worst_ratio, bound5 / bracket12 / k)
                worst_first = max(worst_first, max(rem, 1 / n, series / k) / first)
                ok &= bound5 <= k * bracket12 * (1 + 1e-12)
    check(
        "L4 on the scan range: RHS(5) <= K * bracket(12); remainder, 1/n and series/K are <= 1/(omega sqrt n)",
        ok and worst_first <= 1 + 1e-9,
        f"max RHS(5)/(K bracket) = {worst_ratio:.4f}, max term/first = {worst_first:.6f} (K=1..8, n=1e1..1e12)",
    )


def l4b_wider_grid():
    """Referee remark: the same inequality RHS(5) <= K*bracket(12) holds on the
    WIDER range 1 <= omega <= sqrt(n), so the null control of Section 7 would
    survive a grid up to sqrt(n) (still fewer than 4 sqrt(N) points)."""
    ok, worst = True, 0.0
    for k in range(1, 9):
        for exponent in range(1, 13):
            n = 10.0**exponent
            cap = math.sqrt(n)
            for j in range(0, 201):
                w = math.exp(math.log(cap) * j / 200)
                series = sum(w ** (i - 2) * n ** (-i / 2) for i in range(1, k + 1))
                rem = w**k / n ** ((k + 1) / 2)
                ratio = (series + 1 / n + rem) / (
                    k * (1 / (w * math.sqrt(n)) + 1 / n + rem)
                )
                worst = max(worst, ratio)
                ok &= ratio <= 1 + 1e-12
    check(
        "L4b the same domination holds up to omega = sqrt(n): a grid capped at sqrt(n) keeps the null control",
        ok,
        f"max RHS(5)/(K bracket) = {worst:.4f}; with that grid the binding alpha would be (2K+4)/(K+1) = 2+2/(K+1)",
    )


def l5_exponents():
    # N ~ lambda^-alpha, n ~ N up to logs, omega = 1/lambda, signal ~ lambda^2.
    # signal/term diverges iff alpha exceeds:
    #   lambda/sqrt(n)            : 2
    #   1/n                       : 2
    #   omega^K n^(-(K+1)/2)      : (2K+4)/(K+1)
    #   1/lambda <= n^(K/(2K+2))  : (2K+2)/K        (scan range)
    k = sp.symbols("K", positive=True)
    alpha = sp.symbols("alpha", positive=True)
    lam = sp.symbols("lambda", positive=True)
    n = lam**-alpha

    def exponent(
        expr,
    ):  # expr = lambda^e ; the quantity diverges as lambda -> 0 iff e < 0
        return sp.simplify(lam * sp.diff(sp.log(expr), lam))

    e1 = exponent(lam**2 / (lam / sp.sqrt(n)))
    e2 = exponent(lam**2 * n)
    e3 = exponent(lam**2 * lam**k * n ** ((k + 1) / 2))
    e4 = exponent(lam * n ** (k / (2 * k + 2)))
    roots = [sp.solve(sp.Eq(e, 0), alpha)[0] for e in (e1, e2, e3, e4)]
    expected = [2, 2, (2 * k + 4) / (k + 1), (2 * k + 2) / k]
    same = all(sp.simplify(r - x) == 0 for r, x in zip(roots, expected, strict=False))
    gap = sp.simplify(roots[3] - roots[2])
    check(
        "L5 critical alpha per requirement: 2, 2, (2K+4)/(K+1), (2K+2)/K; scan range is binding by 2/(K(K+1))",
        same and sp.simplify(gap - 2 / (k * (k + 1))) == 0,
        f"thresholds {[str(sp.simplify(r)) for r in roots]} ; range minus remainder = {gap}",
    )
    table = []
    for kk in (1, 2, 3, 4, 8, 16):
        table.append(
            "K=%d: range %s, remainder %s"
            % (kk, Fraction(2 * kk + 2, kk), Fraction(2 * kk + 4, kk + 1))
        )
    print("     " + "; ".join(table))


def scan_state(k, alpha, log_inv_lambda):
    """All quantities in natural logs; n = (N/L - 3 ceil(log(N/L)))/3."""
    log_total = alpha * log_inv_lambda
    blocks = math.ceil(16 * log_total)
    log_block = log_total - math.log(blocks)
    reserve = 3 * math.ceil(log_block) / math.exp(min(log_block, 700.0))
    if reserve >= 1:
        return False, -math.inf, log_total, -math.inf  # block too small: vote is zero
    log_n = log_block - math.log(3.0) + math.log1p(-reserve)
    log_cap = k / (2 * k + 2) * log_n
    in_range = log_inv_lambda < log_cap
    log_w = log_inv_lambda
    terms = [-log_w - 0.5 * log_n, -log_n, k * log_w - (k + 1) / 2 * log_n]
    log_threshold = (
        math.log(log_n)
        + max(terms)
        + math.log(sum(math.exp(t - max(terms)) for t in terms))
    )
    return in_range, -2 * log_inv_lambda - log_threshold, log_total, log_n


def l6_with_logs():
    grid = [j * 0.25 for j in range(8, 4001)]  # log10(1/lambda) from 2 to 1000
    report = []
    ok_above, ok_below = True, True
    for k in (1, 2, 4, 8):
        critical = 2 + 2 / k
        for delta in (0.1, 0.5, 1.0):
            first = None
            for log10 in grid:
                state = scan_state(k, critical + delta, log10 * math.log(10.0))
                if state[0]:
                    first = (log10, state)
                    break
            if first is None:
                ok_above = False
                continue
            log10, state = first
            later = [
                scan_state(k, critical + delta, x * math.log(10.0))
                for x in (log10 + 5, log10 + 20, log10 + 100)
            ]
            ok_above &= all(s[0] for s in later)
            # asymptotic claim only: the ratio must grow and end above 1; it may be < 1 at entry
            ok_above &= (
                later[0][1] > state[1]
                and later[1][1] > later[0][1]
                and later[2][1] > later[1][1]
                and later[2][1] > 0
            )
            beaten = next(
                x
                for x in grid
                if x >= log10
                and scan_state(k, critical + delta, x * math.log(10.0))[1] > 0
            )
            report.append(
                "K=%d alpha=%.3f: 1/lambda enters the scan at lambda=1e-%.2f (N=1e%.1f, n=1e%.1f) with lambda^2/t_n=1e%.1f;"
                " lambda^2 > t_n from lambda=1e-%.2f on"
                % (
                    k,
                    critical + delta,
                    log10,
                    state[2] / math.log(10.0),
                    state[3] / math.log(10.0),
                    state[1] / math.log(10.0),
                    beaten,
                )
            )
        for delta in (0.0, -0.05):
            inside = [
                scan_state(k, critical + delta, x * math.log(10.0))[0] for x in grid
            ]
            ok_below &= not any(inside)
    for line in report:
        print("     " + line)
    check(
        "L6 with logs and block splitting: for alpha > 2+2/K, 1/lambda eventually stays inside the scan and lambda^2/t_n grows",
        ok_above,
        "%d (K, alpha) cases, see lines above" % len(report),
    )
    check(
        "L6 at alpha = 2+2/K and just below, 1/lambda is never inside the scan (lambda from 1e-2 to 1e-1000)",
        ok_below,
        "the condition delivered is the strict inequality alpha > 2+2/K",
    )
    # what the scan cap costs: with the cap alone removed, alpha > (2K+4)/(K+1) would do
    wins = []
    for k in (1, 2, 4):
        a = 0.5 * (
            (2 * k + 4) / (k + 1) + (2 * k + 2) / k
        )  # strictly between the two thresholds
        states = [scan_state(k, a, x * math.log(10.0)) for x in (50.0, 200.0, 800.0)]
        wins.append(
            all(not s[0] for s in states)
            and states[2][1] > states[1][1] > states[0][1]
            and states[2][1] > 0
        )
    check(
        "L6 between (2K+4)/(K+1) and (2K+2)/K the threshold IS beaten by lambda^2 but 1/lambda is outside the scan",
        all(wins),
        "so every term of (12) is beaten on the stated range and the scan cap is the only binding requirement",
    )


if __name__ == "__main__":
    l1_majority()
    l2_l3_union()
    l4_domination()
    l4b_wider_grid()
    l5_exponents()
    l6_with_logs()
    print("SUMMARY failures=%d" % len(FAILS))
    sys.exit(1 if FAILS else 0)
