"""Orchestrator-level arithmetic checks for the 2026-11-24 audit.

Exact (SymPy) checks of the exponent bookkeeping in the three audited
arguments, of the sizing numbers in the design draft D1, and of the
heuristic table in next_question.md. This script checks arithmetic only.
It does not check any analytic lemma; those are the referees' scripts in
the sibling folders.

Run with the system Python that has SymPy:

    python docs/paper/paper_a/order_wave_energy/independent_audit/orchestrator_checks.py
"""

from __future__ import annotations

import math
import sys

import sympy as sp

FAILS = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global FAILS
    print(f"{'PASS' if ok else 'FAIL'}  {name}{('  ' + detail) if detail else ''}")
    if not ok:
        FAILS += 1


def proof1_exponents() -> None:
    """quantized_order_bound.md eqs. (7)-(8) and the 4/3 limit."""
    N, d, K = sp.symbols("N delta K", positive=True)
    h = (d / N**K) ** (sp.Rational(2) / (2 * K + 1))
    # The two terms that the choice of h balances.
    tree = N ** (K + 1) * h**K
    hell = N * d / sp.sqrt(h)
    target = (N ** (3 * K + 1) * d ** (2 * K)) ** (1 / (2 * K + 1))
    for k in (2, 3, 4, 7):
        sub = {K: k}
        r1 = sp.simplify(
            sp.powsimp(
                sp.expand_power_base((tree / target).subs(sub), force=True), force=True
            )
        )
        r2 = sp.simplify(
            sp.powsimp(
                sp.expand_power_base((hell / target).subs(sub), force=True), force=True
            )
        )
        check(
            f"proof1 K={k}: tree term equals (N^(3K+1) delta^(2K))^(1/(2K+1))",
            r1 == 1,
            str(r1),
        )
        check(f"proof1 K={k}: Hellinger term equals the same", r2 == 1, str(r2))
        nh = sp.simplify(
            sp.powsimp(
                sp.expand_power_base(
                    ((N * h) / (N * d**2) ** (sp.Rational(1) / (2 * k + 1))).subs(sub),
                    force=True,
                ),
                force=True,
            )
        )
        check(f"proof1 K={k}: N h equals (N delta^2)^(1/(2K+1))", nh == 1, str(nh))
    # With delta = lambda^2 and N = lambda^-alpha the bound vanishes iff
    # alpha < 4K/(3K+1), which increases to 4/3.
    k = sp.symbols("k", positive=True)
    crit = 4 * k / (3 * k + 1)
    check(
        "proof1: critical exponent 4K/(3K+1) -> 4/3",
        sp.limit(crit, k, sp.oo) == sp.Rational(4, 3),
    )
    check(
        "proof1: 4K/(3K+1) is increasing in K",
        sp.simplify(sp.diff(crit, k)) == 4 / (3 * k + 1) ** 2,
    )
    vals = {2: sp.Rational(8, 7), 3: sp.Rational(6, 5), 4: sp.Rational(16, 13)}
    check(
        "proof1: K=2,3,4 give 8/7, 6/5, 16/13",
        all(crit.subs(k, kk) == v for kk, v in vals.items()),
    )
    check(
        "proof1: K=1 gives exponent 1 (no gain over the coupling bound)",
        crit.subs(k, 1) == 1,
    )
    # Side conditions at N = lambda^-alpha, alpha < 4K/(3K+1): log_lambda of
    # h and of h/delta. h -> 0 needs a positive exponent, h/delta -> infinity
    # a negative one.
    a = sp.symbols("alpha", positive=True)
    for kk in (2, 3, 4, 7):
        e_h = sp.Rational(2, 2 * kk + 1) * (2 + a * kk)
        e_ratio = e_h - 2
        amax = crit.subs(k, kk)
        check(
            f"proof1 K={kk}: h -> 0 for all alpha>0",
            sp.simplify(e_h.subs(a, sp.Rational(1, 100))) > 0,
        )
        check(
            f"proof1 K={kk}: h/delta -> infinity on alpha<4K/(3K+1)",
            sp.simplify(e_ratio.subs(a, amax)) < 0,
            f"exponent at the edge {sp.simplify(e_ratio.subs(a, amax))}",
        )
        e_nh = e_h - a
        check(
            f"proof1 K={kk}: N h -> 0 for alpha<4", sp.solve(sp.Eq(e_nh, 0), a) == [4]
        )


def proof2_exponents() -> None:
    """local_null_challenge.md eqs. (6)-(8)."""
    a = sp.symbols("alpha", positive=True)
    # N = lambda^-alpha, exponents of lambda (logs suppressed).
    pair_count = -2 * a + 6  # N^2 lambda^6
    width = sp.Rational(1, 2) * (1 + a)  # sqrt(lambda/N)
    check(
        "proof2: E W = N^2 lambda^6 diverges iff alpha > 3",
        sp.solve(sp.Eq(pair_count, 0), a) == [3],
    )
    check(
        "proof2: sqrt(lambda/N) << lambda^2 iff alpha > 3",
        sp.solve(sp.Eq(width, 2), a) == [3],
    )
    # The other two width terms are smaller on alpha > 3.
    check("proof2: 1/N << lambda^2 iff alpha > 2", sp.solve(sp.Eq(a, 2), a) == [2])
    check(
        "proof2: lambda sqrt(1/N) << lambda^2 iff alpha > 2",
        sp.solve(sp.Eq(1 + a / 2, 2), a) == [2],
    )
    # Variance ratio Var W / (E W)^2 with Var <= N^2 lambda^6 + N^3 lambda^12.
    check(
        "proof2: N^3 lambda^12 / (N^2 lambda^6)^2 = 1/N",
        sp.simplify((3 * -a + 12) - 2 * (-2 * a + 6)) == a,
    )
    # Separation sweep in local_null_results.md: the two prices cross at h = lambda.
    lam, h = sp.symbols("lambda h", positive=True)
    below = (1 / (lam * h**2), 1 / (lam**2 * h))  # h <= lambda, d = lambda h
    above = (1 / (lam**2 * h), h / lam**4)  # h >= lambda, d = lambda^2
    check(
        "proof2: both regimes give lambda^-3 at h = lambda",
        all(sp.simplify(e.subs(h, lam) * lam**3) == 1 for e in below + above),
    )


def proof3_exponents() -> None:
    """gram_fourier_noise.md eqs. (5), (10)-(12) and the scan cap."""
    a, K = sp.symbols("alpha K", positive=True)
    # exponents of lambda at omega = 1/lambda, n = lambda^-alpha; signal is 2.
    first = 1 + a / 2  # 1/(omega sqrt n)
    probe = a  # 1/n
    rem = -K + a * (K + 1) / 2  # omega^K n^-(K+1)/2
    check(
        "proof3: 1/(omega sqrt n) << lambda^2 iff alpha > 2",
        sp.solve(sp.Eq(first, 2), a) == [2],
    )
    check("proof3: 1/n << lambda^2 iff alpha > 2", sp.solve(sp.Eq(probe, 2), a) == [2])
    sol = sp.solve(sp.Eq(rem, 2), a)
    check(
        "proof3: remainder << lambda^2 iff alpha > 2 + 2/(K+1)",
        sp.simplify(sol[0] - (2 + 2 / (K + 1))) == 0,
        str(sol),
    )
    cap = sp.solve(sp.Eq(a * K / (2 * K + 2), 1), a)
    check(
        "proof3: 1/lambda inside the scan cap iff alpha > 2 + 2/K",
        sp.simplify(cap[0] - (2 + 2 / K)) == 0,
        str(cap),
    )
    check(
        "proof3: the cap condition is the binding one",
        sp.simplify((2 + 2 / K) - (2 + 2 / (K + 1))) == 2 / (K * (K + 1)),
    )
    # Middle Taylor terms omega^(k-2) n^(-k/2), k = 2..K: on omega <= sqrt(n)
    # each is at most 1/n.
    kk = sp.symbols("k", positive=True)
    mid = -(kk - 2) + a * kk / 2
    check(
        "proof3: order-k Taylor term << lambda^2 iff alpha > 2",
        sp.solve(sp.Eq(mid, 2), a) == [2],
    )
    # Majority vote: L = ceil(16 log N), per-block vote probability <= 1/4.
    for n in (10**3, 10**6, 10**12):
        blocks = math.ceil(16 * math.log(n))
        hoeffding = math.exp(-2 * blocks * (0.25) ** 2)
        check(
            f"proof3: exp(-L/8) <= N^-2 at N={n:.0e}",
            hoeffding <= n**-2.0 * (1 + 1e-12),
            f"{hoeffding:.3e} vs {n**-2.0:.3e}",
        )
    # Elementary check: exact variance for H_hat = F_n(X) - X at omega = 2 pi j.
    x, y = sp.symbols("x y", real=True)
    j = sp.symbols("j", integer=True, positive=True)
    w = 2 * sp.pi * j
    # covariance of the uniform empirical process: min(x,y) - x y, per reference.
    cov = sp.integrate(
        sp.integrate(
            (sp.Min(x, y) - x * y) * sp.exp(-sp.I * w * x) * sp.exp(sp.I * w * y),
            (y, 0, 1),
        ),
        (x, 0, 1),
    )
    cov = sp.simplify(cov)
    check(
        "proof3 toy: Fourier coefficient of the Brownian-bridge kernel is 1/omega^2",
        sp.simplify(cov - 1 / w**2) == 0,
        str(cov),
    )
    diag = sp.integrate(x * (1 - x), (x, 0, 1))
    check("proof3 toy: integral of x(1-x) is 1/6", diag == sp.Rational(1, 6))
    # Var S = (1/n)[(1-1/m)/omega^2 + (1/m)(1/6)] follows from the two lines above.


def design_d1_sizing() -> None:
    """next_question.md, draft D1: 95% two-sided, 90% power."""
    z95, z90 = 1.959964, 1.281552
    mult = z95 + z90
    check(
        "D1: 90% power multiplier is 3.2415", abs(mult - 3.2415) < 1e-3, f"{mult:.4f}"
    )
    grid, span = 5, math.log(16.0)
    sxx = span**2 * grid * (grid + 1) / (12 * (grid - 1))
    t, bias = 0.15, 0.05
    se_max = (t - bias) / mult
    reps = (0.5 / (se_max * math.sqrt(sxx))) ** 2
    check("D1: SE cap 0.0308", abs(se_max - 0.0308) < 5e-4, f"{se_max:.4f}")
    check("D1: 90% power needs R >= 55", 54.0 < reps <= 55.0, f"R = {reps:.1f}")
    reps_even = (0.5 / ((t - bias) / z95 * math.sqrt(sxx))) ** 2
    check(
        "D1: break-even is about R = 20 for the same allowance",
        19.0 < reps_even < 21.0,
        f"R = {reps_even:.1f}",
    )
    reps_even0 = (0.5 / (t / z95 * math.sqrt(sxx))) ** 2
    check(
        "D1: break-even with no bias allowance is about R = 9",
        8.0 < reps_even0 < 10.0,
        f"R = {reps_even0:.1f}",
    )
    check(
        "D1: power/break-even ratio is 2.74 at fixed margin",
        abs((mult / z95) ** 2 - 2.735) < 0.01,
        f"{(mult / z95) ** 2:.3f}",
    )


def heuristic_table() -> None:
    """Exponents quoted in next_question.md Section 3 (arithmetic only)."""
    a = sp.symbols("alpha", positive=True)
    # log_lambda of each information count with N = lambda^-alpha, delta = lambda^2;
    # the count diverges when its exponent is negative.
    rows = {
        "cell-limited, phase known: N (N delta)^2": (
            -a + 2 * (-a + 2),
            sp.Rational(4, 3),
        ),
        "cell-limited, d_c = N lambda^2: N^2 d_c (N delta)^4": (
            -2 * a + (-a + 2) + 4 * (-a + 2),
            sp.Rational(10, 7),
        ),
        "cell-limited, per event: N (N delta)^4": (
            -a + 4 * (-a + 2),
            sp.Rational(8, 5),
        ),
        "count-limited, phase known: N^(3/2) (sqrt(N) delta)^2": (
            -sp.Rational(3, 2) * a + 2 * (-a / 2 + 2),
            sp.Rational(8, 5),
        ),
        "count-limited, per event: sqrt(N) (N^(3/4) delta)^2": (
            -a / 2 + 2 * (-sp.Rational(3, 4) * a + 2),
            sp.Integer(2),
        ),
    }
    for name, (expo, want) in rows.items():
        got = sp.solve(sp.Eq(expo, 0), a)
        check(f"heuristic: {name} -> alpha = {want}", got == [want], str(got))
    # windowed filter: 1/(N lambda^2) << N lambda^2 iff N lambda^2 >> 1 iff alpha > 2.
    check(
        "heuristic: window exists iff alpha > 2",
        sp.solve(sp.Eq(2 * (-a + 2), 0), a) == [2],
    )
    # discreteness-length translation: lambda = l^(4/alpha).
    check(
        "units: alpha=2 is lambda ~ l^2, alpha=4/3 is lambda ~ l^3",
        (4 / sp.Integer(2), 4 / sp.Rational(4, 3)) == (2, 3),
    )


def main() -> int:
    proof1_exponents()
    proof2_exponents()
    proof3_exponents()
    design_d1_sizing()
    heuristic_table()
    print(f"\n{'ALL PASS' if FAILS == 0 else f'{FAILS} FAIL'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
