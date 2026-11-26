# ruff: noqa: E501, UP031
"""Referee check of Sections 3-4 of local_null_challenge.md at the level of EXACT
population order moments (independent implementation; see audit_geometry.py).

For a constructed pair X < Y (exact wave relation) the observer-side calibration
F(p,f,a) = V(p,f,a/p^2;q) is fed the exact wave population moments

    p = Pr(Z<X), f = Pr(X<Z), a = Pr(Z1<X, Z2<X, Z1~Z2) = p^2 c,

computed here from the referee's own quadrature of the Jacobi-kernel hierarchy.

a  Background self-test: quadrature moments + calibration return v.
b  Quadrature refinement leaves the result unchanged.
c  Independent replication of the author's frozen 12-row table and r=2048 witness.
e0 Uniform pointwise kernel bound that makes the imported bias expansion uniform.
d  Independent witness grid (other u, tau, direction, phase than the author's):
   exact margin > 0 and F(theta_Y) < F(theta_X).
e  Bias F(theta_X) - v_tilde(X) = O(lambda^4 log(1/lambda)) (ratio test over rungs).
f  |theta_Y - theta_X| = O(lambda), with the nestedness signs.
INFO lines: finite-lambda failures and the size of the constants.

Deterministic.  Prints PASS/FAIL lines; exits nonzero on any FAIL.
Run:  python s4_population_witness.py        (system Python: numpy only; ~1.5 min)
"""

from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_geometry import A_TIP, Wave, calibration_f, sym  # noqa: E402

FAILS = []
SIGNS = np.array([1.0, -1.0, 1.0])


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def leading(w, u, tau):
    m = float(w.m(A_TIP, u))
    return (
        2j
        * w.k
        / math.cos(w.k * u) ** 2
        / (m * (w.Mtot - m))
        * cmath.exp(1j * (u / w.lam + w.phase))
        * (cmath.exp(1j * tau) - 1)
    )


def part_a():
    worst = 0.0
    rng = np.random.default_rng(3)
    for q in (1.0, 0.25):
        bg = Wave(q, 1.0 / 64, background=True)
        gamma = bg.gamma()
        for _ in range(6):
            u = rng.uniform(-0.3, 0.3)
            z = rng.uniform(-0.15, 0.15, 2)
            vt = rng.uniform(0.4, 0.6)
            point = bg.from_canonical(u, vt, z)
            theta = bg.moments(point, gamma=gamma)
            worst = max(
                worst, abs(calibration_f(*theta, q) - point[1]), abs(point[1] - vt)
            )
    check(
        "a  background: quadrature moments + calibration F return the null coordinate v",
        worst < 1e-11,
        f"max |F - v| = {worst:.2e}",
    )


def part_b():
    q, r = 1.0, 256
    w = Wave(q, 1.0 / r)
    point = w.from_canonical(0.2, 0.5, np.array([0.12, 0.16]))
    one = w.moments(point, factor=1.0)
    two = w.moments(point, factor=2.0)
    diff = float(np.max(np.abs(one - two)))
    check(
        "b  doubling the quadrature order changes the wave moments by < 1e-13",
        diff < 1e-13,
        f"max change {diff:.2e}",
    )


PUBLISHED = {  # local_null_results.md: decrease / lambda^2 at phase 0.37, u1=1/5, tau=pi
    0.25: {
        64: -0.271854,
        128: 0.016543,
        256: 0.075523,
        512: 0.090478,
        1024: 0.092035,
        2048: 0.091750,
    },
    1.0: {
        64: -0.179415,
        128: 0.092361,
        256: 0.147062,
        512: 0.160516,
        1024: 0.162671,
        2048: 0.162087,
    },
}


def frozen_pair(q, r):
    """The author's frozen construction, rebuilt from the design text with referee code."""
    w = Wave(q, 1.0 / r)  # phase 0.37 as in the frozen design
    u1, u2 = 0.2, 0.2 + math.pi / r
    diff = sym(w.f_matrix_complex(u1) - w.f_matrix_complex(u2))
    _, vectors = np.linalg.eigh(diff)
    z = 0.2 * vectors[:, -1]
    first = w.from_canonical(u1, 0.5, z)
    second = w.from_canonical(u2, 0.5, z)
    initial = w.margin(first, second)
    second = w.from_canonical(u2, 0.5 - max(initial, 0.0) / 2, z)
    gamma = w.gamma()
    th1, th2 = w.moments(first, gamma=gamma), w.moments(second, gamma=gamma)
    decrease = calibration_f(*th1, q) - calibration_f(*th2, q)
    return w, first, second, th1, th2, decrease, initial


def part_c():
    worst, margins_ok, keep = 0.0, True, None
    lines = []
    for q in (0.25, 1.0):
        for r in (64, 128, 256, 512, 1024, 2048):
            out = frozen_pair(q, r)
            scaled = out[5] * r * r
            worst = max(worst, abs(scaled - PUBLISHED[q][r]))
            margins_ok &= out[6] > 0
            lines.append(f"{scaled:.6f}")
            if (q, r) == (1.0, 2048):
                keep = out
    check(
        "c1 independent replication of the 12-row frozen table (decrease / lambda^2)",
        worst < 1e-5 and margins_ok,
        f"max |referee - published| = {worst:.2e}; all 12 exact margins positive: {margins_ok}; "
        "referee values: " + ", ".join(lines),
    )
    w, first, second, th1, th2, decrease, _ = keep
    margin = w.margin(first, second)
    inc = (th2 - th1) * SIGNS * 2048
    published_margin, published_decrease = 3.86450424167e-8, 3.86446596945e-8
    published_inc = np.array([1.035476, 0.386949, 0.021323])
    check(
        "c2 archived q=1, r=2048 witness: margin, decrease and increments/lambda",
        abs(margin / published_margin - 1) < 1e-6
        and abs(decrease / published_decrease - 1) < 1e-5
        and float(np.max(np.abs(inc - published_inc))) < 2e-6,
        f"margin {margin:.11e} (published {published_margin:.11e}); "
        f"F decrease {decrease:.11e} (published {published_decrease:.11e}); "
        f"increments/lambda {np.array2string(inc, precision=6)}",
    )
    return w, first, second, th1, th2, decrease


def part_e0():
    """Pointwise kernel bound behind the uniform O(lambda^4 log 1/lambda) bias:
    |Kbar_lambda - K_0|(s,t) <= C lambda^2 min(l,lambda)^2 / l, l = t - s, uniformly."""
    for q in (1.0, 0.25):
        sup = {}
        for r in (64, 256, 1024, 4096):
            w = Wave(q, 1.0 / r)
            lam = w.lam
            best = 0.0
            for s0 in np.linspace(-0.5, 0.45, 20):
                ell = np.geomspace(1e-6, 0.5 - s0, 400)
                t = s0 + ell
                m = w.m(s0, t)
                x = (w.eps * np.abs(w.osc(np.full_like(t, s0), t)) / m) ** 2
                k0 = np.sin(w.k * ell) / w.k
                delta = k0 * x / (1 + np.sqrt(1 - x))  # = K_0 - Kbar_lambda >= 0
                best = max(
                    best,
                    float(np.max(delta * ell / (lam**2 * np.minimum(ell, lam) ** 2))),
                )
            sup[r] = best
        vals = list(sup.values())
        check(
            f"e0 q={q}: |Kbar_lambda - K_0| <= C lambda^2 min(l,lambda)^2/l with C independent of lambda",
            max(vals) < 1.3 * min(vals) and max(vals) < 40 * q + 5,
            "sup of the ratio at r=64,256,1024,4096: "
            + ", ".join(f"{v:.3f}" for v in vals),
        )


def part_def():
    rungs = (64, 128, 256, 512, 1024)
    u_values, taus = (-0.15, 0.05, 0.25), (2.0, math.pi, 4.2)
    for q in (1.0, 0.25):
        table = {}
        bias = {}
        for r in rungs:
            w = Wave(q, 1.0 / r, phase=1.1)
            lam = w.lam
            gamma = w.gamma()
            for u1 in u_values:
                seg1 = w.segment_data(u1)
                probe = w.from_canonical(u1, 0.45, np.array([0.12, -0.16]))
                bias[(r, u1)] = calibration_f(
                    *w.moments(probe, seg1, gamma), q
                ) - w.v_tilde(probe)
                for tau in taus:
                    u2 = u1 + tau * lam
                    seg2 = w.segment_data(u2)
                    lead = leading(w, u1, tau)
                    # 20 degrees away from the optimal direction
                    psi = 0.5 * (cmath.phase(lead) + math.pi) + math.radians(20)
                    z = 0.2 * np.array([math.cos(psi), math.sin(psi)])
                    first = w.from_canonical(u1, 0.45, z)
                    second = w.from_canonical(u2, 0.45, z)
                    initial = w.margin(first, second)
                    second = w.from_canonical(u2, 0.45 - 0.5 * max(initial, 0.0), z)
                    th1 = w.moments(first, seg1, gamma)
                    th2 = w.moments(second, seg2, gamma)
                    table[(r, u1, tau)] = (
                        initial / lam**2,
                        w.margin(first, second) / lam**2,
                        (calibration_f(*th1, q) - calibration_f(*th2, q)) / lam**2,
                        (th2 - th1) * SIGNS / lam,
                        w.v_tilde(first) - w.v_tilde(second),
                    )
        keys = [(u1, tau) for u1 in u_values for tau in taus]
        check(
            f"d1 q={q}: exact wave causal margin > 0 for all 45 constructed pairs (r=64..1024)",
            all(
                table[(r, *key)][0] > 0 and table[(r, *key)][1] > 0
                for r in rungs
                for key in keys
            ),
            "min margin/lambda^2 after lowering v_tilde: "
            + f"{min(table[(r, *key)][1] for r in rungs for key in keys):.4f}",
        )
        fails = {r: sum(table[(r, *key)][2] <= 0 for key in keys) for r in rungs}
        print(
            f"INFO q={q}: constructed pairs whose calibrated coordinate does NOT decrease "
            "(finite-lambda bias), out of 9: "
            + ", ".join(f"r={r}: {fails[r]}" for r in rungs)
        )
        check(
            f"d2 q={q}: F(theta_Y) < F(theta_X) for all 9 pairs at r=512 and r=1024",
            fails[512] == 0 and fails[1024] == 0,
            "decrease/lambda^2 at r=1024: "
            + ", ".join(f"{table[(1024, *key)][2]:.4f}" for key in keys),
        )
        dev = {
            r: max(abs(table[(r, *key)][2] / table[(r, *key)][1] - 1) for key in keys)
            for r in (512, 1024)
        }
        check(
            f"d3 q={q}: F decrease -> v_tilde decrease (= exact margin after the shift)",
            dev[1024] < 0.10 and dev[1024] < 0.35 * dev[512],
            f"max relative deviation r=512: {dev[512]:.3f}, r=1024: {dev[1024]:.3f}",
        )
        vt_ok = max(
            abs(table[(r, *key)][4] * r * r - table[(r, *key)][1])
            for r in rungs
            for key in keys
        )
        check(
            f"d4 q={q}: after the shift, v_tilde decrease equals the remaining exact margin",
            vt_ok < 1e-6,
            f"max |v_tilde decrease - margin|/lambda^2 = {vt_ok:.2e}",
        )
        # e: bias scaling.  bias(lambda/2)/bias(lambda) = (1/16) * log ratio * (1+o(1)).
        ratios = [
            16 * bias[(2 * r, u1)] / bias[(r, u1)]
            for r in (128, 256, 512)
            for u1 in u_values
        ]
        coeff = {
            r: max(abs(bias[(r, u1)]) * r**4 / math.log(r) for u1 in u_values)
            for r in rungs
        }
        check(
            f"e  q={q}: F(theta_X) - v_tilde(X) = O(lambda^4 log(1/lambda)) at 3 u values",
            0.9 < min(ratios) and max(ratios) < 1.6 and max(coeff.values()) < 1e5,
            f"16*bias(2r)/bias(r) in [{min(ratios):.3f}, {max(ratios):.3f}]; "
            "sup |bias|/(lambda^4 log r): "
            + ", ".join(f"r={r}: {coeff[r]:.0f}" for r in rungs),
        )
        inc = {r: np.array([table[(r, *key)][3] for key in keys]) for r in rungs}
        signs_ok = all(np.all(inc[r] > 0) for r in rungs)
        l1 = {r: float(inc[r].sum(axis=1).max()) for r in rungs}
        check(
            f"f  q={q}: increments have the nestedness signs and |theta_Y-theta_X|_1 = O(lambda)",
            signs_ok and max(l1.values()) < 3.0 and abs(l1[1024] / l1[256] - 1) < 0.05,
            "max |d|_1/lambda: " + ", ".join(f"r={r}: {l1[r]:.3f}" for r in rungs),
        )


def info_constants(w, first, th1, th2, decrease):
    """Size of the constants hidden in (8): local gradient/Hessian of F by differences."""
    q = w.q
    theta = np.array(th1)

    def fun(t):
        return calibration_f(t[0], t[1], t[2], q)

    grad = np.zeros(3)
    hess = np.zeros((3, 3))
    step = 2e-4 * theta
    for i in range(3):
        e = np.zeros(3)
        e[i] = step[i]
        grad[i] = (fun(theta + e) - fun(theta - e)) / (2 * step[i])
        for j in range(3):
            g = np.zeros(3)
            g[j] = step[j]
            hess[i, j] = (
                fun(theta + e + g)
                - fun(theta + e - g)
                - fun(theta - e + g)
                + fun(theta - e - g)
            ) / (4 * step[i] * step[j])
    d = np.abs(np.array(th2) - theta)
    print(
        "INFO local |grad F| at the r=2048 witness: "
        + np.array2string(np.abs(grad), precision=3)
        + "; max |Hessian| entry: "
        + f"{np.abs(hess).max():.3e}"
    )

    def radius(n):
        b = math.ceil(math.log2(32.0) + 4 * math.log2(n))
        k = np.array([n / 2, n / 2, n / 4])
        r_d = np.sqrt(2 * d * b / k) + 2 * b / (3 * k)
        e_b = np.sqrt(2 * theta * b / k) + 2 * b / (3 * k)
        return float(np.abs(grad) @ r_d + e_b @ np.abs(hess) @ (d + r_d))

    lo, hi = 1e6, 1e40
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if radius(mid) < decrease / 2:
            hi = mid
        else:
            lo = mid
    print(
        "INFO event count at which the ideal-constant version of radius (4) drops below "
        f"half the r=2048 signal: N ~ {hi:.2e} (lambda^-3 = {2048.0**3:.2e}, "
        f"ratio {hi / 2048.0**3:.2e}); asymptotic statement only, no finite-N claim."
    )


if __name__ == "__main__":
    part_a()
    part_b()
    wave, x1, x2, t1, t2, dec = part_c()
    part_e0()
    part_def()
    info_constants(wave, x1, t1, t2, dec)
    print(f"SUMMARY s4_population_witness: {len(FAILS)} FAIL")
    sys.exit(1 if FAILS else 0)
