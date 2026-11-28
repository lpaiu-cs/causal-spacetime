# ruff: noqa: E501, UP031
"""Referee check of the 'elementary check' at the end of Section 3 of
gram_fourier_noise.md.

Claim: X, R uniform on [0,1]; n references, m independent probes;
    S = (1/m) sum_i (F_n(X_i) - X_i) exp(-i omega X_i),  omega = 2 pi j, j != 0.
Then  E|S|^2 = (1 - 1/m)/(n omega^2) + 1/(6 n m).

Derivation checked here:
  diagonal      E (F_n(X)-X)^2 = E X(1-X)/n = 1/(6n)
  off-diagonal  (1/n) int int (min(x,y)-xy) exp(-i omega (x-y)) dx dy
                = (1 - |phi(omega)|^2)/(n omega^2),  phi = E exp(-i omega R),
  and phi(2 pi j) = 0.
Part A does the two integrals exactly with sympy; part B is a Monte Carlo.
Part C checks the general-omega version, to show what the restriction to
omega = 2 pi j buys (only phi = 0).

Deterministic; prints PASS/FAIL; exits nonzero on any FAIL.
"""

from __future__ import annotations

import math
import sys

import numpy as np
import sympy as sp

FAILS = []


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL") + " " + name + (" :: " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def exact_part():
    x, y = sp.symbols("x y", real=True)
    j = sp.symbols("j", integer=True, positive=True)
    w = 2 * sp.pi * j
    diag = sp.integrate(x * (1 - x), (x, 0, 1))
    check("A1 diagonal E[X(1-X)] = 1/6", diag == sp.Rational(1, 6), str(diag))
    # Brownian-bridge covariance min(x,y)-xy, split at x=y; the sine part is odd.
    low = sp.integrate(
        sp.integrate((y - x * y) * sp.cos(w * (x - y)), (y, 0, x)), (x, 0, 1)
    )
    high = sp.integrate(
        sp.integrate((x - x * y) * sp.cos(w * (x - y)), (y, x, 1)), (x, 0, 1)
    )
    total = sp.simplify(low + high)
    check(
        "A2 symbolic off-diagonal integral = 1/omega^2 for integer j",
        sp.simplify(total - 1 / w**2) == 0,
        "integral = " + str(total),
    )
    odd = sp.integrate(
        sp.integrate((y - x * y) * sp.sin(w * (x - y)), (y, 0, x)), (x, 0, 1)
    ) + sp.integrate(
        sp.integrate((x - x * y) * sp.sin(w * (x - y)), (y, x, 1)), (x, 0, 1)
    )
    check(
        "A3 imaginary (sine) part vanishes",
        sp.simplify(odd) == 0,
        str(sp.simplify(odd)),
    )
    # General real omega: (1 - sinc^2(omega/2)) / omega^2.
    om = sp.symbols("omega", positive=True)
    lowg = sp.integrate(
        sp.integrate((y - x * y) * sp.cos(om * (x - y)), (y, 0, x)), (x, 0, 1)
    )
    highg = sp.integrate(
        sp.integrate((x - x * y) * sp.cos(om * (x - y)), (y, x, 1)), (x, 0, 1)
    )
    general = sp.simplify(lowg + highg)
    target = (1 - (sp.sin(om / 2) / (om / 2)) ** 2) / om**2
    gap = max(
        abs(float((general - target).subs(om, val))) for val in (0.7, 3.3, 11.0, 40.5)
    )
    check(
        "A4 general omega: (1-|phi|^2)/omega^2 with phi=sinc",
        gap < 1e-12,
        f"max gap {gap:.2e}",
    )


def formula(n, m, omega):
    phi2 = (math.sin(omega / 2) / (omega / 2)) ** 2
    return (1 - 1 / m) * (1 - phi2) / (n * omega**2) + 1 / (6 * n * m)


def simulate(rng, n, m, omega, reps, chunk=20000):
    values = []
    for start in range(0, reps, chunk):
        size = min(chunk, reps - start)
        refs = rng.random((size, n))
        probes = rng.random((size, m))
        ecdf = (refs[:, None, :] <= probes[:, :, None]).mean(axis=2)
        score = ((ecdf - probes) * np.exp(-1j * omega * probes)).mean(axis=1)
        values.append(np.abs(score) ** 2)
    values = np.concatenate(values)
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(values.size))


def monte_carlo_part():
    rng = np.random.default_rng(20261124)
    cases = [
        (5, 1, 1, 400000),
        (20, 7, 1, 300000),
        (20, 7, 2, 300000),
        (50, 50, 3, 120000),
        (64, 200, 5, 60000),
    ]
    for n, m, j, reps in cases:
        omega = 2 * math.pi * j
        mean, se = simulate(rng, n, m, omega, reps)
        stated = (1 - 1 / m) / (n * omega**2) + 1 / (6 * n * m)
        z = (mean - stated) / se
        check(
            "B Monte Carlo n=%d m=%d j=%d" % (n, m, j),
            abs(z) < 4.0,
            f"stated {stated:.6e}  simulated {mean:.6e} +- {se:.1e}  z={z:+.2f}",
        )
    # Non-vacuity: with m>1 the two terms are distinguishable from either alone.
    n, m, j = 50, 50, 3
    omega = 2 * math.pi * j
    mean, se = simulate(np.random.default_rng(7), n, m, omega, 120000)
    only_first = (1 - 1 / m) / (n * omega**2)
    only_second = 1 / (6 * n * m)
    check(
        "B' both terms are needed (each alone is rejected at >10 sigma)",
        abs(mean - only_first) / se > 10 and abs(mean - only_second) / se > 10,
        f"z(first only)={(mean - only_first) / se:.1f} z(second only)={(mean - only_second) / se:.1f}",
    )
    # General omega (not a multiple of 2 pi): the stated formula is NOT valid,
    # the (1-|phi|^2) version is.
    n, m, omega = 20, 7, 4.0
    mean, se = simulate(np.random.default_rng(11), n, m, omega, 300000)
    z_general = (mean - formula(n, m, omega)) / se
    z_stated = (mean - ((1 - 1 / m) / (n * omega**2) + 1 / (6 * n * m))) / se
    check(
        "C general omega=4.0 follows the (1-|phi|^2) version",
        abs(z_general) < 4.0 and abs(z_stated) > 10,
        f"z(general)={z_general:+.2f}  z(2 pi j formula)={z_stated:+.1f}",
    )


def independent_reference_part():
    """Section 2 says independent noise per probe 'would be the wrong error
    model'.  Check in which direction: give every probe its OWN n references."""
    rng = np.random.default_rng(31)
    n, m, j, reps = 20, 7, 1, 300000
    omega = 2 * math.pi * j
    refs = rng.random((reps, m, n))
    probes = rng.random((reps, m))
    ecdf = (refs <= probes[:, :, None]).mean(axis=2)
    values = np.abs(((ecdf - probes) * np.exp(-1j * omega * probes)).mean(axis=1)) ** 2
    mean, se = float(values.mean()), float(values.std(ddof=1) / math.sqrt(reps))
    indep = 1 / (6 * n * m)
    shared = (1 - 1 / m) / (n * omega**2) + indep
    check(
        "D independent references per probe give 1/(6nm) only: SMALLER than the shared-reference variance",
        abs(mean - indep) < 4 * se and shared > indep and abs(mean - shared) > 10 * se,
        f"independent {indep:.6e} (MC {mean:.6e} +- {se:.1e}), shared {shared:.6e}; the common-mode term (1-1/m)/(n omega^2) is what (3) controls",
    )


if __name__ == "__main__":
    exact_part()
    monte_carlo_part()
    independent_reference_part()
    print("SUMMARY failures=%d" % len(FAILS))
    sys.exit(1 if FAILS else 0)
