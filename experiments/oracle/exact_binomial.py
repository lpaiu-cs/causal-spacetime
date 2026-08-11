"""Exact (Clopper-Pearson) binomial bounds, without scipy.

G2's upper end and G3's cluster bound are both "how large can the
underlying rate be, given this many successes out of this many
independent trials" -- the Clopper-Pearson construction, which is
exact (conservative) for every n and k rather than asymptotic.

    p_upper(k, n, alpha) = the p solving P[Bin(n, p) <= k] = alpha

with p_upper = 1 for k = n. The regularized incomplete beta is
implemented here (Lentz continued fraction) because this host has no
scipy and the freeze may not acquire a new dependency: the O3
environment lock names python/gmpy2/MPFR/GMP, and O4 inherits it.

The zero-count case has the closed form 1 - alpha^(1/n), which the
contract test uses as an independent check of the general path.
"""

from __future__ import annotations

import math

_ITMAX, _EPS, _TINY = 300, 1e-15, 1e-300


def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < _TINY:
        d = _TINY
    d = 1.0 / d
    h = d
    for m in range(1, _ITMAX + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < _TINY:
            d = _TINY
        c = 1.0 + aa / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < _TINY:
            d = _TINY
        c = 1.0 + aa / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < _EPS:
            break
    return h


def betai(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""

    if not 0.0 <= x <= 1.0:
        raise ValueError(f"x must lie in [0, 1], got {x}")
    if x in (0.0, 1.0):
        return x
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def binom_cdf(k: int, n: int, p: float) -> float:
    """P[Bin(n, p) <= k] = I_{1-p}(n - k, k + 1)."""

    if k >= n:
        return 1.0
    if k < 0:
        return 0.0
    return betai(n - k, k + 1.0, 1.0 - p)


def cp_upper(k: int, n: int, alpha: float) -> float:
    """One-sided Clopper-Pearson upper bound on the success rate.

    Solves P[Bin(n, p) <= k] = alpha by bisection; the CDF is
    decreasing in p, so the root is unique."""

    if n <= 0:
        raise ValueError("n must be positive")
    if not 0 <= k <= n:
        raise ValueError(f"k={k} outside [0, {n}]")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must lie in (0, 1), got {alpha}")
    if k == n:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if binom_cdf(k, n, mid) > alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cp_upper_zero(n: int, alpha: float) -> float:
    """Closed form for k = 0: 1 - alpha^(1/n)."""

    if n <= 0:
        raise ValueError("n must be positive")
    return 1.0 - alpha ** (1.0 / n)
