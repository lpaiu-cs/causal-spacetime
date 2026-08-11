"""Empirical Bernstein interval -- Maurer-Pontil (2009) Theorem 4.

The O4 review ruling: coverage is bought by the THEOREM, not by a
simulation. A fixed-seed coverage simulation is a useful regression
diagnostic and is labelled as one; it is never called a certification
and never gates the freeze. What is pinned here is the formula, its
preconditions, and the identity of this implementation with the
theorem's statement.

Maurer & Pontil, "Empirical Bernstein Bounds and Sample-Variance
Penalization" (arXiv:0907.3740), Theorem 4: for i.i.d. Z_1..Z_n in
[0, 1] and delta in (0, 1), with probability at least 1 - delta,

    E[Z] <= Zbar + sqrt(2 V_n ln(2/delta) / n)
                 + 7 ln(2/delta) / (3 (n - 1)),

    V_n = sum_{i<j} (Z_i - Z_j)^2 / (n (n - 1))          (= the
    unbiased sample variance, n-1 denominator).

The bound is ONE-SIDED. Applying it to Z and to 1 - Z -- which share
the same V_n, since the variance is invariant under Z -> 1 - Z --
and taking a union bound gives the two-sided interval

    E[Z] in Zbar +/- half_width,     coverage >= 1 - 2 delta.

O4 uses delta = 0.02 per side, hence a two-sided 96% interval, and
spends the remaining familywise budget on G2.

Preconditions, all fail-closed (a violated precondition voids the
theorem, so degrading to a heuristic would silently void coverage):

  * Z_i in [0, 1] exactly -- values are NOT clipped, because clipping
    biases the mean and voids the i.i.d.-in-[0,1] hypothesis; an
    out-of-range value aborts.
  * n >= 2 (V_n and the 1/(n-1) term are undefined at n = 1).
  * every accumulated value finite.

The i.i.d. unit is the SAMPLE POINT. Streams/partitions are an
implementation detail of generation and must not enter the statistic:
`Accumulator.merge` exists so a partitioned run produces bit-identical
mean and variance to a single-stream run of the same points.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


class EBError(ValueError):
    """A precondition of Theorem 4 failed; no interval is produced."""


@dataclass(frozen=True)
class EBInterval:
    """Two-sided empirical-Bernstein interval on E[Z], in Z units."""

    n: int
    mean: float
    var: float
    half_width: float
    delta: float

    @property
    def lo(self) -> float:
        return self.mean - self.half_width

    @property
    def hi(self) -> float:
        return self.mean + self.half_width

    def rescaled(self, scale: float) -> tuple[float, float]:
        """The interval carried back to the estimand's own units.

        O4 estimates V = scale * E[Z] with scale = dt * B_box, so the
        volume interval is the Z interval times a positive constant --
        a monotone map, which preserves coverage exactly."""

        if not scale > 0.0:
            raise EBError(f"scale must be positive, got {scale}")
        return scale * self.lo, scale * self.hi


class Accumulator:
    """Streaming mean/unbiased-variance (Welford) over sample points.

    Welford rather than sum-of-squares: at O4 sizes (2.6e7 points with
    a mean near 2e-3) the naive form loses most of the variance's
    significant digits to cancellation, and V_n is what the half-width
    is built from."""

    __slots__ = ("n", "_mean", "_m2")

    def __init__(self) -> None:
        self.n = 0
        self._mean = 0.0
        self._m2 = 0.0

    def add(self, z: float) -> None:
        if not math.isfinite(z):
            raise EBError(f"non-finite sample {z}")
        if z < 0.0 or z > 1.0:
            raise EBError(
                f"sample {z!r} outside [0, 1] -- Theorem 4 does not "
                f"apply and clipping would bias the mean")
        self.n += 1
        d = z - self._mean
        self._mean += d / self.n
        self._m2 += d * (z - self._mean)

    def merge(self, other: Accumulator) -> None:
        """Chan-Golub-LeVeque pooling, so that a partitioned run and a
        single-stream run over the same points agree."""

        if other.n == 0:
            return
        if self.n == 0:
            self.n, self._mean, self._m2 = other.n, other._mean, other._m2
            return
        n = self.n + other.n
        d = other._mean - self._mean
        self._m2 += other._m2 + d * d * self.n * other.n / n
        self._mean += d * other.n / n
        self.n = n

    @property
    def mean(self) -> float:
        if self.n == 0:
            raise EBError("no samples")
        return self._mean

    @property
    def var(self) -> float:
        """V_n: the unbiased sample variance (n-1 denominator), which
        is exactly the pairwise form in the theorem."""

        if self.n < 2:
            raise EBError("V_n needs n >= 2")
        return self._m2 / (self.n - 1)


def half_width(n: int, var: float, delta: float) -> float:
    """Theorem 4's deviation term, transcribed literally."""

    if n < 2:
        raise EBError("Theorem 4 needs n >= 2")
    if not 0.0 < delta < 1.0:
        raise EBError(f"delta must lie in (0, 1), got {delta}")
    if var < 0.0 or not math.isfinite(var):
        raise EBError(f"invalid V_n {var}")
    ln = math.log(2.0 / delta)
    return math.sqrt(2.0 * var * ln / n) + 7.0 * ln / (3.0 * (n - 1))


def interval(acc: Accumulator, delta: float) -> EBInterval:
    """Two-sided interval at coverage >= 1 - 2*delta (union bound over
    Theorem 4 applied to Z and to 1 - Z, which share V_n)."""

    return EBInterval(n=acc.n, mean=acc.mean, var=acc.var,
                      half_width=half_width(acc.n, acc.var, delta),
                      delta=delta)


def lower_bound(acc: Accumulator, delta: float) -> float:
    """One-sided lower bound on E[Z] at coverage >= 1 - delta.

    G2 spends its familywise budget as delta/2 on this bound and
    delta/2 on the Clopper-Pearson upper bound, so the two ends of the
    leakage decision carry one shared error probability (review R2:
    the v0.1 draft mixed a 0.01 upper bound with an uncharged lower
    bound)."""

    return acc.mean - half_width(acc.n, acc.var, delta)


def pairwise_variance(values) -> float:
    """V_n by the theorem's literal pairwise formula.

    Quadratic, so this exists to PIN the identity `Accumulator.var ==
    pairwise_variance` on small samples in the contract tests -- never
    to be used at campaign sizes."""

    xs = list(values)
    n = len(xs)
    if n < 2:
        raise EBError("V_n needs n >= 2")
    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = xs[i] - xs[j]
            total += d * d
    return total / (n * (n - 1))
