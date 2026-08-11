"""Directed-rounding interval arithmetic on MPFR -- the floating-point
backend of the volume-oracle certification contract
(schwarzschild_volume_oracle_certification.md section 0).

Every `Iv` is a closed interval [lo, hi] of MPFR numbers: lower
endpoints are computed under RoundDown, upper endpoints under RoundUp,
at 96-bit precision. MPFR documents correct rounding in the requested
direction for every operation used here, including the transcendentals
(log, sin, cos, const_pi) that platform libm does not guarantee --
which is why this module exists instead of `math` + `nextafter`.

Fail-closed: any operation whose precondition cannot be certified
(division by an interval containing zero, sqrt/log reaching a
non-positive endpoint, trig outside its certified monotone domain, an
empty intersection) raises CertificationError. Nothing here degrades
to a heuristic.
"""

from __future__ import annotations

import gmpy2

PRECISION = 96
_DN = gmpy2.context(precision=PRECISION, round=gmpy2.RoundDown)
_UP = gmpy2.context(precision=PRECISION, round=gmpy2.RoundUp)


class CertificationError(ValueError):
    """A certified precondition failed; no enclosure is produced."""


def _mpfr_exact(x) -> gmpy2.mpfr:
    """Embed a Python number into MPFR, refusing anything inexact.

    binary64 floats and small ints embed exactly at 96 bits; mpfr
    values pass through. Anything else (strings, huge ints, mpq)
    must be converted by the caller through an explicit interval
    operation so the rounding direction is visible."""

    if isinstance(x, gmpy2.mpfr):
        v = x
    elif isinstance(x, bool):
        raise CertificationError("bool is not a certified scalar")
    elif isinstance(x, int):
        v = gmpy2.mpfr(x, PRECISION)
        if v != x:
            raise CertificationError(f"int {x} does not embed exactly")
    elif isinstance(x, float):
        v = gmpy2.mpfr(x, PRECISION)
    else:
        raise CertificationError(f"cannot embed {type(x).__name__}")
    if not gmpy2.is_finite(v):
        raise CertificationError(f"non-finite endpoint {v}")
    return v


class Iv:
    """Closed MPFR interval [lo, hi], lo <= hi, both finite."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = _mpfr_exact(lo)
        self.hi = self.lo if hi is None else _mpfr_exact(hi)
        if not self.lo <= self.hi:
            raise CertificationError(f"inverted interval [{lo}, {hi}]")

    # -- construction helpers -----------------------------------------

    @staticmethod
    def pi() -> Iv:
        return _wrap(_DN.const_pi(), _UP.const_pi())

    def __repr__(self) -> str:
        return f"Iv({self.lo!r}, {self.hi!r})"

    # -- arithmetic (directed on both endpoints) ----------------------

    def __add__(self, other) -> Iv:
        o = _as_iv(other)
        return _wrap(_DN.add(self.lo, o.lo), _UP.add(self.hi, o.hi))

    __radd__ = __add__

    def __neg__(self) -> Iv:
        return _wrap(-self.hi, -self.lo)

    def __sub__(self, other) -> Iv:
        o = _as_iv(other)
        return _wrap(_DN.sub(self.lo, o.hi), _UP.sub(self.hi, o.lo))

    def __rsub__(self, other) -> Iv:
        return _as_iv(other) - self

    def __mul__(self, other) -> Iv:
        o = _as_iv(other)
        pairs = ((self.lo, o.lo), (self.lo, o.hi),
                 (self.hi, o.lo), (self.hi, o.hi))
        return _wrap(min(_DN.mul(a, b) for a, b in pairs),
                     max(_UP.mul(a, b) for a, b in pairs))

    __rmul__ = __mul__

    def __truediv__(self, other) -> Iv:
        o = _as_iv(other)
        if o.lo <= 0 <= o.hi:
            raise CertificationError(f"division by {o} containing 0")
        pairs = ((self.lo, o.lo), (self.lo, o.hi),
                 (self.hi, o.lo), (self.hi, o.hi))
        return _wrap(min(_DN.div(a, b) for a, b in pairs),
                     max(_UP.div(a, b) for a, b in pairs))

    def __rtruediv__(self, other) -> Iv:
        return _as_iv(other) / self

    def sq(self) -> Iv:
        """x^2 -- tighter than self * self when the interval spans 0."""

        if self.lo >= 0:
            return _wrap(_DN.mul(self.lo, self.lo),
                         _UP.mul(self.hi, self.hi))
        if self.hi <= 0:
            return _wrap(_DN.mul(self.hi, self.hi),
                         _UP.mul(self.lo, self.lo))
        m = max(-self.lo, self.hi)
        return _wrap(gmpy2.mpfr(0), _UP.mul(m, m))

    def sqrt(self) -> Iv:
        if not self.lo > 0:
            raise CertificationError(f"sqrt of {self} not certified > 0")
        return _wrap(_DN.sqrt(self.lo), _UP.sqrt(self.hi))

    def sqrt0(self) -> Iv:
        """sqrt clamped at 0: certified for intervals whose true value
        is known nonnegative but whose lower endpoint may round below
        zero (e.g. u_t - u at the turning point)."""

        if self.hi < 0:
            raise CertificationError(f"sqrt0 of negative {self}")
        lo = _DN.sqrt(self.lo) if self.lo > 0 else gmpy2.mpfr(0)
        return _wrap(lo, _UP.sqrt(self.hi))

    def log(self) -> Iv:
        if not self.lo > 0:
            raise CertificationError(f"log of {self} not certified > 0")
        return _wrap(_DN.log(self.lo), _UP.log(self.hi))

    def cos(self) -> Iv:
        """cos on the certified monotone-decreasing domain [0, pi]."""

        if not (self.lo >= 0 and self.hi <= _DN.const_pi()):
            raise CertificationError(
                f"cos of {self} outside certified domain [0, pi]")
        return _wrap(_DN.cos(self.hi), _UP.cos(self.lo))

    def sin(self) -> Iv:
        """sin on the certified monotone-increasing domain
        [0, pi/2]."""

        half_pi = _DN.div(_DN.const_pi(), gmpy2.mpfr(2))
        if not (self.lo >= 0 and self.hi <= half_pi):
            raise CertificationError(
                f"sin of {self} outside certified domain [0, pi/2]")
        return _wrap(_DN.sin(self.lo), _UP.sin(self.hi))

    def abs(self) -> Iv:
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return _wrap(gmpy2.mpfr(0), max(-self.lo, self.hi))

    # -- certified comparisons (three-valued) -------------------------

    def certainly_lt(self, other) -> bool:
        return self.hi < _as_iv(other).lo

    def certainly_gt(self, other) -> bool:
        return self.lo > _as_iv(other).hi

    def contains(self, x) -> bool:
        """Point membership; exact for int/float/mpfr/mpq operands
        (gmpy2 compares mixed types exactly)."""

        return self.lo <= x <= self.hi

    # -- lattice operations -------------------------------------------

    def hull(self, other) -> Iv:
        o = _as_iv(other)
        return _wrap(min(self.lo, o.lo), max(self.hi, o.hi))

    def intersect(self, other) -> Iv:
        """Intersection; EMPTY is a certification failure (two
        enclosures of the same true value cannot be disjoint)."""

        o = _as_iv(other)
        lo, hi = max(self.lo, o.lo), min(self.hi, o.hi)
        if lo > hi:
            raise CertificationError(f"empty intersection {self} ^ {o}")
        return _wrap(lo, hi)

    def widen(self, delta: Iv) -> Iv:
        """[lo - max(delta), hi + max(delta)] for a certified
        nonnegative error radius delta."""

        d = _as_iv(delta)
        if d.lo < 0:
            raise CertificationError(f"widen by signed {d}")
        return _wrap(_DN.sub(self.lo, d.hi), _UP.add(self.hi, d.hi))

    # -- directed float conversions (outward: coverage-preserving) ----

    def lo_float(self) -> float:
        """self.lo rounded DOWN to binary64 -- never above the true
        lower endpoint. Plain float() rounds to nearest and can shrink
        an enclosure inward; certified paths must use these."""

        return float_down(self.lo)

    def hi_float(self) -> float:
        """self.hi rounded UP to binary64 -- never below the true
        upper endpoint."""

        return float_up(self.hi)

    # -- non-certified conveniences (diagnostics only) ----------------

    def width(self) -> float:
        return float(_UP.sub(self.hi, self.lo))

    def mid(self) -> float:
        return float(self.lo + self.hi) / 2.0


def _wrap(lo: gmpy2.mpfr, hi: gmpy2.mpfr) -> Iv:
    iv = Iv.__new__(Iv)
    iv.lo, iv.hi = lo, hi
    if not (gmpy2.is_finite(lo) and gmpy2.is_finite(hi) and lo <= hi):
        raise CertificationError(f"bad enclosure [{lo}, {hi}]")
    return iv


def _as_iv(x) -> Iv:
    return x if isinstance(x, Iv) else Iv(x)


def iv_min(a: Iv, b: Iv) -> Iv:
    return _wrap(min(a.lo, b.lo), min(a.hi, b.hi))


def iv_max(a: Iv, b: Iv) -> Iv:
    return _wrap(max(a.lo, b.lo), max(a.hi, b.hi))


_DN53 = gmpy2.context(precision=53, round=gmpy2.RoundDown)
_UP53 = gmpy2.context(precision=53, round=gmpy2.RoundUp)
_ZERO = gmpy2.mpfr(0)


def float_down(x: gmpy2.mpfr) -> float:
    """x rounded toward -inf to binary64 (result <= x, exactly
    representable): a directed 53-bit rounding, then an exact
    float() of the 53-bit value."""

    return float(_DN53.add(x, _ZERO))


def float_up(x: gmpy2.mpfr) -> float:
    """x rounded toward +inf to binary64 (result >= x)."""

    return float(_UP53.add(x, _ZERO))


def iv_sum(terms) -> Iv:
    """Directed-rounding sum of an iterable of intervals -- the L6d
    composition primitive."""

    lo, hi = gmpy2.mpfr(0), gmpy2.mpfr(0)
    for t in terms:
        lo = _DN.add(lo, t.lo)
        hi = _UP.add(hi, t.hi)
    return _wrap(lo, hi)
