"""An order-only null-coordinate consistency test and a separate wave fixture.

The observer consumes a transitive order and known 0<q<=1. Generator coordinates
are used only in the deterministic --output challenge, never by order_test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from order_wave_energy import (
    background_probability,
    relation_probability,
    rosen_kernel,
    rule,
)
from order_wave_resolution import (
    circular_correction,
    circular_primitive,
    convolution,
    pair_change,
)
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/null_coordinate_challenge.md"


@dataclass(frozen=True)
class Box:
    """Rational intervals; no floating rounding is used by the observer."""

    lo: F
    hi: F

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(value):
        return value if isinstance(value, Box) else Box(F(value), F(value))

    def __add__(self, other):
        other = self.point(other)
        return Box(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.point(other)

    def __rsub__(self, other):
        return self.point(other) + -self

    def __mul__(self, other):
        other = self.point(other)
        values = [x * y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Box(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.point(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval includes zero")
        return self * Box(1 / other.hi, 1 / other.lo)

    def sqrt(self, bits):
        if self.lo < 0:
            raise ValueError("negative square root")
        denominator = 1 << bits

        def floor_root(x):
            return math.isqrt(x.numerator * denominator**2 // x.denominator)

        return Box(
            F(floor_root(self.lo), denominator), F(floor_root(self.hi) + 1, denominator)
        )


def series(z, r, bits):
    """Enclose C_r(z); r=0 instead means cos(sqrt(z)), on 0<=z<=1."""
    z = F(z)
    if not 0 <= z <= 1:
        raise ValueError("series certificate is restricted to [0,1]")
    term = total = F(1, math.factorial(2 * r - 1)) if r else F(1)
    n = 0
    while True:
        divisor = 2 * (n + 1) * (2 * r + 2 * n + 1) if r else (2 * n + 2) * (2 * n + 1)
        following = -term * z / divisor
        if abs(following) <= F(1, 1 << bits):
            return Box(min(total, total + following), max(total, total + following))
        total += following
        term = following
        n += 1


def p0_box(z, bits):
    c1, c2, c3 = (series(z, r, bits) for r in (1, 2, 3))
    return c1 * c3 / (3 * c2 * c2)


def monotonicity_certificate():
    bounds = []
    for j in range(16):
        a, b = F(j, 16), F(j + 1, 16)
        ca = [None] + [series(a, r, 100) for r in range(1, 5)]
        cb = [None] + [series(b, r, 100) for r in range(1, 5)]
        upper = -cb[2].lo / ca[1].hi - 3 * cb[4].lo / ca[3].hi + 4 * ca[3].hi / cb[2].lo
        assert upper < -F(1, 50)
        bounds.append(float(upper))
    assert p0_box(F(1), 100).lo > F(9, 100)
    return max(bounds)


def inverse_z(target, q, bits):
    """Enclose the clipped inverse p0(z)=target, with the certified slope bound."""
    if target >= F(1, 10):
        return Box.point(0)
    if target <= p0_box(q, 2 * bits + 16).lo:
        return Box.point(q)
    lo, hi = F(0), q
    for _ in range(bits + 1):
        mid = (lo + hi) / 2
        value = p0_box(mid, 2 * bits + 16)
        if value.lo > target:
            lo = mid
        elif value.hi < target:
            hi = mid
        else:
            radius = 1024 * (value.hi - value.lo)
            return Box(max(lo, mid - radius), min(hi, mid + radius))
    return Box(lo, hi)


def alpha_box(s, q, bits):
    return s * s * series(q * s * s, 2, bits) / series(q * s * s, 1, bits)


def t_box(s, q, bits):
    u = s - F(1, 2)
    half_tangent = F(1, 2) * series(q / 4, 1, bits) / series(q / 4, 0, bits)
    at_u = u * series(q * u * u, 1, bits) / series(q * u * u, 0, bits)
    return (at_u + half_tangent) / (2 * half_tangent)


def null_v_box(p, f, chain, q, bits=64):
    """Background-compatible v interval using only three probability intervals."""
    q = F(q)
    if not 0 < q <= 1:
        raise ValueError("observer requires known 0<q<=1")
    if chain.lo > F(1, 10) or chain.hi < p0_box(q, 2 * bits + 16).lo:
        return None
    lower = (inverse_z(chain.hi, q, bits) / q).sqrt(bits)
    upper = (inverse_z(chain.lo, q, bits) / q).sqrt(bits)
    s = Box(lower.lo, min(F(1), upper.hi))
    if s.lo <= 0 or s.hi >= 1:
        return None
    precision = 2 * bits + 16
    a = Box(alpha_box(s.lo, q, precision).lo, alpha_box(s.hi, q, precision).hi)
    b = Box(alpha_box(1 - s.hi, q, precision).lo, alpha_box(1 - s.lo, q, precision).hi)
    gamma = alpha_box(F(1), q, precision)
    eta, zeta = (p * gamma / a).sqrt(bits), (f * gamma / b).sqrt(bits)
    t = Box(t_box(s.lo, q, precision).lo, t_box(s.hi, q, precision).hi)
    return 1 - t + t * eta - (1 - t) * zeta


def probability_box(count, trials, exponent, bits):
    radius = Box.point(F(exponent, 2 * trials)).sqrt(bits).hi
    mean = F(count, trials)
    return Box(max(F(0), mean - radius), min(F(1), mean + radius))


def order_test(order, q, rng):
    """Randomized test at level N^-2; order must contain transitive relations."""
    order = np.asarray(order)
    if order.dtype != bool or order.ndim != 2 or order.shape[0] != order.shape[1]:
        raise ValueError("square boolean transitive-order matrix required")
    n = len(order)
    if n < 4 or np.any(np.diag(order)) or np.any(order & order.T):
        raise ValueError("at least four events in a strict order required")
    q = F(q)
    if not 0 < q <= 1:
        raise ValueError("observer requires known 0<q<=1")
    shuffled = rng.permutation(n)
    references, probes = shuffled[: n // 2], shuffled[n // 2 :]
    exponent, bits = (6 * n**3).bit_length(), 40 + n.bit_length()
    boxes = {}
    # ponytail: direct O(N^3) chain counting; batch/bitset it before a powered run.
    for i in probes:
        ancestors = references[order[references, i]]
        future_count = int(np.count_nonzero(order[i, references]))
        m = len(ancestors)
        if m < 2:
            continue
        pairs = int(np.count_nonzero(order[np.ix_(ancestors, ancestors)]))
        p = probability_box(m, len(references), exponent, bits)
        f = probability_box(future_count, len(references), exponent, bits)
        radius = Box.point(F(exponent, 2 * (m // 2))).sqrt(bits).hi
        mean = F(pairs, math.comb(m, 2))
        chain = Box(max(F(0), mean - radius), min(F(1), mean + radius))
        boxes[int(i)] = null_v_box(p, f, chain, q, bits)
    for i, j in zip(probes[::2], probes[1::2], strict=False):
        first, second = boxes.get(int(i)), boxes.get(int(j))
        if first is None or second is None:
            continue
        if (order[i, j] and second.hi < first.lo) or (
            order[j, i] and first.hi < second.lo
        ):
            return True
    return False


def population_v(p, f, chain, q):
    """Floating population diagnostic, separate from the certified count observer."""
    s = brentq(lambda at: background_probability(q * at**2) - chain, 0, 1, xtol=1e-14)

    def alpha(length):
        return length**2 * convolution(q * length**2, 2) / convolution(q * length**2, 1)

    eta, zeta = np.sqrt(p * alpha(1) / alpha(s)), np.sqrt(f * alpha(1) / alpha(1 - s))
    k = math.sqrt(q)
    t = (math.tan(k * (s - 0.5)) + math.tan(k / 2)) / (2 * math.tan(k / 2))
    return (1 - t) + t * eta - (1 - t) * zeta, s - 0.5


class WaveFixture:
    """Generator-side exact geometry; never passed to order_test or null_v_box."""

    def __init__(self, q, r, *, tighter=False):
        self.q, self.r, self.k = q, r, math.sqrt(q)
        self.sol = circular_primitive(q, 1 / r, tighter=tighter)
        self.M = 2 * math.tan(self.k / 2) / self.k
        self.factor = 2 * math.sqrt(q) / r / math.sqrt(1 + 4 * q / r**2)

    def parts(self, s, t):
        s, t = np.broadcast_arrays(s, t)
        base = (t - s) * np.sinc(self.k * (t - s) / np.pi)
        m = base / (np.cos(self.k * s) * np.cos(self.k * t))
        integral = (
            (self.sol.sol(t.ravel())[0] - self.sol.sol(s.ravel())[0]) * np.exp(0.37j)
        ).reshape(s.shape)
        return base, m, integral

    def matrix(self, s, t):
        _, m, z = self.parts(s, t)
        j = np.stack((z.real, z.imag, z.imag, -z.real), axis=-1).reshape(
            z.shape + (2, 2)
        )
        return m[..., None, None] * np.eye(2) - self.factor * j

    def correction(self, s, t):
        base, m, integral = self.parts(s, t)
        ratio = np.divide(abs(integral), m, out=np.ones_like(m), where=m > 0)
        z = self.factor**2 * ratio**2
        return base, -z / (1 + np.sqrt(1 - z))

    def canonical_point(self, u, z, v):
        left, right = self.matrix(-0.5, u), self.matrix(u, 0.5)
        m = float(self.parts(-0.5, u)[1])
        values, vectors = np.linalg.eigh(
            m * (self.M - m) / self.M * (np.linalg.inv(left) + np.linalg.inv(right))
        )
        y = (vectors * values**-0.5) @ vectors.T @ z
        eta = v - z @ z / (2 * m)
        return np.r_[u, eta + 0.5 * y @ np.linalg.solve(left, y), y]

    def F_matrix(self, u):
        m = float(self.parts(-0.5, u)[1])
        e = self.matrix(-0.5, u) - m * np.eye(2)
        end = self.matrix(-0.5, 0.5) - self.M * np.eye(2)
        return (e - m / self.M * end) / (m * (self.M - m))

    def margin(self, first, second):
        dy = second[2:] - first[2:]
        return float(
            second[1]
            - first[1]
            - 0.5 * dy @ np.linalg.solve(self.matrix(first[0], second[0]), dy)
        )

    def segment(self, start, end, factor=1, *, pairs=False):
        length = end - start
        n = factor * max(96, math.ceil(self.r * length / 2) + 64)

        def normalized(s, t):
            base, relative = self.correction(
                start + length * (s + 0.5), start + length * (t + 0.5)
            )
            return base / length, relative

        effective_q = self.q * length**2
        if pairs:
            change = pair_change(normalized, effective_q, n)
            dc2 = change["c2_change"]
            probability = background_probability(effective_q) + change["p2_change"]
        else:
            x, w = rule(n)
            left, dl = normalized(-0.5, x - 0.5)
            right, dr = normalized(x - 0.5, 0.5)
            dc2 = float(np.sum(w * left * right * (dl + dr + dl * dr)))
            probability = None
        base, relative = normalized(-0.5, 0.5)
        alpha = (
            length**2 * (convolution(effective_q, 2) + dc2) / (base * (1 + relative))
        )
        return float(alpha), probability

    def moments(self, point, factor=1):
        u, v, x, y = point
        transverse = np.array([x, y])
        eta = v - 0.5 * transverse @ np.linalg.solve(self.matrix(-0.5, u), transverse)
        zeta = (
            1 - v - 0.5 * transverse @ np.linalg.solve(self.matrix(u, 0.5), transverse)
        )
        gamma = self.segment(-0.5, 0.5, factor)[0]
        left, chain = self.segment(-0.5, u, factor, pairs=True)
        right = self.segment(u, 0.5, factor)[0]
        if min(eta, zeta) <= 0:
            raise RuntimeError("fixture left the diamond")
        return float(left * eta**2 / gamma), float(right * zeta**2 / gamma), chain


def choose_pair(geometry):
    best = None
    offsets = (np.arange(16) - 7.5) * 2 * np.pi / (16 * geometry.r)
    for u in offsets:
        for v in 0.25 + offsets:
            _, vectors = np.linalg.eigh(geometry.F_matrix(u) - geometry.F_matrix(v))
            z = 0.2 * vectors[:, -1]
            first, second = (
                geometry.canonical_point(u, z, 0.5),
                geometry.canonical_point(v, z, 0.5),
            )
            margin = geometry.margin(first, second)
            if best is None or margin > best[0]:
                best = (margin, first, second, z)
    margin, first, second, z = best
    if margin <= 0:
        raise RuntimeError("no positive exact margin in the frozen construction")
    second[1] -= margin / 2
    return first, second, z, margin / 2


def self_check():
    bound = monotonicity_certificate()
    q, u, v, radius = 1.0, 0.1, 0.4, 0.2
    k, s = math.sqrt(q), u + 0.5
    m = (math.tan(k * u) + math.tan(k / 2)) / k
    total = 2 * math.tan(k / 2) / k

    def alpha(length):
        return length**2 * convolution(q * length**2, 2) / convolution(q * length**2, 1)

    eta, zeta = v - radius**2 / (2 * m), 1 - v - radius**2 / (2 * (total - m))
    moments = (
        alpha(s) * eta**2 / alpha(1),
        alpha(1 - s) * zeta**2 / alpha(1),
        background_probability(q * s**2),
    )
    recovered, recovered_u = population_v(*moments, q)
    assert abs(recovered - v) < 1e-12 and abs(recovered_u - u) < 1e-12
    boxes = [Box(F(value) - F(1, 10**13), F(value) + F(1, 10**13)) for value in moments]
    interval = null_v_box(*boxes, q)
    assert interval.lo < F(v) < interval.hi
    assert float(interval.hi - interval.lo) < 1e-8
    assert not order_test(np.zeros((8, 8), dtype=bool), q, np.random.default_rng(1))
    fixture = WaveFixture(1, 32)
    independent, _ = circular_correction(1, 1 / 32, tighter=True)
    base, correction = independent(-0.31, 0.27)
    from_matrix = (
        np.cos(-0.31)
        * np.cos(0.27)
        * np.sqrt(np.linalg.det(fixture.matrix(-0.31, 0.27)))
    )
    assert abs(from_matrix - base * (1 + correction)) < 2e-12
    # The scalar correction above shares a primitive implementation. Also use
    # the original five-state Rosen ODE for an independent formulation check.
    rosen, _ = rosen_kernel(1, 1 / 32, "circular", phase=0.37, tighter=True)
    assert abs(from_matrix - rosen(-0.31, 0.27) / math.sqrt(1 + 4 / 32**2)) < 2e-12
    print(
        json.dumps(
            {
                "controls": "passed",
                "certified_log_derivative_upper": bound,
                "background_v_box_width": float(interval.hi - interval.lo),
            }
        )
    )


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("commit source and design, then run from a clean tree")
    stamp = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )
    output.mkdir(parents=True, exist_ok=False)
    files = [
        Path(__file__),
        Path(__file__).with_name("order_wave_resolution.py"),
        Path(__file__).with_name("order_wave_energy.py"),
        DESIGN,
    ]
    manifest = {
        "status": "in_progress",
        "code_version": stamp,
        "utc_started": time.strftime("%FT%TZ", time.gmtime()),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "files": {
            str(path.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in files
        },
        "scope": "constructed population witnesses, not a finite-N power experiment",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    started = time.time()
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:
        for q in (0.25, 1.0):
            for r in (64, 128, 256, 512, 1024, 2048):
                geometry = WaveFixture(q, r)
                first, second, z, decrease = choose_pair(geometry)
                moments = [geometry.moments(point) for point in (first, second)]
                recovered = [population_v(*row, q) for row in moments]
                row = {
                    "q": q,
                    "r": r,
                    "phase": 0.37,
                    "points": [first.tolist(), second.tolist()],
                    "canonical_z": z.tolist(),
                    "canonical_decrease": decrease,
                    "wave_margin": geometry.margin(first, second),
                    "moments": moments,
                    "recovered": recovered,
                    "population_v_decrease": recovered[0][0] - recovered[1][0],
                    "ode_evaluations": geometry.sol.nfev,
                }
                row["scaled_population_decrease"] = row["population_v_decrease"] * r**2
                row["scaled_wave_margin"] = row["wave_margin"] * r**2
                if q == 1 and r == 2048:
                    refined = WaveFixture(q, r, tighter=True)
                    # Same generator coordinates and observed pair, without reselection.
                    refined_moments = [
                        refined.moments(point, 2) for point in (first, second)
                    ]
                    values = [population_v(*item, q) for item in refined_moments]
                    row["refinement"] = {
                        "moments": refined_moments,
                        "population_v_decrease": values[0][0] - values[1][0],
                        "wave_margin": refined.margin(first, second),
                    }
                stream.write(json.dumps(row) + "\n")
                stream.flush()
                print(
                    json.dumps(
                        {
                            key: value
                            for key, value in row.items()
                            if key not in ("points", "moments", "recovered")
                        }
                    ),
                    flush=True,
                )
    manifest.update(
        status="completed",
        elapsed_seconds=time.time() - started,
        rows_sha256=hashlib.sha256((output / "rows.jsonl").read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


def summarize(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    raw = (folder / "rows.jsonl").read_bytes()
    assert manifest["status"] == "completed"
    assert hashlib.sha256(raw).hexdigest() == manifest["rows_sha256"]
    stamp = manifest["code_version"]
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=ROOT, check=True
    )
    for path, digest in manifest["files"].items():
        saved = subprocess.check_output(["git", "show", f"{stamp}:{path}"], cwd=ROOT)
        assert hashlib.sha256(saved).hexdigest() == digest
    rows = [json.loads(line) for line in raw.splitlines()]
    assert [(row["q"], row["r"]) for row in rows] == [
        (q, r) for q in (0.25, 1.0) for r in (64, 128, 256, 512, 1024, 2048)
    ]
    assert all(row["wave_margin"] > 0 for row in rows)
    summary = {
        "source_commit": stamp,
        "raw_sha256": manifest["rows_sha256"],
        "elapsed_seconds": manifest["elapsed_seconds"],
        "positive_population_witnesses": sum(
            row["population_v_decrease"] > 0 for row in rows
        ),
        "rows": [
            {
                key: row[key]
                for key in (
                    "q",
                    "r",
                    "population_v_decrease",
                    "scaled_population_decrease",
                    "scaled_wave_margin",
                )
            }
            for row in rows
        ],
        "refinement": [
            {
                "q": row["q"],
                "r": row["r"],
                "population_decrease_change": row["refinement"]["population_v_decrease"]
                - row["population_v_decrease"],
                "largest_moment_change": float(
                    np.max(
                        abs(
                            np.array(row["refinement"]["moments"])
                            - np.array(row["moments"])
                        )
                    )
                ),
            }
            for row in rows
            if "refinement" in row
        ],
    }
    (folder / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def cross_check(folder):
    """Post-run independent-formulation check of the largest q=1 witness."""
    rows = [
        json.loads(line) for line in (folder / "rows.jsonl").read_text().splitlines()
    ]
    row = next(item for item in rows if item["q"] == 1 and item["r"] == 2048)
    q, r, phase = row["q"], row["r"], row["phase"]
    kernel, details = rosen_kernel(q, 1 / r, "circular", phase=phase, tighter=True)
    k = math.sqrt(q)

    def order_for(length):
        return 2 * max(96, math.ceil(r * length / 2) + 64)

    def alpha(start, end):
        x, w = rule(order_for(end - start))
        u = start + (end - start) * x
        return float(
            (end - start)
            * np.sum(w * kernel(start, u) * kernel(u, end))
            / kernel(start, end)
        )

    def matrix(start, end):
        x, w = rule(order_for(end - start))
        u = start + (end - start) * x
        weight = (end - start) * w / np.cos(k * u) ** 2
        factor = 2 * k / r / math.sqrt(1 + 4 * q / r**2)
        real, imag = weight @ np.cos(r * u + phase), weight @ np.sin(r * u + phase)
        return weight.sum() * np.eye(2) - factor * np.array(
            [[real, imag], [imag, -real]]
        )

    gamma = alpha(-0.5, 0.5)
    moments = []
    for point in row["points"]:
        u, v, *xy = point
        xy = np.array(xy)
        eta = v - 0.5 * xy @ np.linalg.solve(matrix(-0.5, u), xy)
        zeta = 1 - v - 0.5 * xy @ np.linalg.solve(matrix(u, 0.5), xy)
        length = u + 0.5

        def local(s, t, length=length):
            return kernel(-0.5 + length * (s + 0.5), -0.5 + length * (t + 0.5))

        probability = relation_probability(local, order_for(length))
        moments.append(
            (
                alpha(-0.5, u) * eta**2 / gamma,
                alpha(u, 0.5) * zeta**2 / gamma,
                probability,
            )
        )
    coordinates = [population_v(*item, q)[0] for item in moments]
    gap = coordinates[0] - coordinates[1]
    difference = gap - row["population_v_decrease"]
    result = {
        "scope": "post-run independent five-state ODE and direct action quadrature",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "raw_sha256": hashlib.sha256((folder / "rows.jsonl").read_bytes()).hexdigest(),
        "q": q,
        "r": r,
        "moments": moments,
        "population_v_decrease": gap,
        "change_from_primary": difference,
        "relative_change": abs(difference / row["population_v_decrease"]),
        **details,
    }
    (folder / "independent_formulation_check.json").write_text(
        json.dumps(result, indent=2)
    )
    print(json.dumps(result, indent=2))
    assert gap > 0 and abs(difference) < 0.01 * row["population_v_decrease"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summarize", type=Path)
    parser.add_argument("--cross-check", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
    if args.summarize:
        summarize(args.summarize)
    if args.cross_check:
        cross_check(args.cross_check)
