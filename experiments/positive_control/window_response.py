"""Exact-plane-wave window response; see docs/paper/paper_a/window_response_theory.md.

The numerical integrals are development/calibration, never holdout posets.
The existing P14 predicate resolves near-boundary vector decisions.
"""

from __future__ import annotations

import functools
import math

import numpy as np
from p14_plane_wave import PlaneWaveGeometry, Slab, causal_relation
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator
from scipy.stats import qmc

BOXES = {"R+": (2.0, 6.0), "R0": (math.sqrt(12), math.sqrt(12)), "R-": (6.0, 2.0)}
REGIONS = (*BOXES, "D")


def box_coefficients(region: str) -> np.ndarray:
    """Analytic p0,c1,c2 at Du=Dv=1, from polynomial disk integrals."""
    if region == "R0":
        return np.array(
            [
                math.pi / 108 + 1 / 5184 - 64 * math.sqrt(6) / 33075,
                0.0,
                257 * math.pi / 48600 + 73 / 435456 - 137204 * math.sqrt(6) / 30405375,
            ]
        )
    if region not in BOXES:
        raise ValueError("Taylor coefficients are for the three declared boxes")
    return np.array(
        [
            math.pi / 108 + 1 / 5184 - 128 * math.sqrt(2) / 33075,
            (1 if region == "R+" else -1)
            * (math.pi / 54 + 1 / 540 - 352 * math.sqrt(2) / 25515),
            1613 * math.pi / 48600
            + 16733 / 2177280
            - 1316008 * math.sqrt(2) / 30405375,
        ]
    )


def coefficient_quadrature(region: str, order: int = 40) -> np.ndarray:
    """Independent Gauss integration before polynomial moment reduction."""
    a, b = BOXES[region]
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes[:, None, None] + 1) / 2
    r = (nodes[None, :, None] + 1) / 2
    theta = (nodes[None, None, :] + 1) * math.pi / 4
    quadrature = (
        weights[:, None, None]
        * weights[None, :, None]
        * weights[None, None, :]
        * math.pi
        / 16
    )
    s = t * t
    x, y = math.sqrt(2) * t * r * np.cos(theta), math.sqrt(2) * t * r * np.sin(theta)
    m2x, m2y = (a - x) ** 2 / 12, (b - y) ** 2 / 12
    c0 = r * r
    c1 = s * ((m2x - m2y) / 2 + (x * x - y * y) / 24)
    c2 = -(s**3) * ((m2x + m2y) / 24 + (x * x + y * y) / 1440)
    c11 = (
        s
        * s
        * (
            ((a - x) ** 4 + (b - y) ** 4) / 320
            - m2x * m2y / 2
            + (x * x - y * y) * (m2x - m2y) / 24
            + (x * x - y * y) ** 2 / 576
        )
    )
    measure = quadrature * 32 * t**3 * (1 - s) * r * (a - x) * (b - y) / (a * a * b * b)
    return np.array(
        [
            np.sum(measure * value)
            for value in ((1 - c0) ** 2 / 2, -(1 - c0) * c1, -(1 - c0) * c2 + c11 / 2)
        ]
    )


def action(s, xp, xq, yp, yq, amplitude: float):
    """Vector geodesic action plus a conservative floating-point guard.

    Domain: study amplitudes |A|<=2, 0<s<=1. The guard is an arithmetic
    model (like the legacy double path), not a directed-rounding proof.
    Near the guard, order decisions use the legacy certified predicate.
    """
    if not math.isfinite(amplitude) or abs(amplitude) > 2:
        raise ValueError("study amplitude must be finite and in [-2,2]")
    s = np.asarray(s, dtype=float)
    if np.any((s <= 0) | (s > 1) | ~np.isfinite(s)):
        raise ValueError("separations must be finite and in (0,1]")
    if amplitude < 0:
        xp, yp, xq, yq = yp, xp, yq, xq
    if amplitude == 0:
        terms = ((xq - xp) ** 2 / (2 * s), (yq - yp) ** 2 / (2 * s))
    else:
        w = math.sqrt(abs(amplitude))
        z = w * s / 2
        terms = (
            w * np.tanh(z) * ((xp + xq) / 2) ** 2,
            w / 4 / np.tanh(z) * (xq - xp) ** 2,
            -w * np.tan(z) * ((yp + yq) / 2) ** 2,
            w / 4 / np.tan(z) * (yq - yp) ** 2,
        )
    return sum(terms), 256 * np.finfo(float).eps * (1 + sum(abs(v) for v in terms))


def v_survival(cost):
    c = np.clip(cost, -1, 1)
    return np.where(c >= 0, (1 - c) ** 2 / 2, 1 - (1 + c) ** 2 / 2)


def box_pair_values(unit: np.ndarray, region: str, amplitude: float):
    """Five-dimensional exact conditional mean, with all v integrated out."""
    a, b = BOXES[region]
    s = 1 - np.sqrt(1 - unit[:, 0])
    dx = a * (1 - np.sqrt(1 - unit[:, 1]))
    dy = b * (1 - np.sqrt(1 - unit[:, 2]))
    mx, my = (a - dx) * unit[:, 3] / 2, (b - dy) * unit[:, 4] / 2
    cost, _ = action(s, mx - dx / 2, mx + dx / 2, my - dy / 2, my + dy / 2, amplitude)
    return v_survival(cost)


def diamond_b(u, amplitude: float):
    u = np.asarray(u)
    if amplitude == 0:
        b = 1 / (2 * u * (1 - u))
        return b, b
    w = math.sqrt(abs(amplitude))
    bx = w / 2 * (1 / np.tanh(w * u) + 1 / np.tanh(w * (1 - u)))
    by = w / 2 * (1 / np.tan(w * u) + 1 / np.tan(w * (1 - u)))
    return (bx, by) if amplitude >= 0 else (by, bx)


def diamond_slice(u, amplitude: float):
    u = np.asarray(u)
    bx, by = diamond_b(np.clip(u, 1e-14, 1 - 1e-14), amplitude)
    return np.where((u > 0) & (u < 1), math.pi / (2 * np.sqrt(bx * by)), 0)


def diamond_moments(amplitude: float, order: int = 64):
    """One/two-dimensional integrals from interval volumes; independent of QMC."""
    if not 0 <= amplitude <= 2:
        raise ValueError("study amplitude outside [0,2]")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    u, weights = (nodes + 1) / 2, weights / 2
    w = math.sqrt(amplitude)

    def factor(z):
        return np.sinc(z / math.pi) * np.where(
            z == 0, 1.0, np.sinh(z) / np.where(z == 0, 1.0, z)
        )

    def j(s):
        z = w * np.asarray(s)[..., None]
        values = (
            2 * u * (1 - u) * np.sqrt(factor(z * u) * factor(z * (1 - u)) / factor(z))
        )
        return np.asarray(s) ** 2 * np.sum(values * weights, axis=-1)

    j1, ju, jv = float(j(1.0)), j(u), j(1 - u)
    density = 2 / math.pi * diamond_slice(u, amplitude)
    p = float(np.sum(weights * density * (ju + jv)) / (6 * j1 * j1))
    both = float(
        np.sum(weights * density * ((ju * ju + jv * jv) / 15 + ju * jv / 45)) / j1**3
    )
    return {"p": p, "zeta": both - p * p, "volume": math.pi * j1 / 2}


@functools.lru_cache(maxsize=32)
def diamond_cdf(amplitude: float, grid: int = 4097):
    u = np.linspace(0, 1, grid)
    cdf = cumulative_simpson(diamond_slice(u, amplitude), x=u, initial=0)
    volume = float(cdf[-1])
    cdf /= volume
    return PchipInterpolator(cdf, u), volume


def diamond_backbone(u, radius_unit, angle_unit, amplitude: float):
    bx, by = diamond_b(u, amplitude)
    r = np.sqrt(1 - np.sqrt(1 - radius_unit))
    theta = 2 * math.pi * angle_unit
    x, y = r * np.cos(theta) / np.sqrt(bx), r * np.sin(theta) / np.sqrt(by)
    c_p, _ = action(u, 0, x, 0, y, amplitude)
    c_q, _ = action(1 - u, x, 0, y, 0, amplitude)
    return np.column_stack((u, x, y, -0.5 + c_q, 0.5 - c_p))


def backbone(unit: np.ndarray, region: str, amplitude: float, grid: int = 4097):
    """Return u,x,y,v_low,v_high; v is analytically marginalized."""
    if region == "D":
        inv, _ = diamond_cdf(amplitude, grid)
        return diamond_backbone(inv(unit[:, 0]), unit[:, 1], unit[:, 2], amplitude)
    a, b = BOXES[region]
    return np.column_stack(
        (
            unit[:, 0],
            a * (unit[:, 1] - 0.5),
            b * (unit[:, 2] - 0.5),
            np.full(len(unit), -0.5),
            np.full(len(unit), 0.5),
        )
    )


def sample_points(region: str, amplitude: float, n: int, rng):
    """Fixed-N physical-volume sample; D uses rejection, not inverse CDF."""
    if n < 2 or region not in REGIONS or not 0 <= amplitude <= 2:
        raise ValueError("invalid study population")
    if region == "D":
        accepted, count = [], 0
        ceiling = float(diamond_slice(0.5, amplitude))
        while count < n:
            u = rng.random(max(16, 2 * (n - count)))
            keep = rng.random(len(u)) * ceiling < diamond_slice(u, amplitude)
            accepted.append(u[keep])
            count += int(keep.sum())
        u = np.concatenate(accepted)[:n]
        z = rng.random((n, 3))
        bb = diamond_backbone(u, z[:, 0], z[:, 1], amplitude)
        v = bb[:, 3] + z[:, 2] * (bb[:, 4] - bb[:, 3])
    else:
        z = rng.random((n, 4))
        bb = backbone(z[:, :3], region, amplitude)
        v = z[:, 3] - 0.5
    return np.column_stack((bb[:, 0], v, bb[:, 1], bb[:, 2]))


def conditional_line(anchor, other, amplitude: float):
    cost, _ = action(
        np.abs(other[:, 0] - anchor[:, 0]),
        anchor[:, 1],
        other[:, 1],
        anchor[:, 2],
        other[:, 2],
        amplitude,
    )
    future = other[:, 0] > anchor[:, 0]
    width = other[:, 4] - other[:, 3]
    slope = np.where(future, 1, -1) / width
    intercept = np.where(future, -cost - other[:, 3], other[:, 4] - cost) / width
    return slope, intercept


def integrate_v(lo, hi, a, b, c, d):
    """Means of two clipped lines and their product on uniform [lo,hi]."""
    if np.any(np.asarray(hi) <= lo) or np.any(a == 0) or np.any(c == 0):
        raise ValueError("positive intervals and nonzero slopes required")
    roots = np.column_stack((lo, hi, -b / a, (1 - b) / a, -d / c, (1 - d) / c))
    knots = np.sort(np.clip(roots, np.asarray(lo)[:, None], np.asarray(hi)[:, None]))
    left, right = knots[:, :-1], knots[:, 1:]
    mid = (left + right) / 2
    f = lambda x: np.clip(a[:, None] * x + b[:, None], 0, 1)  # noqa: E731
    g = lambda x: np.clip(c[:, None] * x + d[:, None], 0, 1)  # noqa: E731
    weights = (right - left) / (hi - lo)[:, None]
    p = np.sum(weights * (f(left) + f(right) + g(left) + g(right)) / 4, axis=1)
    both = np.sum(
        weights * (f(left) * g(left) + 4 * f(mid) * g(mid) + f(right) * g(right)) / 6,
        axis=1,
    )
    return p, both


def triple_values(unit: np.ndarray, region: str, amplitude: float, grid: int = 4097):
    anchor, second, third = [
        backbone(unit[:, i : i + 3], region, amplitude, grid) for i in (0, 3, 6)
    ]
    a, b = conditional_line(anchor, second, amplitude)
    c, d = conditional_line(anchor, third, amplitude)
    return integrate_v(anchor[:, 3], anchor[:, 4], a, b, c, d)


def quadrature_replicates(
    region: str, amplitude: float, power: int, reps: int, seed: int, grid: int = 4097
):
    result = []
    for rep in range(reps):
        engine = qmc.Sobol(
            9, scramble=True, bits=52, rng=np.random.default_rng([*seed, rep])
        )
        values = engine.random_base2(power)
        # Chunking retains every Sobol point and keeps memory bounded.
        totals = np.zeros(2)
        for chunk in np.array_split(values, max(1, len(values) // 32768)):
            p, both = triple_values(chunk, region, amplitude, grid)
            totals += [p.sum(), both.sum()]
        result.append(totals / len(values))
    return np.array(result)


@functools.lru_cache(maxsize=8)
def _pair_indices(n):
    return np.triu_indices(n, 1)


def order_matrix(points: np.ndarray, amplitude: float):
    """u-sorted transitive relation matrix, legacy fallback and ambiguity count."""
    if points.ndim != 2 or points.shape[1] != 4 or not np.isfinite(points).all():
        raise ValueError("finite (n,4) points required")
    pts = points[np.lexsort((-points[:, 1], points[:, 0]))]
    i, j = _pair_indices(len(pts))
    s = pts[j, 0] - pts[i, 0]
    cost, error = action(
        np.maximum(s, 1e-300), pts[i, 2], pts[j, 2], pts[i, 3], pts[j, 3], amplitude
    )
    dv = pts[j, 1] - pts[i, 1]
    margin = -(dv + cost)
    guard = error + 256 * np.finfo(float).eps * (1 + abs(dv))
    near = (abs(margin) <= guard) | (s == 0)
    related = margin > 0
    geometry = PlaneWaveGeometry(Slab(1, 8, 16, 16), math.sqrt(abs(amplitude)))
    for k in np.flatnonzero(near):
        p, q = pts[i[k]].copy(), pts[j[k]].copy()
        if amplitude < 0:
            p[[2, 3]], q[[2, 3]] = p[[3, 2]], q[[3, 2]]
        verdict = causal_relation(geometry, p, q)
        if verdict.related is None:
            raise ArithmeticError("unresolved pair: stop rather than drop it")
        related[k] = verdict.related
    matrix = np.zeros((len(pts), len(pts)), dtype=bool)
    matrix[i, j] = related
    return matrix, int(near.sum())


def statistics(points: np.ndarray, amplitude: float):
    matrix, escalations = order_matrix(points, amplitude)
    n = len(matrix)
    return {
        "f": float(matrix.sum() / math.comb(n, 2)),
        "chain3": float(
            np.dot(matrix.sum(axis=0), matrix.sum(axis=1)) / math.comb(n, 3)
        )
        if n >= 3
        else 0.0,
        "escalations": escalations,
    }
