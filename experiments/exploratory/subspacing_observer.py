"""Blind exploratory observer: only strict-order callbacks and optional q.

No generator imports, point labels, truth coordinates or per-probe chain
U-statistic. Reconstructed coordinates are nuisance parameters of the background
model. Multistart initialization mitigates, but cannot certify absence of, local
minima. Scores and stability diagnostics are not confirmatory p-values.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import czt, fftconvolve
from scipy.special import zeta

PAST = np.array([-0.5, 0.0, 0.0, 0.0])
FUTURE = np.array([0.5, 1.0, 0.0, 0.0])
# alpha(L) = L^2 sum_n zeta(2n) (q L^2)^(n-1) / pi^(2n).
_ALPHA_POWERS = 2 * np.arange(1, 13)
_ALPHA_SERIES = zeta(_ALPHA_POWERS) / np.pi**_ALPHA_POWERS


def _alpha_q(model, length):
    """Analytic d/dq of unnormalized alpha(length)."""
    length = np.asarray(length)
    x2 = model.q * length**2
    powers = np.arange(1, len(_ALPHA_SERIES))
    series = length**4 * np.polynomial.polynomial.polyval(
        x2, powers * _ALPHA_SERIES[1:]
    )
    if model.q == 0:
        return series
    return np.where(
        x2 < 0.25,
        series,
        (length * model.alpha_prime(length) - 2 * model.alpha(length)) / (2 * model.q),
    )


def _m_q(model, u):
    """Analytic d/dq of m(u) = tan(sqrt(q)*u)/sqrt(q)."""
    u = np.asarray(u)
    x2 = model.q * u**2
    series = u**3 * np.polynomial.polynomial.polyval(
        x2, [1 / 3, 4 / 15, 17 / 105, 248 / 2835, 1382 / 31185]
    )
    if model.q == 0:
        return series
    return np.where(
        x2 < 1e-4,
        series,
        (u / np.cos(model.k * u) ** 2 - model.m(u)) / (2 * model.q),
    )


def _interval_dq(model, first, second):
    """Analytic d/dq of normalized interval probability model.interval."""
    first, second = np.broadcast_arrays(first, second)
    length = np.maximum(second[..., 0] - first[..., 0], 0.0)
    dm = model.m(second[..., 0]) - model.m(first[..., 0])
    safe = np.where(length > 0, dm, 1.0)
    radius = np.sum((second[..., 2:] - first[..., 2:]) ** 2, axis=-1)
    kap = np.where(
        length > 0,
        np.maximum(second[..., 1] - first[..., 1] - radius / (2 * safe), 0.0),
        0.0,
    )
    dq = _m_q(model, second[..., 0]) - _m_q(model, first[..., 0])
    alpha = model.alpha(length)
    return (
        _alpha_q(model, length) / model.gamma
        - alpha * _alpha_q(model, 1.0) / model.gamma**2
    ) * kap**2 + (alpha / model.gamma) * kap * radius * dq / safe**2


class BackgroundModel:
    def __init__(self, q):
        if not np.isfinite(q) or not 0 <= q < np.pi**2:
            raise ValueError("q must lie in [0,pi^2)")
        self.q, self.k = float(q), np.sqrt(q)
        self.gamma = self.alpha(1.0)

    def alpha(self, length):
        """(1 - x cot(x))/(2q); 12 series terms avoid cancellation for |x|<.5."""
        length = np.asarray(length)
        if not self.q:
            return length**2 / 6
        x = self.k * length
        small = np.abs(x) < 0.5
        series = length**2 * np.polynomial.polynomial.polyval(x * x, _ALPHA_SERIES)
        safe = np.where(small, 1.0, x)
        exact = (1 - safe / np.tan(safe)) / (2 * np.where(small, 1.0, self.q))
        return np.where(small, series, exact)

    def alpha_prime(self, length):
        length = np.asarray(length)
        if not self.q:
            return length / 3
        x = self.k * length
        small = np.abs(x) < 0.5
        series = length * np.polynomial.polynomial.polyval(
            x * x, _ALPHA_POWERS * _ALPHA_SERIES
        )
        safe = np.where(small, 1.0, x)
        exact = (safe - np.sin(safe) * np.cos(safe)) / (2 * self.k * np.sin(safe) ** 2)
        return np.where(small, series, exact)

    def m(self, u):
        return np.tan(self.k * np.asarray(u)) / self.k if self.q else np.asarray(u)

    def interval(self, first, second, jacobian=False):
        """Probability and derivatives with respect to BOTH endpoint positions."""
        first, second = np.broadcast_arrays(first, second)
        length = second[..., 0] - first[..., 0]
        dm = self.m(second[..., 0]) - self.m(first[..., 0])
        safe = np.where(length > 0, dm, 1.0)
        dy = second[..., 2:] - first[..., 2:]
        radius = np.sum(dy * dy, axis=-1)
        margin = second[..., 1] - first[..., 1] - radius / (2 * safe)
        active = (length > 0) & (margin > 0)
        kap = np.where(active, margin, 0)
        alpha = self.alpha(np.maximum(length, 0)) / self.gamma
        value = alpha * kap**2
        if not jacobian:
            return value
        left, right = np.zeros(first.shape), np.zeros(second.shape)
        grad_k = 2 * alpha * kap
        grad_l = self.alpha_prime(np.maximum(length, 0)) / self.gamma * kap**2
        left[..., 0] = (
            -grad_l
            - grad_k * radius / (2 * safe**2) / np.cos(self.k * first[..., 0]) ** 2
        )
        right[..., 0] = (
            grad_l
            + grad_k * radius / (2 * safe**2) / np.cos(self.k * second[..., 0]) ** 2
        )
        left[..., 1], right[..., 1] = -grad_k, grad_k
        left[..., 2:] = grad_k[..., None] * dy / safe[..., None]
        right[..., 2:] = -left[..., 2:]
        return value, left, right

    def moments(self, points, early, late, jacobian=False):
        """Interval probabilities with only the moving probe endpoint differentiated."""
        points = np.atleast_2d(points)
        anchors = np.vstack((PAST, FUTURE, early, late))
        sign = np.r_[1, -1, np.ones(len(early)), -np.ones(len(late))]
        length = sign * (points[:, :1] - anchors[:, 0])
        dm = sign * (self.m(points[:, :1]) - self.m(anchors[:, 0]))
        safe = np.where(length > 0, dm, 1.0)
        dy = points[:, None, 2:] - anchors[:, 2:]
        radius = np.sum(dy * dy, axis=-1)
        margin = sign * (points[:, 1:2] - anchors[:, 1]) - radius / (2 * safe)
        kap = np.where(length > 0, np.maximum(margin, 0), 0)
        alpha = self.alpha(np.maximum(length, 0)) / self.gamma
        values = alpha * kap**2
        if not jacobian:
            return values
        grad_k = 2 * alpha * kap
        grad_l = self.alpha_prime(np.maximum(length, 0)) / self.gamma * kap**2
        jac = np.empty(values.shape + (4,))
        jac[..., 0] = sign * (
            grad_l
            + grad_k * radius / (2 * safe**2) / np.cos(self.k * points[:, :1]) ** 2
        )
        jac[..., 1] = sign * grad_k
        jac[..., 2:] = -grad_k[..., None] * dy / safe[..., None]
        return values, jac

    def shear_response(self, points, early, late):
        """Exact differential of interval probabilities to endpoint STF shear.

        beta is represented by (B11+i B12). The local oscillatory endpoint
        carrier is factored out; the remaining sec^2(ku) is its principal symbol.
        Thus the Fourier matcher uses a high-frequency symbol, while truth
        population probabilities can use the full finite-lambda integral.
        """
        points = np.atleast_2d(points)
        responses = []
        for anchor, endpoint in (
            [(PAST, 1), (FUTURE, -1)]
            + [(a, 1) for a in early]
            + [(b, -1) for b in late]
        ):
            first, second = (anchor, points) if endpoint == 1 else (points, anchor)
            length = second[..., 0] - first[..., 0]
            dm = self.m(second[..., 0]) - self.m(first[..., 0])
            dy = second[..., 2:] - first[..., 2:]
            safe = np.where(length > 0, dm, 1)
            kap = np.maximum(
                second[..., 1] - first[..., 1] - np.sum(dy * dy, axis=-1) / (2 * safe),
                0,
            )
            coeff = -self.alpha(np.maximum(length, 0)) / self.gamma * kap / safe**2
            responses.append(
                endpoint
                * coeff
                * (dy[..., 0] - 1j * dy[..., 1]) ** 2
                / np.cos(self.k * points[:, 0]) ** 2
            )
        return np.stack(responses, axis=-1)


def order_counts(
    query, references, probes, early=(), late=(), tile=4096, count_query=None
):
    """Tiled counts of tip intervals and anchor intervals, with no a_X."""
    refs, probes = np.asarray(references), np.asarray(probes)
    early, late = np.asarray(early, dtype=int), np.asarray(late, dtype=int)
    if not len(refs) or np.intersect1d(refs, probes, assume_unique=True).size:
        raise ValueError("nonempty reference split must be disjoint from probes")
    if count_query is not None:
        return count_query(refs, probes, early, late)
    sums = np.zeros((len(probes), 2 + len(early) + len(late)))
    for start in range(0, len(refs), tile):
        rr = refs[start : start + tile]
        ar = query(early, rr)[0].astype(float) if len(early) else None
        rb = query(rr, late)[0].astype(float) if len(late) else None
        for j in range(0, len(probes), 512):
            fw, bw = query(rr, probes[j : j + 512])
            target = sums[j : j + fw.shape[1]]
            target[:, 0] += fw.sum(axis=0)
            target[:, 1] += bw.sum(axis=0)
            if ar is not None:
                target[:, 2 : 2 + len(early)] += (ar @ fw).T
            if rb is not None:
                target[:, 2 + len(early) :] += bw.T @ rb
    return sums / len(refs)


def select_anchors(query, pilot, count=4, count_query=None, calibration_count=None):
    """Ordered rings with compact ``count``-anchor cores first on each side.

    Without calibration_count, return only the two cores, as before.
    Auxiliary anchors improve calibration without restricting the probe pool.
    """
    pilot = np.asarray(pilot)
    candidates, reference = pilot[: len(pilot) // 3], pilot[len(pilot) // 3 :]
    pf = order_counts(query, reference, candidates, count_query=count_query)
    p, f = pf.T
    early = candidates[(p >= 0.0005) & (p <= 0.04) & (f >= 0.25) & (f <= 0.80)][:96]
    late = candidates[(f >= 0.0005) & (f <= 0.04) & (p >= 0.25) & (p <= 0.80)][:96]
    if min(len(early), len(late)) < count:
        raise ValueError("insufficient pilot anchor candidates")
    cross = query(early, late)[0]
    volumes = order_counts(
        query, reference, late, early=early, count_query=count_query
    )[:, 2:].T
    ar = query(early, reference)[0]
    rb = query(reference, late)[0]
    support = np.ones(len(reference), dtype=bool)
    chosen, available = [], np.ones(len(late), dtype=bool)
    for _ in range(count):
        eligible = [
            i
            for i in range(len(early))
            if i not in chosen and np.count_nonzero(available & cross[i]) >= count
        ]
        if not eligible:
            raise ValueError("no complete early/late anchor family")
        if not chosen:
            pick = max(eligible, key=lambda i: np.count_nonzero(ar[i]))
        else:
            pick = max(
                eligible,
                key=lambda i: (
                    min(np.linalg.norm(volumes[i] - volumes[j]) for j in chosen)
                    * np.count_nonzero(support & ar[i]) ** 3
                ),
            )
        chosen.append(pick)
        available &= cross[pick]
        support &= ar[pick]
    eligible = list(np.flatnonzero(available))
    selected = [max(eligible, key=lambda i: np.count_nonzero(support & rb[:, i]))]
    support &= rb[:, selected[0]]
    while len(selected) < count:
        selected.append(
            max(
                (i for i in eligible if i not in selected),
                key=lambda i: (
                    min(
                        np.linalg.norm(volumes[chosen, i] - volumes[chosen, j])
                        for j in selected
                    )
                    * np.count_nonzero(support & rb[:, i]) ** 3
                ),
            )
        )
        support &= rb[:, selected[-1]]
    core_a, core_b = early[chosen], late[selected]
    if calibration_count is None or calibration_count <= count:
        return core_a, core_b
    # Broad auxiliary ring improves conditioning without shrinking core support.
    broad_a = candidates[(p >= 0.0005) & (p <= 0.12) & (f >= 0.12) & (f <= 0.80)][:192]
    broad_b = candidates[(f >= 0.0005) & (f <= 0.12) & (p >= 0.12) & (p <= 0.80)][:192]
    broad_a = np.unique(np.r_[core_a, broad_a])
    broad_b = np.unique(np.r_[core_b, broad_b])
    cross_mass = order_counts(
        query, reference, broad_b, early=broad_a, count_query=count_query
    )[:, 2:].T
    a_index = {int(idx): i for i, idx in enumerate(broad_a)}
    b_index = {int(idx): i for i, idx in enumerate(broad_b)}
    aa, bb = list(core_a), list(core_b)

    while len(aa) < calibration_count or len(bb) < calibration_count:
        changed = False
        for ring, other, pool in ((aa, bb, broad_a), (bb, aa, broad_b)):
            if len(ring) >= calibration_count:
                continue
            eligible = np.array([idx for idx in pool if idx not in ring], dtype=int)
            ordered = (
                query(eligible, other)[0].all(axis=1)
                if ring is aa
                else query(other, eligible)[0].all(axis=0)
            )
            supported = np.array(
                [
                    min(
                        cross_mass[a_index[int(idx)], b_index[int(j)]]
                        if ring is aa
                        else cross_mass[a_index[int(j)], b_index[int(idx)]]
                        for j in other
                    )
                    > 0.0005
                    for idx in eligible
                ],
                dtype=bool,
            )
            eligible = eligible[ordered & supported]
            if len(eligible):
                profiles = np.sqrt(cross_mass if ring is aa else cross_mass.T)
                index = a_index if ring is aa else b_index
                pick = max(
                    eligible,
                    key=lambda idx: min(
                        np.linalg.norm(
                            profiles[index[int(idx)]] - profiles[index[int(j)]]
                        )
                        for j in ring
                    ),
                )
                ring.append(int(pick))
                changed = True
        if not changed:
            break
    return np.asarray(aa), np.asarray(bb)


def flat_gram_initialization(pf, cross):
    """Minkowski cross-Gram SVD plus spatial metric completion from norms."""
    na, nb = cross.shape
    t = (1 + np.sqrt(pf[:, 0]) - np.sqrt(pf[:, 1])) / 2
    radii = np.maximum(t * t - np.sqrt(pf[:, 0]), 1e-8)
    gram = (
        radii[:na, None]
        + radii[None, na:]
        - (t[:na, None] - t[None, na:]) ** 2
        + np.sqrt(cross)
    ) / 2
    u, singular, vh = np.linalg.svd(gram, full_matrices=False)
    a = u[:, :3] * np.sqrt(singular[:3])
    b = vh[:3].T * np.sqrt(singular[:3])

    def unpack(z):
        return np.array(
            [[np.exp(z[0]), 0, 0], [z[3], np.exp(z[1]), 0], [z[4], z[5], np.exp(z[2])]]
        )

    def residual(z):
        transform = unpack(z)
        aa, bb = a @ transform, np.linalg.solve(transform, b.T).T
        return np.r_[np.sum(aa * aa, axis=1), np.sum(bb * bb, axis=1)] - radii

    fit = least_squares(
        residual,
        np.zeros(6),
        max_nfev=100,
        bounds=(-5, 5),
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    transform = unpack(fit.x)
    return t, np.vstack((a @ transform, np.linalg.solve(transform, b.T).T))


def _gauge(points):
    points = points.copy()
    angle = np.arctan2(points[0, 3], points[0, 2])
    c, s = np.cos(angle), np.sin(angle)
    points[:, 2:] = points[:, 2:] @ np.array([[c, -s], [s, c]])
    if points[1, 3] < 0:
        points[:, 3] *= -1
    points[0, 3] = 0
    return points


@dataclass
class Calibration:
    q: float
    early: np.ndarray
    late: np.ndarray
    rms: float
    condition: float
    evaluations: int


def _calibration_system(pf, cross, q):
    """Shared weighted equations, exact Jacobian, bounds, and pack/unpack."""
    na, nb = cross.shape
    target = np.r_[pf.ravel(), cross.ravel()]
    inv_w = np.minimum(1 / np.sqrt(np.maximum(target * (1 - target), 1e-6)), 100.0)
    keep = np.arange(4 * (na + nb)) != 3

    def pack(points, qq):
        flat = points.ravel()[keep]
        return np.r_[flat, qq] if q is None else flat

    def unpack(z):
        points = np.zeros(4 * (na + nb))
        points[keep] = z[: np.count_nonzero(keep)]
        return points.reshape(-1, 4), float(z[-1]) if q is None else float(q)

    def equations(z, jac=False):
        points, qq = unpack(z)
        model = BackgroundModel(qq)
        if not jac:
            pred = np.column_stack(
                (model.interval(PAST, points), model.interval(points, FUTURE))
            )
            return (
                np.r_[
                    pred.ravel(),
                    model.interval(points[:na, None], points[None, na:]).ravel(),
                ]
                - target
            ) * inv_w
        _, _, dp = model.interval(PAST, points, True)
        _, df, _ = model.interval(points, FUTURE, True)
        _, da, db = model.interval(points[:na, None], points[None, na:], True)
        full = np.zeros((len(target), 4 * (na + nb)))
        idx_a = np.arange(na + nb)
        cols_a = 4 * idx_a[:, None] + np.arange(4)
        full[2 * idx_a[:, None], cols_a] = dp
        full[2 * idx_a[:, None] + 1, cols_a] = df
        rows_c = 2 * (na + nb) + np.arange(na * nb)
        ci = np.repeat(np.arange(na), nb)
        cj = np.tile(np.arange(nb), na)
        full[rows_c[:, None], 4 * ci[:, None] + np.arange(4)] = da.reshape(-1, 4)
        full[rows_c[:, None], 4 * (na + cj)[:, None] + np.arange(4)] = db.reshape(-1, 4)
        result = full[:, keep] * inv_w[:, None]
        if q is None:
            dqp = np.column_stack(
                (_interval_dq(model, PAST, points), _interval_dq(model, points, FUTURE))
            ).ravel()
            dqc = _interval_dq(model, points[:na, None], points[None, na:]).ravel()
            result = np.column_stack((result, np.r_[dqp, dqc] * inv_w))
        return result

    lo = np.tile([-0.4999, -0.5, -2, -2], na + nb)[keep]
    hi = np.tile([0.4999, 1.5, 2, 2], na + nb)[keep]
    if q is None:
        lo, hi = np.r_[lo, 0], np.r_[hi, np.pi**2 - 1e-4]
    return inv_w, pack, unpack, equations, lo, hi


def calibrate(pf, cross, q=None, orientations=64, starts=16, max_nfev=180):
    """Joint exact-background fit; q=None never receives a truth q initializer."""
    pf, cross = np.asarray(pf), np.asarray(cross)
    na, nb = cross.shape
    if pf.shape != (na + nb, 2) or min(na, nb) < 4:
        raise ValueError("at least four anchors per side and matching counts required")
    if np.any(pf <= 0) or np.any(cross <= 0):
        raise ValueError("calibration requires positive interval counts")
    t, spatial = flat_gram_initialization(pf, cross)
    inv_w, pack, unpack, equations, lo, hi = _calibration_system(pf, cross, q)
    guesses = [0.05, 0.5, 1.5, 4.0, 7.0] if q is None else [q]
    buckets = {guess: [] for guess in guesses}
    for i in range(orientations):
        zaxis = 1 - 2 * (i + 0.5) / orientations
        azimuth = i * np.pi * (3 - np.sqrt(5))
        axis = np.array(
            [
                np.sqrt(1 - zaxis * zaxis) * np.cos(azimuth),
                np.sqrt(1 - zaxis * zaxis) * np.sin(azimuth),
                zaxis,
            ]
        )
        transverse = np.cross(axis, [0, 0, 1])
        transverse /= np.linalg.norm(transverse)
        for sign in (1.0, -1.0):
            rotated = spatial @ np.column_stack(
                (axis, transverse, sign * np.cross(axis, transverse))
            )
            points = _gauge(
                np.column_stack(
                    (
                        t - 0.5 + rotated[:, 0],
                        t - rotated[:, 0],
                        np.sqrt(2) * rotated[:, 1:],
                    )
                )
            )
            if q is None or q == 0:
                flat_seed = np.clip(pack(points, 0), lo, hi)
                flat_error = equations(flat_seed) / inv_w
                if np.max(np.abs(flat_error)) < 1e-12:
                    points, _ = unpack(flat_seed)
                    singular = np.linalg.svd(
                        equations(flat_seed, True), compute_uv=False
                    )
                    return Calibration(
                        0.0,
                        points[:na],
                        points[na:],
                        float(np.sqrt(np.mean(flat_error**2))),
                        float(singular[0] / max(singular[-1], 1e-30)),
                        0,
                    )
            for guess in guesses:
                seed = np.clip(pack(points, guess), lo + 1e-8, hi - 1e-8)
                # _gauge can identify the two transverse reflections exactly.
                if not any(
                    np.allclose(seed, old, atol=1e-12, rtol=0)
                    for _, old in buckets[guess][-1:]
                ):
                    buckets[guess].append((np.linalg.norm(equations(seed)), seed))
    # Equal opportunity per q bucket; a global flat-biased ranking starves large q.
    per_bucket = max(1, starts)
    for bucket in buckets.values():
        bucket.sort(key=lambda item: item[0])
    seeds = [
        buckets[guess][rank]
        for rank in range(min(per_bucket, min(map(len, buckets.values()))))
        for guess in guesses
    ]
    best, evaluations = None, 0
    for _, seed in seeds:
        fit = least_squares(
            equations,
            seed,
            jac=lambda z: equations(z, True),
            bounds=(lo, hi),
            max_nfev=max_nfev,
            ftol=1e-14,
            xtol=1e-14,
            gtol=1e-15,
        )
        evaluations += fit.nfev
        if best is None or np.linalg.norm(fit.fun) < np.linalg.norm(best.fun):
            best = fit
        if np.max(np.abs(fit.fun)) < 1e-14:
            break
    points, qq = unpack(best.x)
    points = _gauge(points)
    singular = np.linalg.svd(best.jac, compute_uv=False)
    return Calibration(
        qq,
        points[:na],
        points[na:],
        float(np.sqrt(np.mean((best.fun / inv_w) ** 2))),
        float(singular[0] / max(singular[-1], 1e-30)),
        evaluations,
    )


def _polish_outer_anchors(pf_sub, opp_counts, opp_anchors, opp_sign, qq):
    """Initialize outer anchors via lifted 1D u-grid against core and run damped GN."""
    model = BackgroundModel(qq)
    all_a = np.vstack((PAST, FUTURE, opp_anchors))
    sgn = np.r_[1, -1, np.full(len(opp_anchors), opp_sign)]
    mu = np.column_stack((pf_sub, opp_counts))
    best = np.full(len(mu), np.inf)
    points = np.empty((len(mu), 4))
    u_lo = (
        -0.49 if opp_sign < 0 else min(0.47, float(np.max(opp_anchors[:, 0])) + 0.005)
    )
    u_hi = (
        max(-0.47, float(np.min(opp_anchors[:, 0])) - 0.005) if opp_sign < 0 else 0.49
    )
    if u_lo >= u_hi:
        u_lo, u_hi = (-0.49, 0.0) if opp_sign < 0 else (0.0, 0.49)
    for u in np.linspace(u_lo, u_hi, 65):
        length = sgn * (u - all_a[:, 0])
        dm = sgn * (model.m(u) - model.m(all_a[:, 0]))
        safe = np.where(length > 0, dm, 1.0)
        kap = np.where(
            length > 0,
            np.sqrt(
                np.maximum(
                    model.gamma
                    * mu
                    / np.maximum(model.alpha(np.maximum(length, 0)), 1e-30),
                    0.0,
                )
            ),
            0.0,
        )
        mat = np.column_stack((sgn, all_a[:, 2:] / safe[:, None], -0.5 / safe))
        con = sgn * all_a[:, 1] + np.sum(all_a[:, 2:] ** 2, axis=1) / (2 * safe)
        valid = length > 1e-4
        if np.count_nonzero(valid) < 4:
            continue
        lifted = (kap[:, valid] + con[valid]) @ np.linalg.pinv(mat[valid]).T
        cand = np.column_stack(
            (np.full(len(mu), u), np.clip(lifted[:, 0], 1e-4, 1 - 1e-4), lifted[:, 1:3])
        )
        rad = np.sum((cand[:, None, 2:] - all_a[:, 2:]) ** 2, axis=-1)
        mar = sgn * (cand[:, 1:2] - all_a[:, 1]) - rad / (2 * safe)
        pred = np.where(
            length > 0,
            model.alpha(np.maximum(length, 0)) / model.gamma * np.maximum(mar, 0) ** 2,
            0.0,
        )
        loss = np.sum((mu[:, valid] - pred[:, valid]) ** 2, axis=1)
        good = loss < best
        points[good], best[good] = cand[good], loss[good]
    inv_w = np.minimum(1 / np.sqrt(np.maximum(mu * (1 - mu), 1e-6)), 100.0)
    early_a = opp_anchors if opp_sign > 0 else np.empty((0, 4))
    late_a = opp_anchors if opp_sign < 0 else np.empty((0, 4))
    for _ in range(25):
        pred, jac = model.moments(points, early_a, late_a, True)
        res = (mu - pred) * inv_w
        jw = jac * inv_w[..., None]
        gram = np.einsum("nki,nkj->nij", jw, jw)
        rhs = np.einsum("nki,nk->ni", jw, res)
        damp = 1e-10 + 1e-7 * np.trace(gram, axis1=1, axis2=2)
        d = np.linalg.solve(gram + damp[:, None, None] * np.eye(4), rhs[..., None])[
            ..., 0
        ]
        d /= np.maximum(1, np.max(np.abs(d), axis=1) / 0.15)[:, None]
        loss = np.sum(res**2, axis=1)
        for st in (1.0, 0.5, 0.25, 0.125):
            pend = np.flatnonzero(np.any(d != 0, axis=1))
            if not len(pend):
                break
            prop = points[pend] + st * d[pend]
            prop[:, 0] = np.clip(prop[:, 0], -0.499, 0.499)
            prop[:, 1] = np.clip(prop[:, 1], 0.001, 0.999)
            npd = model.moments(prop, early_a, late_a)
            good = np.sum(((mu[pend] - npd) * inv_w[pend]) ** 2, axis=1) < loss[pend]
            acc = pend[good]
            points[acc] = prop[good]
            d[acc] = 0
    return points


def calibrate_progressive(pf, cross, q=None, n_ref=None, core=5, **kwargs):
    """Run calibration, escalating to hierarchical core->ring fit on high RMS."""
    pf, cross = np.asarray(pf), np.asarray(cross)
    cal = calibrate(pf, cross, q=q, **kwargs)
    bad = (cal.rms * np.sqrt(n_ref) >= 0.35) if n_ref is not None else (cal.rms > 2e-4)
    na, nb = cross.shape
    if not bad or min(na, nb) <= core:
        return cal
    sub_idx = np.r_[:core, na : na + core]
    c0 = calibrate(
        pf[sub_idx],
        cross[:core, :core],
        q=q,
        starts=max(int(kwargs.get("starts", 1)), 4),
        max_nfev=max(int(kwargs.get("max_nfev", 180)), 300),
    )
    ee = _polish_outer_anchors(pf[core:na], cross[core:na, :core], c0.late, -1.0, c0.q)
    el = _polish_outer_anchors(
        pf[na + core :], cross[:core, core:].T, c0.early, +1.0, c0.q
    )
    pts0 = _gauge(np.vstack((c0.early, ee, c0.late, el)))
    inv_w, pack, unpack, equations, lo, hi = _calibration_system(pf, cross, q)
    z0 = pack(pts0, c0.q)
    fit = least_squares(
        equations,
        np.clip(z0, lo + 1e-8, hi - 1e-8),
        jac=lambda z: equations(z, True),
        bounds=(lo, hi),
        max_nfev=300,
        ftol=1e-14,
        xtol=1e-14,
        gtol=1e-15,
    )
    rms = float(np.sqrt(np.mean((fit.fun / inv_w) ** 2)))
    if rms < cal.rms:
        pts, qq = unpack(fit.x)
        pts = _gauge(pts)
        singular = np.linalg.svd(fit.jac, compute_uv=False)
        return Calibration(
            qq,
            pts[:na],
            pts[na:],
            rms,
            float(singular[0] / max(singular[-1], 1e-30)),
            c0.evaluations + cal.evaluations + fit.nfev,
        )
    return cal


def linear_probe_initialization(mu, calibration, grid=33):
    """Exact lifted [v,y1,y2,|y|^2] linear solve over a one-dimensional u grid."""
    model = BackgroundModel(calibration.q)
    early, late = calibration.early, calibration.late
    anchors = np.vstack((PAST, FUTURE, early, late))
    sign = np.r_[1, -1, np.ones(len(early)), -np.ones(len(late))]
    lower = max(np.max(early[:, 0]) + 1e-4, -0.49)
    upper = min(np.min(late[:, 0]) - 1e-4, 0.49)
    if lower >= upper:
        raise ValueError("calibrated anchors leave no ordered probe u interval")
    best = np.full(len(mu), np.inf)
    points = np.empty((len(mu), 4))
    for u in np.linspace(lower, upper, grid):
        length = sign * (u - anchors[:, 0])
        dm = sign * (model.m(u) - model.m(anchors[:, 0]))
        kappa = np.sqrt(
            np.maximum(model.gamma * mu / np.maximum(model.alpha(length), 1e-30), 0)
        )
        matrix = np.column_stack((sign, anchors[:, 2:] / dm[:, None], -0.5 / dm))
        const = sign * anchors[:, 1] + np.sum(anchors[:, 2:] ** 2, axis=1) / (2 * dm)
        lifted = (kappa + const) @ np.linalg.pinv(matrix).T
        candidate = np.column_stack(
            (np.full(len(mu), u), np.clip(lifted[:, 0], 1e-4, 1 - 1e-4), lifted[:, 1:3])
        )
        # Exact model.moments at fixed u: reuse its channel coefficients for all probes.
        radius = np.sum((candidate[:, None, 2:] - anchors[:, 2:]) ** 2, axis=-1)
        margin = sign * (candidate[:, 1:2] - anchors[:, 1]) - radius / (2 * dm)
        predicted = model.alpha(length) / model.gamma * np.maximum(margin, 0) ** 2
        loss = np.sum((mu - predicted) ** 2, axis=1)
        good = loss < best
        points[good], best[good] = candidate[good], loss[good]
    return points


def triangulate(mu, calibration, iterations=10, initial=None):
    """Batched damped Gauss-Newton, using analytic four-coordinate Jacobians."""
    model = BackgroundModel(calibration.q)
    mu = np.asarray(mu)
    if initial is None:
        points = linear_probe_initialization(mu, calibration)
    else:
        points = np.asarray(initial).copy()
    for _ in range(iterations):
        pred, jac = model.moments(points, calibration.early, calibration.late, True)
        residual = mu - pred
        gram = np.einsum("nki,nkj->nij", jac, jac)
        rhs = np.einsum("nki,nk->ni", jac, residual)
        damping = 1e-10 + 1e-7 * np.trace(gram, axis1=1, axis2=2)
        delta = np.linalg.solve(
            gram + damping[:, None, None] * np.eye(4), rhs[..., None]
        )[..., 0]
        delta /= np.maximum(1, np.max(np.abs(delta), axis=1) / 0.15)[:, None]
        loss = np.sum(residual**2, axis=1)
        for step in (1, 0.5, 0.25, 0.125):
            pending = np.flatnonzero(np.any(delta != 0, axis=1))
            if not len(pending):
                break
            proposed = points[pending] + step * delta[pending]
            proposed[:, 0] = np.clip(proposed[:, 0], -0.499, 0.499)
            proposed[:, 1] = np.clip(proposed[:, 1], 0.001, 0.999)
            new = model.moments(proposed, calibration.early, calibration.late)
            good = np.sum((mu[pending] - new) ** 2, axis=1) < loss[pending]
            accepted = pending[good]
            points[accepted] = proposed[good]
            delta[accepted] = 0
    pred, jac = model.moments(points, calibration.early, calibration.late, True)
    return points, mu - pred, jac


def tangent_project(values, jac):
    # SVD handles nearly singular geometries without an arbitrary inverse ridge.
    u, s, _ = np.linalg.svd(jac, full_matrices=False)
    u = u * (s > s[:, :1] * 1e-10)[:, None, :]
    return values - np.einsum("nki,ni->nk", u, np.einsum("nki,nk->ni", u, values))


def anchor_theta_jacobian(model, points, early, late, free_q=False):
    """Core-anchor differential in the calibration gauge; optional analytic q column."""
    anchors = np.vstack((early, late))
    full = np.zeros((len(points), 2 + len(anchors), 4 * len(anchors)))
    for i, anchor in enumerate(anchors):
        if i < len(early):
            _, derivative, _ = model.interval(anchor, points, True)
        else:
            _, _, derivative = model.interval(points, anchor, True)
        full[:, i + 2, 4 * i : 4 * i + 4] = derivative
    full = full[..., np.arange(4 * len(anchors)) != 3]
    if free_q:
        all_anchors = np.vstack((PAST, FUTURE, anchors))
        sign = np.r_[1, -1, np.ones(len(early)), -np.ones(len(late))]
        length = np.maximum(sign * (points[:, :1] - all_anchors[:, 0]), 0)

        def alpha_q(val):
            return _alpha_q(model, val)

        def m_q(val):
            return _m_q(model, val)

        dm = sign * (model.m(points[:, :1]) - model.m(all_anchors[:, 0]))
        safe = np.where(length > 0, dm, 1)
        radius = np.sum((points[:, None, 2:] - all_anchors[:, 2:]) ** 2, axis=-1)
        kap = np.where(
            length > 0,
            np.maximum(
                sign * (points[:, 1:2] - all_anchors[:, 1]) - radius / (2 * safe), 0
            ),
            0,
        )
        dq = sign * (m_q(points[:, :1]) - m_q(all_anchors[:, 0]))
        alpha = model.alpha(length)
        derivative = (
            alpha_q(length) / model.gamma - alpha * alpha_q(1.0) / model.gamma**2
        ) * kap**2 + alpha / model.gamma * kap * radius * dq / safe**2
        full = np.concatenate((full, derivative[..., None]), axis=-1)
    return full


def neyman_orthogonalize(points, residual, jac, calibration, free_q=False):
    """Exploratory weighted nuisance projection, returning residual and step norm."""
    u, singular, _ = np.linalg.svd(jac, full_matrices=False)
    u *= (singular > singular[:, :1] * 1e-10)[:, None, :]
    r = residual - np.einsum("nki,ni->nk", u, np.einsum("nki,nk->ni", u, residual))
    nuisance = anchor_theta_jacobian(
        BackgroundModel(calibration.q),
        points,
        calibration.early,
        calibration.late,
        free_q,
    )
    nuisance -= u @ (u.swapaxes(1, 2) @ nuisance)
    weights = np.repeat(bump_window(points), residual.shape[1])
    design = nuisance.reshape(-1, nuisance.shape[-1])
    sqrt_w = np.sqrt(weights)
    u_n, s_n, vt_n = np.linalg.svd(sqrt_w[:, None] * design, full_matrices=False)
    if len(s_n) and s_n[0] > 0:
        keep = s_n > s_n[0] * 1e-7
        proj = u_n[:, keep].T @ (sqrt_w * r.ravel())
        step = vt_n[keep].T @ (-proj / s_n[keep])
    else:
        step = np.zeros(design.shape[1])
    return r + nuisance @ step, float(np.linalg.norm(step))


def bump_window(points):
    u, v = points[:, 0], points[:, 1]
    radius = np.sum(points[:, 2:] ** 2, axis=1)
    arguments = np.column_stack(
        (1 - (u / 0.35) ** 2, 1 - ((v - 0.5) / 0.35) ** 2, 1 - radius / 0.04)
    )
    valid = np.all(arguments > 0, axis=1)
    w = np.zeros(len(points))
    w[valid] = np.exp(np.sum(1 - 1 / arguments[valid], axis=1))
    return w


def local_pair_correlator(u, weights_z, frequencies, dc):
    """Exploratory diagonal-free pair energy; finite dc uses a fine 1D FFT grid.

    weights_z already contains the complex symbol residual and bump weight.
    Frequencies must be equally spaced for the finite-dc chirp-z evaluation.
    Infinite dc uses the exact unbinned identity, including signed frequencies.
    """
    u, z, frequencies = (
        np.asarray(u, float),
        np.asarray(weights_z, complex),
        np.asarray(frequencies, float),
    )
    if (
        u.ndim != 1
        or z.shape != u.shape
        or not len(u)
        or frequencies.ndim != 1
        or not len(frequencies)
        or not np.all(np.isfinite(u))
        or not np.all(np.isfinite(z))
        or not np.all(np.isfinite(frequencies))
        or np.isnan(dc)
        or dc <= 0
    ):
        raise ValueError("finite matching vectors and positive dc required")
    diagonal = float(np.sum(abs(z) ** 2))
    if np.isinf(dc):
        return abs(z @ np.exp(-1j * np.outer(u, frequencies))) ** 2 - diagonal
    step = frequencies[1] - frequencies[0] if len(frequencies) > 1 else 0.0
    if not np.allclose(np.diff(frequencies), step, atol=1e-12, rtol=1e-12):
        raise ValueError("finite dc requires equally spaced frequencies")
    span = float(np.ptp(u))
    if span == 0:
        return np.full(len(frequencies), abs(z.sum()) ** 2 - diagonal)
    # Linear assignment makes the binning error quadratic in the grid spacing.
    spacing = min(
        span / 65536, dc / 256, 1 / (max(float(np.max(abs(frequencies))), 1) * 256)
    )
    bins = int(np.ceil(span / spacing)) + 2
    if bins > 2**22:
        raise ValueError("dc/frequency resolution exceeds finite-grid memory limit")
    position = (u - u.min()) / spacing
    index = position.astype(int)
    fraction = position - index
    grid = np.zeros(bins, complex)
    np.add.at(grid, index, z * (1 - fraction))
    np.add.at(grid, index + 1, z * fraction)
    correlation = fftconvolve(grid, grid[::-1].conj())
    lags = np.arange(-bins + 1, bins) * spacing
    # Remove each point's entire gridded self-correlation, not distinct bin mates.
    cross = float(np.sum(abs(z) ** 2 * fraction * (1 - fraction)))
    correlation[bins - 1] -= diagonal - 2 * cross
    correlation[bins - 2] -= cross
    correlation[bins] -= cross
    weighted = correlation * np.exp(-0.5 * (lags / dc) ** 2)
    transformed = czt(
        weighted,
        m=len(frequencies),
        w=np.exp(-1j * step * spacing),
        a=np.exp(1j * frequencies[0] * spacing),
    )
    return (transformed * np.exp(1j * frequencies * (bins - 1) * spacing)).real


def spectrum(points, residual, jac, calibration, frequencies):
    """Globally normalized matched filter, scanning both transverse chiralities.

    For num_i = <s_perp_i,r_perp_i>, S(w) = sum(W_i num_i exp(-iwu_i)) /
    sum(W_i |s_perp_i|^2). No individual near-zero symbol is inverted.
    Causal order is invariant under y2 -> -y2: this conjugates the circular
    polarization symbol and exchanges +omega and -omega (up to conjugation).
    The O(2) gauge, including _gauge's reflection, thus requires both hands;
    max(|S(+omega)|, |S(-omega)|) is invariant under that choice.
    """
    frequencies = np.asarray(frequencies, dtype=float)
    if (
        frequencies.ndim != 1
        or not len(frequencies)
        or not np.all(np.isfinite(frequencies))
    ):
        raise ValueError("frequencies must be a nonempty finite vector")
    model = BackgroundModel(calibration.q)
    shear = model.shear_response(points, calibration.early, calibration.late)
    perpendicular = tangent_project(shear, jac)
    residual = tangent_project(residual, jac)
    weights = bump_window(points)
    if weights.sum() <= 1e-12:
        raise ValueError("no active window mass")
    num = np.sum(np.conj(perpendicular) * residual, axis=1)
    norm_sq = np.sum(np.abs(perpendicular) ** 2, axis=1)
    norm_sum = float(np.sum(weights * norm_sq))
    denom = max(norm_sum, 1e-30)
    carrier = np.exp(-1j * np.outer(points[:, 0], frequencies))
    score = (weights * num) @ carrier / denom
    score_neg = (weights * num) @ carrier.conj() / denom
    amplitude_pos, amplitude_neg = np.abs(score), np.abs(score_neg)
    amplitude = np.maximum(amplitude_pos, amplitude_neg)
    z = weights * num / denom
    var_shot = float(np.sum(np.abs(z) ** 2))
    pwr = 0.5 * (amplitude_pos**2 + amplitude_neg**2)
    if var_shot <= 1e-28 and float(np.max(pwr)) <= 1e-28:
        c_ref_sq = 0.0
    else:
        a0 = max(var_shot, 1e-30)
        inv_w2 = 1.0 / np.maximum(frequencies**2, 1.0)
        rel_pwr = pwr / a0
        rel_hi = max(float(np.max(frequencies**2 * rel_pwr)) * 4.0, 1e6, 1e10)

        def whittle_obj(rel_c1):
            rel_mu = 1.0 + rel_c1 * inv_w2
            ratio = rel_pwr / rel_mu
            excess = np.where(
                ratio <= 4.5,
                ratio,
                4.5 + 2.0 * np.sqrt(np.maximum(ratio - 4.5 + 1e-12, 0.0)) - 2e-6,
            )
            return float(np.sum(np.log(rel_mu) + excess))

        opt = minimize_scalar(
            lambda lc: whittle_obj(float(np.exp(lc))),
            bounds=(float(np.log(1e-8)), float(np.log(rel_hi))),
            method="bounded",
            options={"xatol": 1e-12},
        )
        rel_cand = float(np.exp(opt.x))
        c_ref_sq = a0 * rel_cand if whittle_obj(rel_cand) < whittle_obj(0.0) else 0.0
    sigma_bg = np.sqrt(
        np.maximum(var_shot + c_ref_sq / np.maximum(frequencies**2, 1), 1e-30)
    )
    whitened = (
        amplitude.copy()
        if var_shot <= 1e-28 and c_ref_sq <= 1e-28
        else amplitude / sigma_bg
    )
    w_peak = int(np.argmax(whitened))
    w_off = abs(frequencies - frequencies[w_peak]) > 20
    w_floor = (
        float(np.sqrt(np.mean(whitened[w_off] ** 2))) if w_off.any() else float("nan")
    )
    peak = int(np.argmax(amplitude))
    positive = amplitude_pos[peak] >= amplitude_neg[peak]
    blocks = []
    for block in (np.arange(len(points)) % 2 == 0, np.arange(len(points)) % 2 == 1):
        mass = np.sum(weights[block] * norm_sq[block])
        phase = carrier[block] if positive else carrier[block].conj()
        blocks.append((weights[block] * num[block]) @ phase / max(mass, 1e-30))
    off = np.abs(frequencies - frequencies[peak]) > 20
    high = off & (frequencies >= 44)
    floor = float(np.sqrt(np.mean(amplitude[off] ** 2))) if off.any() else float("nan")
    floor_high = (
        float(np.sqrt(np.mean(amplitude[high] ** 2))) if high.any() else float("nan")
    )
    mean = np.mean(amplitude[off]) if off.any() else float("nan")
    sd = np.std(amplitude[off]) if off.any() else float("nan")
    return {
        "whitened_amplitude": whitened,
        "whitened_sigma_bg": sigma_bg,
        "whitened_peak_frequency": float(frequencies[w_peak]),
        "whitened_peak_signed_frequency": float(
            frequencies[w_peak]
            * (1 if amplitude_pos[w_peak] >= amplitude_neg[w_peak] else -1)
        ),
        "whitened_peak_z": float(whitened[w_peak]),
        "whitened_snr": float(whitened[w_peak] / max(w_floor, 1e-30)),
        "shot_sd": float(np.sqrt(var_shot)),
        "c_ref_estimate": float(np.sqrt(c_ref_sq)),
        "score": score,
        "score_neg": score_neg,
        "amplitude_pos": amplitude_pos,
        "amplitude_neg": amplitude_neg,
        "amplitude": amplitude,
        "peak_signed_frequency": float(frequencies[peak] * (1 if positive else -1)),
        "peak_frequency": float(frequencies[peak]),
        "peak_amplitude": float(amplitude[peak]),
        "off_peak_rms": floor,
        "off_peak_rms_high": floor_high,
        "exploratory_snr": float(amplitude[peak] / max(floor, 1e-30)),
        "exploratory_snr_high": float(amplitude[peak] / max(floor_high, 1e-30)),
        "exploratory_z": float((amplitude[peak] - mean) / max(sd, 1e-30)),
        "block_relative_difference": float(
            abs(blocks[0][peak] - blocks[1][peak]) / max(amplitude[peak], 1e-30)
        ),
        "s_perp_norm_sq_sum": norm_sum,
        "window_mass": float(weights.sum()),
        "active": int(np.count_nonzero(weights)),
        "tangent_leakage": float(
            np.max(np.abs(np.einsum("nki,nk->ni", jac, residual)))
        ),
        "label": "exploratory spectrum; no calibrated significance or p-value",
    }


def observe(
    query,
    pilot,
    references,
    probes,
    frequencies,
    q=None,
    b_min=1e-5,
    max_probes=32768,
    count_query=None,
    diagnostics=False,
    orthogonalize=True,
):
    """Full blind split protocol; all indices are arbitrary shuffled labels."""
    splits = [np.asarray(s) for s in (pilot, references, probes)]
    joined = np.concatenate(splits)
    if np.any(np.diff(np.sort(joined)) == 0):
        raise ValueError("pilot, reference and probe splits must be disjoint")
    core = 4
    early, late = select_anchors(
        query, pilot, count=core, calibration_count=10, count_query=count_query
    )
    anchor_ids = np.r_[early, late]
    counts = order_counts(query, references, anchor_ids, count_query=count_query)
    late_counts = order_counts(
        query, references, late, early=early, count_query=count_query
    )
    calibration = calibrate_progressive(
        counts, late_counts[:, 2:].T, q=q, n_ref=len(references), starts=1
    )
    ring_early, ring_late = early, late
    early, late = early[:core], late[:core]
    calibration = Calibration(
        calibration.q,
        calibration.early[:core],
        calibration.late[:core],
        calibration.rms,
        calibration.condition,
        calibration.evaluations,
    )
    pool = len(probes)
    # Anchor comparabilities cheaply discard probes whose interval counts are zero.
    comparable = np.zeros(pool, dtype=bool)
    for offset in range(0, pool, 65536):
        sub = probes[offset : offset + 65536]
        comparable[offset : offset + len(sub)] = np.all(
            query(early, sub)[0], axis=0
        ) & np.all(query(sub, late)[0], axis=1)
    probes = np.asarray(probes)[comparable][:max_probes]
    mu = order_counts(query, references, probes, early, late, count_query=count_query)
    selected = np.all(mu > b_min, axis=1)
    if np.count_nonzero(selected) < 8:
        raise ValueError("fewer than eight probes pass the order-only active filter")
    points, residual, jac = triangulate(mu[selected], calibration)
    dtheta_norm = 0.0
    if orthogonalize:
        residual, dtheta_norm = neyman_orthogonalize(
            points, residual, jac, calibration, free_q=q is None
        )
    result = spectrum(points, residual, jac, calibration, frequencies)
    result.update(
        dtheta_norm=dtheta_norm,
        orthogonalize=orthogonalize,
        q_hat=calibration.q,
        calibration_rms=calibration.rms,
        calibration_condition=calibration.condition,
        calibration_evaluations=calibration.evaluations,
        accepted=int(selected.sum()),
        n_ref=len(references),
        n_probe=len(probes),
        probe_pool=pool,
        comparable=int(comparable.sum()),
        comparable_fraction=float(comparable.mean()) if pool else 0.0,
        calibration_early_ids=ring_early.tolist(),
        calibration_late_ids=ring_late.tolist(),
        early_ids=early.tolist(),
        late_ids=late.tolist(),
    )
    if diagnostics:
        result["diagnostics"] = {
            "probe_ids": probes[selected],
            "points": points,
            "moments": mu[selected],
            "residual": residual,
            "jac": jac,
            "calibration": calibration,
        }
    return result
