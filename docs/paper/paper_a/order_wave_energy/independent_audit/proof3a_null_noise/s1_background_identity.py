# ruff: noqa: E501, UP031
"""Referee checks of Sections 1 and 2 of gram_fourier_noise.md in the TRUE
3+1 effective background  ds^2 = -2 du dv + a(u)^2 dy^2,  a = cos(sqrt(q) u).

Independent re-implementation (numpy only) of the calibration defined in
null_coordinate_challenge.md / gram_fourier_challenge.md:
    p = alpha(s) eta^2/gamma,  f = alpha(1-s) zeta^2/gamma,  c = a/p^2 = p0(q s^2),
    b_AX = alpha(u_X-u_A) slack_AX^2/gamma,
    g_AX = (R_A^2+R_X^2)/2 - (m_X-m_A)(v_X-v_A-sqrt(b gamma/alpha(u_X-u_A))).

I1  p0 is the stated function, monotone on [0,1].
I2  the interval-volume law depends only on the u-separation (quadrature).
I3  p, f, b ARE the order probabilities (Monte Carlo sprinkling).
I4  c = p0(q s^2) IS the comparable-pair probability (Monte Carlo).
I5  g_AX = y_A.y_X and H(X) = 0 exactly, random chains and the frozen anchors.
I6  the identity fails outside the regular set (so the cutoff is needed).
I7  raw kernels are monotone in u along fixed-(v,y) lines (null generators).
I8  the frozen weight support is strictly inside the diamond.
INFO lines quantify how noisy u_hat is at finite n (constants, not a verdict).

Deterministic; prints PASS/FAIL; exits nonzero on any FAIL.
"""

from __future__ import annotations

import math
import sys

import numpy as np

FAILS = []


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL") + " " + name + (" :: " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- geometry
def hyp0f1(b, x, terms=80):
    x = np.asarray(x, dtype=float)
    total = np.ones_like(x)
    term = np.ones_like(x)
    for k in range(terms):
        term = term * x / ((b + k) * (k + 1))
        total = total + term
    return total


def c_fun(r, z):
    return hyp0f1(r + 0.5, -np.asarray(z, dtype=float) / 4) / math.factorial(2 * r - 1)


def p0(z):
    return c_fun(1, z) * c_fun(3, z) / (3 * c_fun(2, z) ** 2)


def alpha(length, q):
    x = math.sqrt(q) * np.asarray(length, dtype=float)
    return (1 - x / np.tan(x)) / (2 * q)


def mcoord(u, q):
    k = math.sqrt(q)
    return (np.tan(k * np.asarray(u, dtype=float)) + math.tan(k / 2)) / k


def precedes(first, second, q):
    """first < second for arrays of events (u, v, y1, y2) broadcast together."""
    dm = mcoord(second[..., 0], q) - mcoord(first[..., 0], q)
    dv = second[..., 1] - first[..., 1]
    dy2 = (second[..., 2] - first[..., 2]) ** 2 + (second[..., 3] - first[..., 3]) ** 2
    return (dm > 0) & (2 * dm * dv > dy2)


def slack(first, second, q):
    dm = mcoord(second[..., 0], q) - mcoord(first[..., 0], q)
    dy2 = (second[..., 2] - first[..., 2]) ** 2 + (second[..., 3] - first[..., 3]) ** 2
    return second[..., 1] - first[..., 1] - dy2 / (2 * dm)


def population(point, q):
    """Exact (p, f, a) of an event from the closed forms of the cited documents."""
    u, v = point[..., 0], point[..., 1]
    r2 = point[..., 2] ** 2 + point[..., 3] ** 2
    m, total, gamma = mcoord(u, q), mcoord(0.5, q), alpha(1.0, q)
    eta = v - r2 / (2 * m)
    zeta = 1 - v - r2 / (2 * (total - m))
    s = u + 0.5
    p = alpha(s, q) * eta**2 / gamma
    f = alpha(1 - s, q) * zeta**2 / gamma
    return p, f, p0(q * s**2) * p**2


def interval_probability(first, second, q):
    return (
        alpha(second[..., 0] - first[..., 0], q)
        * slack(first, second, q) ** 2
        / alpha(1.0, q)
    )


def state(p, f, a, q):
    """Observer-side inverse calibration: (u, m, v, R^2) from (p, f, a) and q."""
    p, f, a = (np.asarray(x, dtype=float) for x in (p, f, a))
    chain = a / p**2
    lo, hi = np.zeros_like(chain), np.ones_like(chain)
    for _ in range(70):  # p0 is decreasing: bisection on s in (0, 1)
        mid = 0.5 * (lo + hi)
        above = p0(q * mid**2) > chain
        lo = np.where(above, mid, lo)
        hi = np.where(above, hi, mid)
    s = 0.5 * (lo + hi)
    u = s - 0.5
    m, total, gamma = mcoord(u, q), mcoord(0.5, q), alpha(1.0, q)
    eta = np.sqrt(p * gamma / alpha(s, q))
    zeta = np.sqrt(f * gamma / alpha(1 - s, q))
    t = m / total
    v = (1 - t) + t * eta - (1 - t) * zeta
    r2 = 2 * m * (total - m) / total * (1 - eta - zeta)
    return np.stack(np.broadcast_arrays(u, m, v, r2), axis=-1)


def gram(first, second, b, q):
    gap = second[..., 0] - first[..., 0]
    eta0 = np.sqrt(b * alpha(1.0, q) / alpha(gap, q))
    return (first[..., 3] + second[..., 3]) / 2 - (second[..., 1] - first[..., 1]) * (
        second[..., 2] - first[..., 2] - eta0
    )


def residual(anchors, probes, q):
    """Return (max |g - y.y| over all pairs, H(X) per probe) from exact moments only."""
    st_a = [state(*population(a, q), q) for a in anchors]
    st_x = state(*population(probes, q), q)
    worst = 0.0
    gmat = np.zeros((3, 3))
    for i in range(3):
        gmat[i, i] = st_a[i][3]
        for j in range(i + 1, 3):
            b = interval_probability(anchors[i], anchors[j], q)
            gmat[i, j] = gmat[j, i] = gram(st_a[i], st_a[j], b, q)
            worst = max(worst, abs(gmat[i, j] - anchors[i][2:] @ anchors[j][2:]))
        worst = max(worst, abs(gmat[i, i] - anchors[i][2:] @ anchors[i][2:]))
    ell = np.r_[-np.linalg.solve(gmat[:2, :2], gmat[:2, 2]), 1.0]
    h = np.zeros(probes.shape[0])
    for i in range(3):
        b = interval_probability(anchors[i][None, :], probes, q)
        g = gram(st_a[i][None, :], st_x, b, q)
        worst = max(worst, float(np.max(np.abs(g - probes[:, 2:] @ anchors[i][2:]))))
        h += ell[i] * g
    return worst, h


def sprinkle(rng, count, q):
    """iid events from normalised physical volume a(u)^2 du dv d^2y in the diamond."""
    total = mcoord(0.5, q)
    radius = math.sqrt(total / 2)
    out = []
    have = 0
    while have < count:
        size = 4 * count
        u = rng.random(size) - 0.5
        v = rng.random(size)
        y = (2 * rng.random((size, 2)) - 1) * radius
        r2 = np.sum(y * y, axis=1)
        m = mcoord(u, q)
        keep = (v - r2 / (2 * m) > 0) & (1 - v - r2 / (2 * (total - m)) > 0)
        keep &= rng.random(size) < np.cos(math.sqrt(q) * u) ** 2
        out.append(np.column_stack([u[keep], v[keep], y[keep]]))
        have += int(keep.sum())
    return np.concatenate(out)[:count]


def pair_fraction(points, q):
    """Fraction of comparable pairs and its standard error (U-statistic h_1 formula)."""
    n = points.shape[0]
    row = np.zeros(n)
    for start in range(0, n, 400):
        block = points[start : start + 400, None, :]
        rel = precedes(block, points[None, :, :], q) | precedes(
            points[None, :, :], block, q
        )
        row[start : start + 400] = rel.sum(axis=1)
    frac = row / (n - 1)
    return (
        float(frac.mean()),
        float(2 * frac.std(ddof=1) / math.sqrt(n)),
        float(frac.std(ddof=1)),
    )


# ------------------------------------------------------------------ checks
def i1_p0():
    z = np.linspace(0, 1, 2001)
    val = p0(z)
    dlog = np.gradient(np.log(val), z)
    check(
        "I1 p0(0)=1/10, p0 strictly decreasing on [0,1], (log p0)' < -1/50, p0(1) > 0.09",
        abs(val[0] - 0.1) < 1e-14
        and np.all(np.diff(val) < 0)
        and dlog.max() < -0.02
        and val[-1] > 0.09,
        f"p0(1)={val[-1]:.6f}, (log p0)' in [{dlog.min():.4f}, {dlog.max():.4f}]",
    )


def i2_interval_law():
    x, w = np.polynomial.legendre.leggauss(80)
    rng = np.random.default_rng(1)
    worst = 0.0
    for q in (1.0, 0.25, 0.01):
        for _ in range(40):
            ua, ux = np.sort(rng.uniform(-0.5, 0.5, 2))
            if ux - ua < 1e-3:
                continue
            u = 0.5 * (ua + ux) + 0.5 * (ux - ua) * x
            m1, m2 = mcoord(u, q) - mcoord(ua, q), mcoord(ux, q) - mcoord(u, q)
            integral = (
                0.5
                * (ux - ua)
                * np.sum(w * np.cos(math.sqrt(q) * u) ** 2 * m1 * m2 / (m1 + m2))
            )
            worst = max(worst, abs(integral / alpha(ux - ua, q) - 1))
    check(
        "I2 int a^2 m1 m2/(m1+m2) du = alpha(u_X-u_A): interval volume depends on the u-gap only",
        worst < 1e-10,
        f"worst relative gap {worst:.1e}",
    )


def i3_i4_monte_carlo():
    q = 1.0
    rng = np.random.default_rng(2)
    pts = sprinkle(rng, 600000, q)
    n = pts.shape[0]
    tests = [
        np.array([0.20, 0.70, 0.10, 0.10]),
        np.array([-0.15, 0.35, 0.08, 0.08]),
        np.array([0.30, 0.75, 0.05, 0.00]),
    ]
    zmax = 0.0
    for x in tests:
        p, f, _ = population(x, q)
        for prob, count in (
            (p, precedes(pts, x[None, :], q).sum()),
            (f, precedes(x[None, :], pts, q).sum()),
        ):
            zmax = max(zmax, abs(count / n - prob) / math.sqrt(prob * (1 - prob) / n))
    pairs = [
        (tests[1], tests[0]),
        (tests[1], tests[2]),
        (np.array([-0.35, 0.15, 0.08, 0.0]), tests[1]),
    ]
    for first, second in pairs:
        assert slack(first, second, q) > 0
        prob = float(interval_probability(first, second, q))
        count = (
            precedes(first[None, :], pts, q) & precedes(pts, second[None, :], q)
        ).sum()
        zmax = max(zmax, abs(count / n - prob) / math.sqrt(prob * (1 - prob) / n))
    check(
        "I3 Monte Carlo sprinkling (6e5 events, q=1): p, f, b formulas are the order probabilities",
        zmax < 4.5,
        f"largest |z| over 9 probabilities = {zmax:.2f}",
    )

    frac, se, _ = pair_fraction(pts[:6000], q)
    target = float(p0(q))
    check(
        "I4a comparable-pair probability of the whole diamond equals p0(q), q=1",
        abs(frac - target) < 4.5 * se,
        f"MC {frac:.5f} +- {se:.5f}, p0(1)={target:.5f}, flat value 0.1 is {abs(frac - 0.1) / se:.1f} sigma away",
    )
    x = tests[2]
    past = pts[precedes(pts, x[None, :], q)][:6000]
    frac, se, spread = pair_fraction(past, q)
    target = float(p0(q * (x[0] + 0.5) ** 2))
    check(
        "I4b comparable-pair probability of a past subdiamond equals p0(q s^2), s=0.8",
        abs(frac - target) < 4.5 * se,
        "MC %.5f +- %.5f over %d ancestors, p0(q s^2)=%.5f"
        % (frac, se, past.shape[0], target),
    )
    # stress test of the p0 formula far from the flat value (outside the documents' q range)
    q6 = 6.0
    pts6 = sprinkle(np.random.default_rng(3), 6000, q6)
    frac, se, _ = pair_fraction(pts6, q6)
    target = float(p0(q6))
    check(
        "I4c formula stress test at q=6 (whole diamond): pair probability equals p0(6)",
        abs(frac - target) < 4.5 * se and abs(frac - 0.1) > 10 * se,
        f"MC {frac:.5f} +- {se:.5f}, p0(6)={target:.5f}, flat value 0.1 is {abs(frac - 0.1) / se:.0f} sigma away",
    )

    # INFO: how noisy is u_hat?  delta c = (2/(n p)) sum_{R<X} (pi(R) - c), so
    # sd(c_hat) = 2 sd(pi)/sqrt(n p) and sd(u_hat) = sd(c_hat)/|dc/du|.
    x = np.array(
        [0.20, 0.70, 0.10, 0.10]
    )  # centre of the frozen probe window (R^2 = 0.02)
    past = pts[precedes(pts, x[None, :], q)][:6000]
    _, _, spread = pair_fraction(past, q)
    p, _, _ = population(x, q)
    s = x[0] + 0.5
    eps = 1e-5
    dcds = float((p0(q * (s + eps) ** 2) - p0(q * (s - eps) ** 2)) / (2 * eps))
    coeff = 2 * spread / math.sqrt(float(p)) / abs(dcds)
    print(
        f"INFO constants at the window centre, q=1: p={float(p):.4f}, sd(pi)={spread:.4f}, dc/du={dcds:.5f}, sd(u_hat) ~ {coeff:.1f}/sqrt(n)"
    )
    for target_sd, label in (
        (0.05, "half the half-width of the u-window"),
        (1 / 64, "omega*sd(u_hat)=1 at omega=64"),
    ):
        print(
            f"INFO references needed for sd(u_hat)={target_sd:.4f} ({label}): n ~ {(coeff / target_sd) ** 2:.1e}"
        )


def random_chains(rng, pts, q, want):
    out = []
    while len(out) < want:
        idx = rng.integers(0, pts.shape[0], size=(400000, 4))
        quad = pts[idx]
        order = np.argsort(quad[:, :, 0], axis=1)
        quad = np.take_along_axis(quad, order[:, :, None], axis=1)
        ok = (
            precedes(quad[:, 0], quad[:, 1], q)
            & precedes(quad[:, 1], quad[:, 2], q)
            & precedes(quad[:, 2], quad[:, 3], q)
        )
        det = quad[:, 0, 2] * quad[:, 1, 3] - quad[:, 0, 3] * quad[:, 1, 2]
        ok &= np.abs(det) > 1e-3
        ok &= np.min(quad[:, 1:, 0] - quad[:, :-1, 0], axis=1) > 0.02
        out.extend(list(quad[ok]))
    return out[:want]


def frozen_anchors(theta):
    rot = np.array(
        [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
    )
    y = 0.08 * np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]) @ rot.T
    return np.column_stack([[-0.35, -0.25, -0.15], [0.15, 0.25, 0.35], y])


def window_probes():
    us = np.linspace(0.1, 0.3, 9)
    vs = np.linspace(0.65, 0.75, 5)
    r2 = np.linspace(0.0025, 0.04, 5)
    ang = np.arange(8) * math.pi / 4 + 0.1
    u, v, r, a = np.meshgrid(us, vs, np.sqrt(r2), ang, indexing="ij")
    return np.column_stack(
        [u.ravel(), v.ravel(), (r * np.cos(a)).ravel(), (r * np.sin(a)).ravel()]
    )


def i5_i6_identity():
    for q in (1.0, 0.25):
        rng = np.random.default_rng(int(100 * q) + 5)
        pts = sprinkle(rng, 60000, q)
        worst_g, worst_h, scale = 0.0, 0.0, 0.0
        chains = random_chains(rng, pts, q, 200)
        for quad in chains:
            wg, h = residual(quad[:3], quad[3:4], q)
            worst_g, worst_h = max(worst_g, wg), max(worst_h, float(np.abs(h).max()))
            scale = max(scale, float(np.abs(quad[:, 2:]).max()) ** 2)
        check(
            f"I5a q={q:g}: 200 random timelike chains A1<A2<A3<X: g_AX = y_A.y_X and H(X)=0 from exact moments",
            worst_g < 1e-9 and worst_h < 1e-8,
            f"max|g-y.y|={worst_g:.1e}  max|H|={worst_h:.1e}",
        )
    probes = window_probes()
    for q in (1.0, 0.25, 0.01):
        worst_g, worst_h, min_slack = 0.0, 0.0, 1.0
        for theta in (0.0, math.pi / 8):
            anchors = frozen_anchors(theta)
            wg, h = residual(anchors, probes, q)
            worst_g, worst_h = max(worst_g, wg), max(worst_h, float(np.abs(h).max()))
            for i in range(3):
                min_slack = min(
                    min_slack, float(slack(anchors[i][None, :], probes, q).min())
                )
                for j in range(i + 1, 3):
                    min_slack = min(min_slack, float(slack(anchors[i], anchors[j], q)))
        check(
            f"I5b q={q:g}: frozen anchors (theta=0, pi/8), 1800 probes on the closed weight support: H=0, all slacks > 0",
            worst_g < 1e-10 and worst_h < 1e-9 and min_slack > 0,
            f"max|g-y.y|={worst_g:.1e}  max|H|={worst_h:.1e}  min causal slack={min_slack:.4f}",
        )
    # outside the regular set: a probe that is NOT in the future of anchor 3
    q = 1.0
    anchors = frozen_anchors(0.0)
    probe = np.array([[0.20, 0.36, 0.0, 0.0]])
    related = [bool(precedes(a[None, :], probe, q)[0]) for a in anchors]
    st_a = [state(*population(a, q), q) for a in anchors]
    st_x = state(*population(probe, q), q)
    gmat = np.array([[a[2:] @ b[2:] for b in anchors] for a in anchors])
    ell = np.r_[-np.linalg.solve(gmat[:2, :2], gmat[:2, 2]), 1.0]
    h = 0.0
    for i in range(3):
        b = (
            interval_probability(anchors[i][None, :], probe, q)
            if related[i]
            else np.zeros(1)
        )
        h += ell[i] * float(gram(st_a[i][None, :], st_x, b, q)[0])
    check(
        "I6 non-vacuity: a probe outside the future of an anchor has H != 0 (b=0 there), so (2) needs the cutoff",
        related == [True, True, False] and abs(h) > 1e-3,
        f"anchor<probe: {related}, H = {h:.4f}",
    )


def i7_monotone():
    q = 1.0
    rng = np.random.default_rng(7)
    refs = sprinkle(rng, 1500, q)
    lines = sprinkle(rng, 60, q)[:, 1:]  # (v, y1, y2) of interior events
    us = np.linspace(-0.49, 0.49, 400)
    anchor = np.array([-0.35, 0.15, 0.08, 0.0])
    after_anchor = precedes(anchor[None, :], refs, q)
    comparable = precedes(refs[:300, None, :], refs[None, :300, :], q)
    comparable = comparable | comparable.T
    bad = {"past": 0, "future": 0, "interval": 0, "pair": 0}
    jumps = 0
    for line in lines:
        probe = np.column_stack([us, np.tile(line, (us.size, 1))])
        past = precedes(refs[None, :, :], probe[:, None, :], q).astype(np.int8)
        fut = precedes(probe[:, None, :], refs[None, :, :], q).astype(np.int8)
        bad["past"] += int(np.count_nonzero(np.diff(past, axis=0) < 0))
        bad["future"] += int(np.count_nonzero(np.diff(fut, axis=0) > 0))
        bad["interval"] += int(
            np.count_nonzero(np.diff(past * after_anchor[None, :], axis=0) < 0)
        )
        sub = past[:, :300]
        pair = sub[:, :, None] * sub[:, None, :] * comparable[None, :, :]
        bad["pair"] += int(np.count_nonzero(np.diff(pair, axis=0) < 0))
        jumps += int(np.count_nonzero(np.diff(past, axis=0)))
    check(
        "I7a along fixed-(v,y) lines: past / anchor-interval / ancestor-pair kernels never decrease in u, future never increases",
        sum(bad.values()) == 0 and jumps > 1000,
        "violations %s over 60 lines x 1500 references (%d past-kernel jumps seen)"
        % (bad, jumps),
    )
    # control: along a transverse (spacelike) line the past kernel is NOT monotone
    ys = np.linspace(-0.3, 0.3, 400)
    probe = np.column_stack(
        [np.full(ys.size, 0.2), np.full(ys.size, 0.7), ys, np.zeros(ys.size)]
    )
    past = precedes(refs[None, :, :], probe[:, None, :], q).astype(np.int8)
    d = np.diff(past, axis=0)
    both = int(np.count_nonzero((d > 0).any(axis=0) & (d < 0).any(axis=0)))
    check(
        "I7b control: along a transverse line the past kernel does switch on and off (monotonicity is special to u)",
        both > 10,
        "%d of 1500 references enter and leave the past" % both,
    )


def i8_support():
    worst = {}
    probes = window_probes()
    for q in (1.0, 0.25, 0.01):
        total = mcoord(0.5, q)
        m = mcoord(probes[:, 0], q)
        r2 = probes[:, 2] ** 2 + probes[:, 3] ** 2
        eta = probes[:, 1] - r2 / (2 * m)
        zeta = 1 - probes[:, 1] - r2 / (2 * (total - m))
        worst[q] = (float(eta.min()), float(zeta.min()))
    check(
        "I8 frozen weight support (u in [.1,.3], v in [.65,.75], R^2 in [.0025,.04]) is strictly inside the diamond",
        all(min(pair) > 0.05 for pair in worst.values()),
        "min (eta, zeta) per q: "
        + ", ".join(f"q={q:g}: ({a:.3f}, {b:.3f})" for q, (a, b) in worst.items()),
    )


if __name__ == "__main__":
    i1_p0()
    i2_interval_law()
    i3_i4_monte_carlo()
    i5_i6_identity()
    i7_monotone()
    i8_support()
    print("SUMMARY failures=%d" % len(FAILS))
    sys.exit(1 if FAILS else 0)
