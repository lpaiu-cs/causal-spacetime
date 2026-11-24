# ruff: noqa: E501, UP031
"""Check 03: kernel mismatch probability, eq. (8) of gram_fourier_noise.md.

Coupling under test: a wave reference and a background reference share
the same canonical coordinates (u, v_tilde, z); the probe X and the
anchors are FIXED canonical points.  For each raw indicator kernel we
compute Pr{K_lam != K_0} and compare with lam^2.

  past      1{R < X}        deterministic quadrature (exact in v_tilde)
  future    1{X < R}        deterministic quadrature
  interval  1{A < R < X}    deterministic quadrature, frozen anchors
  pair      1{R ~ R'}       Monte Carlo, exact in the later point's v_tilde
  ancestor  1{R<X}1{R'<X}1{R~R'}  Monte Carlo part + past-kernel bound

The author gives no constant.  The test is boundedness of Pr / lam^2 over
a factor 16 in lam, uniformly over probes at the corners of the active
set.  Pr / lam^2 depends on the relative phases r (u_X + 1/2) and r, so
each scale S in {32, 128, 512} uses a cluster r = S * {1, 1.07, ...} and
the envelope (max over cluster, probes, anchors) is compared across
scales: PASS if the envelope at the smallest wavelength (S = 512) is at
most 1.25 times the larger of the envelopes at S = 32 and S = 128.  The
cluster design and the factor 1.25 were fixed after exploratory runs
showed the phase dependence (an exploratory interval-kernel run up to
S = 2048 gave envelopes 1.52, 1.87, 1.74, 1.86); see report.md.
Also checks the L^p step (9) on the exact fourth moment of a centred
Bernoulli(lam^2) average.

Deterministic (fixed seed).  ASCII only.  Exit code 1 on any FAIL.
"""

from __future__ import annotations

import math

import numpy as np
from wavegeom import Reporter, Wave, gauss01, sqf

rep = Reporter()
Q, PHI = 1.0, 0.37
SCALES = (32, 128, 512)
MULTS = (1.0, 1.07, 1.19, 1.31, 1.43, 1.57)
MC_MULTS = (1.0, 1.31)
NMC = 200000
GX4, GW4 = gauss01(4)
RHO_IN, W_IN = gauss01(16)
RHO_OUT, W_OUT = gauss01(6)
NANG = 32
ANG = 2 * np.pi * (np.arange(NANG) + 0.5) / NANG
UNIT = np.exp(1j * ANG)

PROBES = [  # canonical (u, v_tilde, |z|, angle) on the edge of the active set
    (0.20, 0.70, 0.10, 0.3),
    (0.10, 0.65, 0.20, 1.1),
    (0.30, 0.75, 0.20, 2.5),
    (0.10, 0.75, 0.05, 4.0),
    (0.30, 0.65, 0.05, 5.2),
]
ANCHORS = [
    (-0.35, 0.15, 0.08 + 0.0j),
    (-0.25, 0.25, 0.08j),
    (-0.15, 0.35, 0.08 + 0.08j),
]


def graded_edges(lo, hi, lam):
    """Panel edges on [lo, hi], geometric towards both ends, <= lam/2 inside."""
    length = hi - lo
    s, out = 1e-8, [0.0]
    while s < length / 2:
        out.append(s)
        s += min(max(0.3 * s, 1e-8), lam / 2)
    half = np.array(out)
    return np.unique(np.concatenate((lo + half, hi - half, [lo + length / 2])))


def u_rule(lo, hi, lam):
    edges = graded_edges(lo, hi, lam)
    nodes = (edges[:-1, None] + np.diff(edges)[:, None] * GX4).ravel()
    wts = (np.diff(edges)[:, None] * GW4).ravel()
    return nodes, wts


def physical(wave, u, vt, zeta):
    sl = wave.slice(np.asarray(u, dtype=float))
    v, y = wave.to_physical(sl, vt, zeta)
    return sl, v, y


def clip(x, lo, hi):
    return np.minimum(np.maximum(x, lo), hi)


def disc_points(center, rho_max):
    """Polar rule around center covering rho in [0, 1.25 rho_max]."""
    rho = np.concatenate(
        (rho_max[:, None] * RHO_IN, rho_max[:, None] * (1 + 0.25 * RHO_OUT)), axis=1
    )
    wr = np.concatenate(
        (rho_max[:, None] * W_IN, 0.25 * rho_max[:, None] * W_OUT), axis=1
    )
    z = center[:, None, None] + np.sqrt(rho)[:, :, None] * UNIT[None, None, :]
    wz = 0.5 * wr[:, :, None] * (2 * np.pi / NANG) * np.ones(NANG)
    return z, wz


def thresholds(wave, sl, u, z, other, earlier):
    """v_tilde thresholds of the relation between R=(u,.,z) and a fixed point.

    other = (u_o, vt_o, zeta_o, v_ph_o, y_o).  If earlier is True the fixed
    point precedes R (relation other < R, lower threshold); else R < other
    (upper threshold).  Returns (theta_wave, theta_background).
    """
    uo, vto, zo, vpo, yo = other
    m = sl["m"][:, None, None]
    y = sl["c0"][:, None, None] * z + sl["tb"][:, None, None] * np.conj(z)
    dv = 0.5 * sqf(sl["aL"][:, None, None], sl["eL"][:, None, None], y) - np.abs(
        z
    ) ** 2 / (2 * m)
    if earlier:
        a, b = wave.inv_action(uo * np.ones_like(u), u)
        dm = wave.dm(uo, u)
        act = 0.5 * sqf(a[:, None, None], b[:, None, None], y - yo)
        return vpo + act - dv, vto + np.abs(z - zo) ** 2 / (2 * dm[:, None, None])
    a, b = wave.inv_action(u, uo * np.ones_like(u))
    dm = wave.dm(u, uo)
    act = 0.5 * sqf(a[:, None, None], b[:, None, None], yo - y)
    return vpo - act - dv, vto - np.abs(z - zo) ** 2 / (2 * dm[:, None, None])


def fiber(sl, z):
    lo = np.abs(z) ** 2 / (2 * sl["m"][:, None, None])
    hi = 1.0 - np.abs(z) ** 2 / (2 * sl["mr"][:, None, None])
    return lo, hi


def one_sided(wave, probe, past, zbg):
    """Pr{1{R<X} mismatch} (past) or Pr{1{X<R} mismatch} (future)."""
    ux, vx, zx, vpx, yx = probe
    lam = wave.lam
    if past:
        u, wu = u_rule(-0.5, ux, lam)
    else:
        u, wu = u_rule(ux, 0.5, lam)
    sl = wave.slice(u)
    mx = float(wave.m(ux))
    if past:
        eta_x = vx - abs(zx) ** 2 / (2 * mx)
        center = zx * sl["m"] / mx
        rho_max = 2 * sl["m"] * wave.dm(u, ux) / mx * eta_x
    else:
        mrx = wave.M - mx
        zeta_x = 1 - vx - abs(zx) ** 2 / (2 * mrx)
        center = zx * sl["mr"] / mrx
        rho_max = 2 * sl["mr"] * wave.dm(ux, u) / mrx * zeta_x
    z, wz = disc_points(center, rho_max)
    inside = np.abs(z) ** 2 <= 2 * sl["W0"][:, None, None]
    tw, tb = thresholds(wave, sl, u, z, probe, earlier=not past)
    lo, hi = fiber(sl, z)
    length = np.abs(clip(tw, lo, hi) - clip(tb, lo, hi)) * inside
    per_u = np.sum(length * wz, axis=(1, 2)) * wave.a2(u)
    s = np.minimum(u + 0.5, 0.5 - u)
    total = np.sum(wu * per_u) / zbg
    tip = np.sum((wu * per_u)[s < 4 * lam]) / zbg
    return total, tip


def interval_kernel(wave, anchor, probe, zbg):
    """Pr{1{A<R<X} mismatch} with A and X fixed."""
    ua, va, za, vpa, ya = anchor
    ux, vx, zx, vpx, yx = probe
    u, wu = u_rule(ua, ux, wave.lam)
    sl = wave.slice(u)
    d1, d2, d = wave.dm(ua, u), wave.dm(u, ux), float(wave.dm(ua, ux))
    eta = vx - va - abs(zx - za) ** 2 / (2 * d)
    center = za + (zx - za) * d1 / d
    rho_max = 2 * d1 * d2 / d * eta
    z, wz = disc_points(center, rho_max)
    lo, hi = fiber(sl, z)
    aw, ab = thresholds(wave, sl, u, z, anchor, earlier=True)
    xw, xb = thresholds(wave, sl, u, z, probe, earlier=False)
    aw, ab, xw, xb = (clip(t, lo, hi) for t in (aw, ab, xw, xb))
    len_w = np.maximum(xw - aw, 0.0)
    len_b = np.maximum(xb - ab, 0.0)
    both = np.maximum(np.minimum(xw, xb) - np.maximum(aw, ab), 0.0)
    length = len_w + len_b - 2 * both
    return np.sum(wu * np.sum(length * wz, axis=(1, 2)) * wave.a2(u)) / zbg


def sample_law(wave, rng, n):
    """iid canonical/background points: (u, rho=|z|^2, angle, uniform in fiber)."""
    out = []
    need = n
    bound = wave.M / 4
    while need > 0:
        u = rng.uniform(-0.5, 0.5, size=2 * need + 16)
        m = wave.m(u)
        keep = rng.uniform(size=u.size) * bound < wave.a2(u) * m * (wave.M - m) / wave.M
        out.append(u[keep][:need])
        need -= out[-1].size
    u = np.concatenate(out)
    m = wave.m(u)
    w0 = m * (wave.M - m) / wave.M
    rho = 2 * w0 * (1 - np.sqrt(1 - rng.uniform(size=n)))
    z = np.sqrt(rho) * np.exp(2j * np.pi * rng.uniform(size=n))
    frac = rng.uniform(size=n)
    return u, z, frac


def pair_and_ancestor(wave, probe, rng, n):
    """Monte Carlo for 1{R~R'} and the X-ancestor pair kernel."""
    ux, vx, zx, vpx, yx = probe
    u1, z1, f1 = sample_law(wave, rng, n)
    u2, z2, _ = sample_law(wave, rng, n)
    swap = u2 < u1
    u1, u2 = np.where(swap, u2, u1), np.where(swap, u1, u2)
    z1, z2 = np.where(swap, z2, z1), np.where(swap, z1, z2)
    s1, s2 = wave.slice(u1), wave.slice(u2)
    lo1 = np.abs(z1) ** 2 / (2 * s1["m"])
    hi1 = 1 - np.abs(z1) ** 2 / (2 * s1["mr"])
    vt1 = lo1 + f1 * (hi1 - lo1)
    v1, y1 = wave.to_physical(s1, vt1, z1)
    y2 = s2["c0"] * z2 + s2["tb"] * np.conj(z2)
    dv2 = 0.5 * sqf(s2["aL"], s2["eL"], y2) - np.abs(z2) ** 2 / (2 * s2["m"])
    lo2 = np.abs(z2) ** 2 / (2 * s2["m"])
    hi2 = 1 - np.abs(z2) ** 2 / (2 * s2["mr"])
    a, b = wave.inv_action(u1, u2)
    low_w = v1 + 0.5 * sqf(a, b, y2 - y1) - dv2  # R < R' iff vt2 >= low
    low_b = vt1 + np.abs(z2 - z1) ** 2 / (2 * wave.dm(u1, u2))
    cw, cb = clip(low_w, lo2, hi2), clip(low_b, lo2, hi2)
    width = hi2 - lo2
    pair = np.abs(cw - cb) / width
    # ancestor-pair kernel: needs R < X, R' < X, R ~ R'
    ok = u2 < ux
    u2s = np.where(ok, u2, ux - 0.01)  # dummy value where not ok
    a, b = wave.inv_action(u2s, ux * np.ones_like(u2s))
    up_w = vpx - 0.5 * sqf(a, b, yx - y2) - dv2  # R' < X iff vt2 <= up
    up_b = vx - np.abs(zx - z2) ** 2 / (2 * wave.dm(u2s, ux))
    uw, ub = clip(up_w, lo2, hi2), clip(up_b, lo2, hi2)
    u1s = np.where(ok, u1, ux - 0.01)  # ok implies u1 < u2 < ux
    a, b = wave.inv_action(u1s, ux * np.ones_like(u1s))
    r_w = (v1 <= vpx - 0.5 * sqf(a, b, yx - y1)) & ok
    r_b = (vt1 <= vx - np.abs(zx - z1) ** 2 / (2 * wave.dm(u1s, ux))) & ok
    len_w = np.maximum(uw - cw, 0.0) * r_w
    len_b = np.maximum(ub - cb, 0.0) * r_b
    both = np.maximum(np.minimum(uw, ub) - np.maximum(cw, cb), 0.0) * (r_w & r_b)
    anc = (len_w + len_b - 2 * both) / width
    return pair.mean(), pair.std() / math.sqrt(n), anc.mean(), anc.std() / math.sqrt(n)


keys = ("past", "future", "interval", "pair", "ancestor", "tip")
env = {key: [] for key in keys}
print(
    "       r   past/lam^2  future/lam^2  interval/lam^2  tip-part/lam^2   pair/lam^2 (se)   ancestor/lam^2 (se)"
)
for scale in SCALES:
    best = {key: 0.0 for key in keys}
    for mult in MULTS:
        r = scale * mult
        wave = Wave(Q, r, PHI)
        lam = wave.lam
        nodes, wts = u_rule(-0.5, 0.5, 1 / 64)
        mm = wave.m(nodes)
        zbg = math.pi * np.sum(wts * wave.a2(nodes) * mm * (wave.M - mm) / wave.M)
        anchors = []
        for ua, va, za in ANCHORS:
            _, vpa, ya = physical(wave, ua, va, za)
            anchors.append((ua, va, za, float(vpa), complex(ya)))
        row = dict(past=0.0, future=0.0, interval=0.0, tip=0.0)
        probes = []
        for ux, vx, rad, ang in PROBES:
            zx = rad * np.exp(1j * ang)
            _, vpx, yx = physical(wave, ux, vx, zx)
            probe = (ux, vx, zx, float(vpx), complex(yx))
            probes.append(probe)
            p, tip_p = one_sided(wave, probe, True, zbg)
            f, tip_f = one_sided(wave, probe, False, zbg)
            row["past"] = max(row["past"], p / lam**2)
            row["future"] = max(row["future"], f / lam**2)
            row["tip"] = max(row["tip"], (tip_p + tip_f) / lam**2)
            for anchor in anchors:
                row["interval"] = max(
                    row["interval"], interval_kernel(wave, anchor, probe, zbg) / lam**2
                )
        text = "{:8.2f} {:12.4f} {:13.4f} {:15.4f} {:15.5f}".format(
            r, row["past"], row["future"], row["interval"], row["tip"]
        )
        if mult in MC_MULTS:
            rng = np.random.default_rng(1000 + int(round(r)))
            pm, pse, am, ase = pair_and_ancestor(wave, probes[0], rng, NMC)
            row["pair"], row["ancestor"] = pm / lam**2, am / lam**2
            text += f" {pm / lam**2:12.4f} ({pse / lam**2:.4f}) {am / lam**2:13.4f} ({ase / lam**2:.4f})"
        print(text, flush=True)
        for key, value in row.items():
            best[key] = max(best[key], value)
    for key in keys:
        env[key].append(best[key])

print(f"envelopes at scales {SCALES}:")
for key in keys:
    print("  %-9s %s" % (key, " ".join(f"{v:.5f}" for v in env[key])))
for key, label in (
    ("past", "(8) past kernel 1{R<X}: envelope of Pr/lam^2 bounded across scales"),
    ("future", "(8) future kernel 1{X<R}: envelope of Pr/lam^2 bounded across scales"),
    (
        "interval",
        "(8) anchor-interval kernel 1{A<R<X}, anchors fixed: envelope bounded",
    ),
    ("pair", "(8) reference-pair comparability 1{R~R'}: Pr/lam^2 bounded"),
    ("ancestor", "(8) ancestor-pair kernel at fixed X: Pr/lam^2 bounded"),
):
    vals = env[key]
    rep.check(
        label,
        vals[-1] <= 1.25 * max(vals[:-1]),
        "envelopes " + " ".join(f"{v:.4f}" for v in vals),
    )
rep.check(
    "tip region (s < 4 lam) contributes o(lam^2) to the one-sided kernels",
    env["tip"][-1] < 0.25 * env["tip"][0] and env["tip"][-1] < 0.05 * env["past"][-1],
    "envelopes " + " ".join(f"{v:.5f}" for v in env["tip"]),
)

# (9): exact L^4 norm of a centred Bernoulli(pi) mean, pi = c lam^2.
ok_in, ok_out = True, False
for lam in (1 / 32, 1 / 128, 1 / 512):
    pi = lam**2
    for nl2 in (1, 4, 100, 1e4):
        n = nl2 / lam**2
        mu4 = n * pi * (1 - pi) * (1 + 3 * (n - 2) * pi * (1 - pi))
        norm = mu4**0.25 / n
        ok_in &= norm <= 2**0.5 * lam / math.sqrt(n)
    n = 1e-4 / lam**2
    mu4 = n * pi * (1 - pi) * (1 + 3 * (n - 2) * pi * (1 - pi))
    ok_out |= mu4**0.25 / n > 5 * lam / math.sqrt(n)
rep.check(
    "(9) L^4 norm of centred mismatch mean <= sqrt(2) lam/sqrt(n) whenever n lam^2 >= 1",
    ok_in,
)
rep.check(
    "(9) the restriction n lam^2 >= 1 is needed (bound fails at n lam^2 = 1e-4)", ok_out
)
rep.finish()
