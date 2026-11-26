# ruff: noqa: E501, UP031
"""Check 02: canonical slice map, eqs. (6)-(7) of gram_fourier_noise.md, and
the transverse / longitudinal shift bounds, including the tips.

  (6)  |det T - 1| <= C e(s)^2,  e(s) = lam min(1, lam/s)
  (7)  || nu_lam - mu_0 ||_TV <= C lam^4 log(1/lam)
       |delta y| <= C lam min(sqrt(s), lam/sqrt(s)),  L1 average O(lam^2)
       |delta v| <= C lam^2 uniformly

The constants C are not given by the author.  The a priori constant used
for (6) comes from |J(s,t)| <= min(l f_max, lam (2 f_max + Var f)) and
m >= l:  eps <= CJ lam min(1, lam/l),  CJ = 2 sqrt(q) (4 f_max - 2),
|det T - 1| <= 2 CJ^2 e^2 (crude).  For the other quantities the test is
that the normalised supremum stays bounded as lam -> 0 over a factor 32
in lam (threshold: within a factor 1.5 of the r = 32 value).  These
thresholds test boundedness, not a predicted constant, and were fixed
after one exploratory run of this script; see report.md.

Deterministic.  ASCII only.  Exit code 1 on any FAIL.
"""

from __future__ import annotations

import math

import numpy as np
from wavegeom import Reporter, Wave, gauss01

rep = Reporter()
Q = 1.0
RS = (32, 64, 128, 256, 512, 1024)
PHIS = (0.37, 2.1)


def sym(a, b):
    """Real 2x2 arrays for the pairs (a, b)."""
    b = np.asarray(b)
    out = np.empty(b.shape + (2, 2))
    out[..., 0, 0] = a + b.real
    out[..., 1, 1] = a - b.real
    out[..., 0, 1] = out[..., 1, 0] = b.imag
    return out


def u_grid(lam):
    """Points resolving the oscillation everywhere and both tips down to 1e-7."""
    tip = 10.0 ** np.linspace(-7.0, math.log10(0.5), 3000)
    uniform = np.arange(-0.5 + lam / 40, 0.5, lam / 40)
    return np.unique(np.concatenate((-0.5 + tip, 0.5 - tip, uniform)))


def graded_quadrature(lam):
    """Composite Gauss rule on [-1/2, 1/2], graded at the tips, panels <= lam/4."""
    edges_s = [0.0]
    s = 1e-9
    while s < 0.5:
        edges_s.append(s)
        s = min(0.5, s + min(max(s * 0.25, 1e-9), lam / 4))
    edges_s.append(0.5)
    edges_s = np.unique(np.array(edges_s))
    edges = np.unique(np.concatenate((-0.5 + edges_s, 0.5 - edges_s)))
    x, w = gauss01(6)
    nodes = (edges[:-1, None] + np.diff(edges)[:, None] * x).ravel()
    wts = (np.diff(edges)[:, None] * w).ravel()
    return nodes, wts


rows = []
for r in RS:
    lam = 1.0 / r
    best = dict(det=0.0, dy=0.0, dv=0.0, ident=0.0, tv=0.0, l1=0.0)
    for phi in PHIS:
        wave = Wave(Q, r, phi)
        u = u_grid(lam)
        sl = wave.slice(u)
        s = np.minimum(u + 0.5, 0.5 - u)
        e = lam * np.minimum(1.0, lam / s)
        best["det"] = max(best["det"], float(np.max(np.abs(sl["detT"] - 1.0) / e**2)))
        tmat = sym(sl["c0"], sl["tb"])
        left = sym(sl["aL"], sl["eL"])
        right = sym(sl["aR"], sl["eR"])
        eye = np.eye(2)
        # identity T (M_L^-1 + M_R^-1) T = I / W0
        ident = tmat @ (left + right) @ tmat * sl["W0"][:, None, None] - eye
        best["ident"] = max(best["ident"], float(np.max(np.abs(ident))))
        # delta v = (1/2) z^T [T M_L^-1 T - I/m] z, sup over |z|^2 <= 2 W0
        dvm = tmat @ left @ tmat - eye / sl["m"][:, None, None]
        dv = sl["W0"] * np.max(np.abs(np.linalg.eigvalsh(dvm)), axis=-1)
        # same from the right expression (must agree by the identity)
        dvr = -(tmat @ right @ tmat - eye / sl["mr"][:, None, None])
        dv_r = sl["W0"] * np.max(np.abs(np.linalg.eigvalsh(dvr)), axis=-1)
        # use the expression with the long denominator (better conditioned)
        dv_use = np.where(u < 0, dv_r, dv)
        best["dv"] = max(best["dv"], float(np.max(dv_use) / lam**2))
        # delta y = (T - I) z, sup over |z| <= sqrt(2 W0)
        top = np.max(np.abs(np.linalg.eigvalsh(tmat - eye)), axis=-1)
        dy = top * np.sqrt(2 * sl["W0"])
        best["dy"] = max(
            best["dy"],
            float(np.max(dy / (lam * np.minimum(np.sqrt(s), lam / np.sqrt(s))))),
        )
        # total variation of the canonical one-point laws and L1 mean of delta y
        nodes, wts = graded_quadrature(lam)
        q = wave.slice(nodes)
        base = wave.a2(nodes) * q["W0"]
        z0 = np.sum(wts * base)
        zl = np.sum(wts * base * q["detT"])
        tv = 0.5 * np.sum(wts * np.abs(base * q["detT"] / zl - base / z0))
        best["tv"] = max(best["tv"], float(tv / (lam**4 * math.log(r))))
        qt = sym(q["c0"], q["tb"])
        qtop = np.max(np.abs(np.linalg.eigvalsh(qt - eye)), axis=-1)
        # E|z| on a slice <= sqrt(2 W0); gives an upper bound for the L1 mean
        l1 = np.sum(wts * base * qtop * np.sqrt(2 * q["W0"])) / z0
        best["l1"] = max(best["l1"], float(l1 / lam**2))
    rows.append((r, best))

k = math.sqrt(Q)
fmax = 1.0 / math.cos(k / 2) ** 2
cj = 2 * k * (4 * fmax - 2)
print("a priori constant for (6): 2 CJ^2 = %.2f" % (2 * cj * cj))
print(
    "     r  sup|detT-1|/e^2  TV/(lam^4 log r)  sup|dv|/lam^2  sup|dy|/bound  L1(dy)/lam^2  identity_err"
)
for r, b in rows:
    print(
        "%6d %15.4f %17.4f %14.4f %14.4f %13.4f %13.2e"
        % (r, b["det"], b["tv"], b["dv"], b["dy"], b["l1"], b["ident"])
    )

first = rows[0][1]
rep.check(
    "identity T (M_L^-1 + M_R^-1) T = I / W0 on every slice",
    max(b["ident"] for _, b in rows) < 1e-9,
)
rep.check(
    "(6) |det T - 1| <= 2 CJ^2 e(s)^2 for all s in [1e-7, 1/2], both tips, all lam",
    max(b["det"] for _, b in rows) <= 2 * cj * cj,
    "max ratio {:.4f}".format(max(b["det"] for _, b in rows)),
)
for key, label in (
    ("det", "(6) normalised sup does not grow as lam -> 0"),
    ("dv", "longitudinal shift sup |dv| / lam^2 bounded, including tips"),
    ("dy", "transverse shift / (lam min(sqrt s, lam/sqrt s)) bounded, including tips"),
    ("l1", "L1 mean of transverse shift / lam^2 bounded"),
):
    vals = [b[key] for _, b in rows]
    rep.check(
        label,
        max(vals) <= 1.5 * first[key],
        "values " + " ".join(f"{v:.4f}" for v in vals),
    )

# (7): TV / lam^4 must be affine in log(1/lam): constant increment per
# doubling of r (rate lam^4 log(1/lam), no faster growth, and the log is real).
tv4 = [b["tv"] * math.log(r) for r, b in rows]
steps = [tv4[i + 1] - tv4[i] for i in range(len(tv4) - 1)]
mean_step = sum(steps) / len(steps)
print("TV / lam^4 per r: " + " ".join(f"{v:.2f}" for v in tv4))
print("increments per doubling of r: " + " ".join(f"{v:.2f}" for v in steps))
rep.check(
    "(7) TV / lam^4 is affine in log(1/lam): increments per doubling within 20% of their mean",
    all(0.8 * mean_step <= v <= 1.2 * mean_step for v in steps),
    f"mean increment {mean_step:.2f}, i.e. TV ~ {mean_step / math.log(2):.1f} lam^4 log(1/lam)",
)
rep.check(
    "the logarithm in (7) is genuine (increment per doubling > 0)", min(steps) > 1.0
)
rep.finish()
