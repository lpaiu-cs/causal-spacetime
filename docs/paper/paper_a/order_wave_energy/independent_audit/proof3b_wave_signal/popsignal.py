"""Referee-side population Fourier coefficient of the anchor Gram residual.

Independent re-implementation of the quantity in gram_fourier_challenge.md
(eq. (2) at population level, eq. (4) for the prediction).  Nothing is
imported from experiments/positive_control.  Weight and Fourier phase are
applied to RECONSTRUCTED coordinates (observer-side inverse calibration of
the exact wave order probabilities), exactly as the challenge specifies.

ASCII only.
"""

from __future__ import annotations

import math

import numpy as np
from wavegeom import Wave, bump, gauss01, sqf

U_LO, U_HI = 0.1, 0.3
V_LO, V_HI = 0.65, 0.75
T_LO, T_HI = 0.0025, 0.04


def frozen_anchors(theta=0.0, reflect=False, d=0.08):
    z = d * np.array([1.0 + 0.0j, 1.0j, 1.0 + 1.0j]) * np.exp(1j * theta)
    if reflect:
        z = np.conj(z)
    return [(-0.35, 0.15, z[0]), (-0.25, 0.25, z[1]), (-0.15, 0.35, z[2])]


def composite(lo, hi, panels, order):
    x, w = gauss01(order)
    edges = np.linspace(lo, hi, panels + 1)
    nodes = (edges[:-1, None] + np.diff(edges)[:, None] * x).ravel()
    wts = (np.diff(edges)[:, None] * w).ravel()
    return nodes, wts


class Probes:
    """Wave-side data on the probe u nodes, shared by all anchor choices."""

    def __init__(self, q, r, phi, exact_u=True, u_panels=None, u_order=6):
        self.wave = wave = Wave(q, r, phi)
        self.r = r
        self._ax = {}
        panels = u_panels or int(math.ceil((U_HI - U_LO) * r / 1.5))
        self.u, self.wu = composite(U_LO, U_HI, panels, u_order)
        self.sl = wave.slice(self.u)
        self.gamma_ratio = wave.chain(-0.5, 0.5)[0]
        nu = self.u.size
        self.al = np.empty(nu)
        self.ar = np.empty(nu)
        self.cr = np.ones(nu)
        for j, u in enumerate(self.u):
            self.al[j], c = wave.chain(-0.5, u, pairs=exact_u)
            if exact_u:
                self.cr[j] = c
            self.ar[j] = wave.chain(u, 0.5)[0]
        s = self.u + 0.5
        # observer-side u: invert the background pair-probability calibration
        self.s_hat = wave.invert_p0(wave.p0(s) * self.cr) if exact_u else s
        self.u_hat = self.s_hat - 0.5
        # eta_hat / eta and zeta_hat / zeta (exact wave volume hierarchy)
        self.eta_scale = np.sqrt(
            self.al * wave.alpha0(s) / (self.gamma_ratio * wave.alpha0(self.s_hat))
        )
        self.zeta_scale = np.sqrt(
            self.ar
            * wave.alpha0(1 - s)
            / (self.gamma_ratio * wave.alpha0(1 - self.s_hat))
        )

    def anchor_state(self, anchor, reconstructed=True):
        """Physical point and observer-side state of one canonical anchor."""
        wave = self.wave
        ua, va, za = anchor
        sl = wave.slice(np.array(ua))
        vph, y = wave.to_physical(sl, va, za)
        m, mr = float(sl["m"]), float(sl["mr"])
        if not reconstructed:  # canonical coordinates, no calibration at all
            return {
                "u": ua,
                "v": float(vph),
                "y": complex(y),
                "state": (ua, m, va, abs(za) ** 2),
            }
        eta = va - abs(za) ** 2 / (2 * m)
        zeta = 1 - va - abs(za) ** 2 / (2 * mr)
        al, cr = wave.chain(-0.5, ua, pairs=True)
        ar = wave.chain(ua, 0.5)[0]
        s = ua + 0.5
        p = al * wave.alpha0(s) * eta**2 / (self.gamma_ratio * wave.gamma0)
        f = ar * wave.alpha0(1 - s) * zeta**2 / (self.gamma_ratio * wave.gamma0)
        state = tuple(float(x) for x in wave.recon(p, f, wave.p0(s) * cr))
        return {"u": ua, "v": float(vph), "y": complex(y), "state": state}

    def evaluate(self, anchors, omegas, nv=24, nt=24, nang=8, reconstructed=True):
        """Return (S(omega)/mass for each omega, diagnostics).

        reconstructed=True is the quantity of the challenge: every state,
        the weight and the Fourier phase come from the observer-side
        inverse calibration of exact wave order probabilities.
        reconstructed=False is a diagnostic: canonical coordinates are
        used directly and interval volumes are replaced by the exact wave
        slack (pure first-order geometry, no calibration error).
        """
        wave = self.wave
        data = [self.anchor_state(a, reconstructed) for a in anchors]
        gram = np.diag([d["state"][3] for d in data])
        for i in range(3):
            for j in range(i + 1, 3):
                slack = float(
                    wave.margin(
                        data[i]["u"],
                        data[i]["v"],
                        data[i]["y"],
                        data[j]["u"],
                        data[j]["v"],
                        data[j]["y"],
                    )
                )
                if slack <= 0:
                    raise RuntimeError("anchor pair not timelike")
                if reconstructed:
                    ratio = wave.chain(data[i]["u"], data[j]["u"])[0]
                    b = (
                        ratio
                        * wave.alpha0(data[j]["u"] - data[i]["u"])
                        * slack**2
                        / (self.gamma_ratio * wave.gamma0)
                    )
                    gram[i, j] = gram[j, i] = float(
                        wave.gram(data[i]["state"], data[j]["state"], b)
                    )
                else:
                    si, sj = data[i]["state"], data[j]["state"]
                    gram[i, j] = gram[j, i] = 0.5 * (si[3] + sj[3]) - (
                        sj[1] - si[1]
                    ) * (sj[2] - si[2] - slack)
        ell = np.r_[-np.linalg.solve(gram[:2, :2], gram[:2, 2]), 1.0]
        vt, wv = composite(V_LO, V_HI, nv // 6, 6)
        tt, wt = composite(T_LO, T_HI, nt // 6, 6)
        ang = 2 * np.pi * np.arange(nang) / nang
        zeta = np.sqrt(tt)[None, :, None] * np.exp(1j * ang)[None, None, :]
        vgrid = vt[:, None, None]
        tgrid = tt[None, :, None]
        area = wv[:, None] * wt[None, :]
        totals = np.zeros(len(omegas), dtype=complex)
        mass = 0.0
        min_slack = np.inf
        for j, u in enumerate(self.u):
            sl = {k: v[j] for k, v in self.sl.items()}
            m, mr = sl["m"], sl["mr"]
            eta = vgrid - tgrid / (2 * m)
            zet = 1 - vgrid - tgrid / (2 * mr)
            if reconstructed:
                uh = self.u_hat[j]
                mh = float(wave.m(uh))
                th = mh / wave.M
                eh, zh = eta * self.eta_scale[j], zet * self.zeta_scale[j]
                vh = 1 - th + th * eh - (1 - th) * zh
                rh = 2 * mh * (wave.M - mh) / wave.M * (1 - eh - zh)
            else:  # canonical coordinates used directly (diagnostic only)
                uh, mh, vh, rh = u, m, vgrid + 0 * tgrid, tgrid + 0 * vgrid
            xu = (uh - 0.2) / 0.1
            wu_bump = math.exp(-1.0 / (1.0 - xu * xu)) if abs(xu) < 1 else 0.0
            wgt = wu_bump * bump((vh - 0.7) / 0.05) * bump((rh - 0.02125) / 0.01875)
            wgt = wgt[..., 0]
            if not np.any(wgt):
                continue
            y = sl["c0"] * zeta + sl["tb"] * np.conj(zeta)
            vph = eta + 0.5 * sqf(sl["aL"], sl["eL"], y)
            h = 0.0
            for i, d in enumerate(data):
                a, b = wave.inv_action(d["u"], u)
                slack = vph - d["v"] - 0.5 * sqf(a, b, y - d["y"])
                min_slack = min(min_slack, float(np.min(slack[wgt > 0])))
                if reconstructed:
                    key = (d["u"], j)
                    if key not in self._ax:
                        self._ax[key] = wave.chain(d["u"], u)[0]
                    ratio = self._ax[key]
                    bax = (
                        ratio
                        * wave.alpha0(u - d["u"])
                        * slack**2
                        / (self.gamma_ratio * wave.gamma0)
                    )
                    slack0 = np.sqrt(
                        bax * wave.gamma0 / wave.alpha0(uh - d["state"][0])
                    )
                else:
                    slack0 = slack
                g = 0.5 * (d["state"][3] + rh) - (mh - d["state"][1]) * (
                    vh - d["state"][2] - slack0
                )
                h = h + ell[i] * g
            hbar = np.mean(h, axis=-1)
            measure = self.wu[j] * wave.a2(u) * sl["detT"] * area * wgt
            mass += float(np.sum(measure))
            base = np.sum(measure * hbar)
            for n, om in enumerate(omegas):
                totals[n] += base * np.exp(-1j * om * uh)
        return totals / mass, {
            "ell": ell,
            "gram": gram,
            "min_slack": min_slack,
            "mass": mass,
        }


def d_anchor(wave, anchors, u):
    """D_A(u) of gram_fourier_noise.md Section 6, exact canonical vectors."""
    z = np.array([a[2] for a in anchors])
    g = (np.conj(z)[:, None] * z[None, :]).real
    ell = np.r_[-np.linalg.solve(g[:2, :2], g[:2, 2]), 1.0]
    ma = wave.m(np.array([a[0] for a in anchors]))
    gaps = wave.m(u)[:, None] - ma[None, :]
    return -0.5 * np.sum(ell[None, :] * np.conj(z)[None, :] ** 2 / gaps, axis=1)


def prediction(wave, anchors, offset=0.0):
    """Boxed limit (4) of the challenge, with a frequency offset omega - 1/lam."""
    u, w = composite(U_LO, U_HI, 200, 6)
    w = w * bump((u - 0.2) / 0.1)
    coef = np.sum(w * d_anchor(wave, anchors, u) * np.exp(-1j * offset * u))
    return -1j * wave.k * np.exp(1j * wave.phi) * coef / np.sum(w * wave.a2(u))


def anchors_from_gram(
    g11,
    g22,
    g12,
    g13,
    g23,
    du=(0.0, 0.0, 0.0),
    dv=(0.0, 0.0, 0.0),
    angle=0.0,
    reflect=False,
):
    """Anchor triple with prescribed Gram entries; returns (anchors, g33).

    Standard orientation a1 = (sqrt(g11), 0), a2 with positive second
    coordinate, then rotated by angle and optionally reflected.
    """
    a1 = complex(math.sqrt(g11), 0.0)
    a2x = g12 / a1.real
    a2 = complex(a2x, math.sqrt(g22 - a2x * a2x))
    a3x = g13 / a1.real
    a3 = complex(a3x, (g23 - a2.real * a3x) / a2.imag)
    z = np.array([a1, a2, a3]) * np.exp(1j * angle)
    if reflect:
        z = np.conj(z)
    base_u, base_v = (-0.35, -0.25, -0.15), (0.15, 0.25, 0.35)
    anchors = [(base_u[i] + du[i], base_v[i] + dv[i], z[i]) for i in range(3)]
    return anchors, abs(a3) ** 2
