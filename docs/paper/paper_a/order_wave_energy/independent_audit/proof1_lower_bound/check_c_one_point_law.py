"""Referee check (c): one-point laws differ by O(lambda^2) in total variation.

The wave diamond D_lam and the background diamond D_0 share tips (u=-1/2, v=0,
y=0) and (u=1/2, v=h, y=0) in the y chart; both one-point laws have density
proportional to a(u)^2 (the constant c cancels), so
    TV(mu_lam, mu_0) <= max(A(D_lam \\ D_0), A(D_0 \\ D_lam)) / V_0,
with A the a^2-weighted volume.  Referee's explicit bound (referee report):
    A(D_lam symdiff D_0) <= B(lam) := pi h^2 (1+t)^2/(1-t) int a^2 (e_- + e_+) du
                        <= 2 pi h^2 (1+t)^2/(1-t) C_E lam^2,   t = tanh(2 beta),
with e_-(u) = ||E(-1/2,u)||, e_+(u) = ||E(u,1/2)||.

For each (u, y) the v-sections of both diamonds are intervals, so the section
of the symmetric difference is computed exactly; only (u, y) is Monte Carlo
(common random numbers for both diamonds, fixed seed).

Tests
  C0  MC value of V_0 against the closed form pi h^2 C2/K.
  C1  A(symdiff) <= B(lam) and B(lam)/lam^2 <= explicit constant.
  C2  exact-section TV estimate <= B(lam)/V_0, and TV/lam^2 stays bounded as
      lam decreases (no growth by more than 25 percent between rungs).

Run: python check_c_one_point_law.py
"""

from __future__ import annotations

import numpy as np
from audit_common import Family, Reporter

rep = Reporter("check (c): one-point law, domain-strip O(lambda^2) bound")

NS = 1500000
rng = np.random.default_rng(31415926)

for q in (0.25, 1.0, 4.0):
    prev = None
    for lam in (1 / 8, 1 / 16, 1 / 32, 1 / 64, 1 / 128):
        fam = Family(q, lam, phi=float(rng.uniform(0, 2 * np.pi)))
        h = fam.h
        t = fam.tanh2b
        R = np.sqrt(0.5 * h * fam.Mtot * (1 + t))
        box = 1.0 * (2 * R) ** 2
        u = rng.uniform(-0.5, 0.5, NS)
        y1 = rng.uniform(-R, R, NS)
        y2 = rng.uniform(-R, R, NS)
        r2 = y1 * y1 + y2 * y2
        a2 = fam.a(u) ** 2
        m_lo = fam.m(u)
        m_hi = fam.Mtot - m_lo
        P = fam.prim(u)
        Pend = fam.prim(np.array([0.5]))[0]
        G_lo = np.exp(1j * fam.phi) * P  # segment (-1/2, u)
        G_hi = np.exp(1j * fam.phi) * (Pend - P)  # segment (u, 1/2)
        L_lo_w = 0.5 * fam.wave_form(m_lo, G_lo, y1, y2)
        L_hi_w = 0.5 * fam.wave_form(m_hi, G_hi, y1, y2)
        L_lo_0 = 0.5 * r2 / m_lo
        L_hi_0 = 0.5 * r2 / m_hi
        # v-sections [L_lo, h - L_hi]
        lo_w, hi_w = L_lo_w, h - L_hi_w
        lo_0, hi_0 = L_lo_0, h - L_hi_0
        len_w = np.maximum(hi_w - lo_w, 0.0)
        len_0 = np.maximum(hi_0 - lo_0, 0.0)
        inter = np.maximum(np.minimum(hi_w, hi_0) - np.maximum(lo_w, lo_0), 0.0)
        inter = np.where((len_w > 0) & (len_0 > 0), inter, 0.0)
        A_int = box * np.mean(a2 * inter)
        A_w0 = box * np.mean(a2 * (len_w - inter))  # in wave diamond only
        A_0w = box * np.mean(a2 * (len_0 - inter))  # in background diamond only
        V0_mc = A_int + A_0w
        Vw_mc = A_int + A_w0
        V0 = fam.reference_volume()
        tv = 0.5 * (
            A_int * abs(V0_mc - Vw_mc) / (V0_mc * Vw_mc) + A_w0 / Vw_mc + A_0w / V0_mc
        )
        # analytic bound B(lam) by deterministic quadrature in u
        ug = -0.5 + (np.arange(20000) + 0.5) / 20000
        Pg = fam.prim(ug)
        e_lo = t * np.abs(Pg)
        e_hi = t * np.abs(Pend - Pg)
        B = (
            np.pi
            * h
            * h
            * (1 + t) ** 2
            / (1 - t)
            * np.mean(fam.a(ug) ** 2 * (e_lo + e_hi))
        )
        Bconst = 2 * np.pi * h * h * (1 + t) ** 2 / (1 - t) * fam.C_E
        tag = f"q={q:.2f} lam=1/{round(1 / lam)}"
        rep.check(
            "C0 V_0 (MC) vs pi h^2 C2/K " + tag,
            abs(V0_mc / V0 - 1) < 5e-3,
            f"MC {V0_mc:.5f} exact {V0:.5f}",
        )
        sym = A_w0 + A_0w
        rep.check(
            "C1 A(symdiff) <= B(lam) <= const*lam^2 " + tag,
            sym <= B and B / lam**2 <= Bconst,
            f"A(symdiff)/lam^2={sym / lam**2:.4f} B/lam^2={B / lam**2:.4f} "
            f"const={Bconst:.2f}",
        )
        grow_ok = True if prev is None else (tv / lam**2 <= 1.25 * prev)
        rep.check(
            "C2 TV(mu_lam,mu_0) <= B/V_0 and TV/lam^2 not growing " + tag,
            tv <= B / V0 and grow_ok,
            f"TV={tv:.3e} TV/lam^2={tv / lam**2:.4f} bound/lam^2={B / V0 / lam**2:.3f}",
        )
        prev = tv / lam**2

rep.finish()
