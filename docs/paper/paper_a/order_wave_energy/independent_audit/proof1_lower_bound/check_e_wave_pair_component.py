"""Referee check (e), part 3: Section 4 on a real two-event component.

End-to-end test of "sandwich (4) + smoothing lemma => |sqrt p_lam - sqrt p_0|
<= C delta/h" in the EXACT wave family, for one uncertain pair of interior
cubes, conditional on the spatial coordinates (as in the proof), including the
case where the reference probability is exactly zero.

Setup: cube side hc in the (T,Z,y) chart; the two events have fixed spatial
coordinates z_1, z_2 (distance d) and T coordinates with conditional density
proportional to a^4 on intervals of length hc.  With s = (max of T_2 - T_1) - d,
the background relation probability is p_0 = P(T_2 - T_1 >= d): zero for s<=0,
about s^2/(2 hc^2) for small s>0.

For every sampled pair the script computes the exact wave relation.  Let
eps = 2 sqrt(2) max e(u_1,u_2) over the sampled pairs (the pointwise form of
the sandwich proved in the referee report; eps <= c0 delta).  Tests:
  G1  sample-wise sandwich: {T_2-T_1 >= d+eps} subset {wave relation} subset
      {T_2-T_1 >= d-eps}  (zero violations, hence the same for probabilities).
  G2  |sqrt p_lam - sqrt p_0| <= sqrt(L/2) eps, L = (2 rbar + rbar')/hc^2 with
      rbar, rbar' from the actual conditional densities (referee's constant).
  G3  zero-reference-probability rows: p_0 = 0 exactly while p_lam > 0 occurs
      for at least one tested geometry, and still sqrt p_lam <= sqrt(L/2) eps.
  G4  the naive (TV-type) change |p_lam - p_0| is of order eps/hc while the
      squared-root change is of the same order: the Hellinger cost of the
      component is O((delta/hc)^2), as (7) needs.  Reported, with the
      requirement max |sqrt p_lam - sqrt p_0| * hc/delta <= sqrt(L hc^2 / 2) c0.

Monte Carlo with common random numbers and a fixed seed (the inequalities G1
are exact sample-wise; G2/G3 use weighted frequencies and a stated 5-sigma
binomial allowance).  Run: python check_e_wave_pair_component.py
"""

from __future__ import annotations

import numpy as np
from audit_common import SQRT2, Family, Reporter

rep = Reporter("check (e3): two-event component in the exact wave family")
rng = np.random.default_rng(161803)

NS = 500000
any_zero_ref_positive = False
worst_ratio = 0.0

for q, lam, hc in ((1.0, 1 / 64, 0.006), (1.0, 1 / 128, 0.01), (4.0, 1 / 128, 0.04)):
    fam = Family(q, lam, phi=0.9)
    c0d = fam.c0 * fam.delta
    assert hc >= c0d, "cell must satisfy hc >= c0 delta"
    # centre of the background diamond in the Minkowski chart
    Tm = (fam.Mtot + fam.h) / (2 * SQRT2)
    Zm = (fam.Mtot - fam.h) / (2 * SQRT2)
    for label, dz in (
        ("dy along e1", np.array([0.05, 0.22, 0.0])),
        ("dy along e2", np.array([0.05, 0.0, 0.22])),
        ("dy diagonal", np.array([-0.08, 0.15, 0.15])),
    ):
        d = float(np.linalg.norm(dz))
        z1 = np.array([Zm, 0.0, 0.0]) - 0.5 * dz
        z2 = z1 + dz
        for s_rel in (-0.5, -0.01, 0.0, 0.05, 0.2, 1.0):
            # s = s_rel * hc: maximal background margin available in the two cubes
            s = s_rel * hc
            dTc = d - hc + s
            Tc1 = Tm - 0.5 * dTc
            Tc2 = Tc1 + dTc
            T1 = Tc1 + hc * (rng.uniform(0, 1, NS) - 0.5)
            T2 = Tc2 + hc * (rng.uniform(0, 1, NS) - 0.5)
            m1, v1 = (T1 + z1[0]) / SQRT2, (T1 - z1[0]) / SQRT2
            m2, v2 = (T2 + z2[0]) / SQRT2, (T2 - z2[0]) / SQRT2
            u1, u2 = fam.u_of_m(m1), fam.u_of_m(m2)
            assert (
                np.all(u2 > u1)
                and np.all(np.abs(u1) < 0.5)
                and np.all(np.abs(u2) < 0.5)
            )
            w = fam.a(u1) ** 4 * fam.a(u2) ** 4
            b = m2 - m1
            G = fam.pair_G(fam.prim(u1), fam.prim(u2))
            wave = 2.0 * (v2 - v1) >= fam.wave_form(b, G, dz[1], dz[2])
            e = fam.tanh2b * np.abs(G)
            eps = 2.0 * SQRT2 * float(e.max())
            g = (T2 - T1) - d
            back = g >= 0.0
            tight = g >= eps
            relax = g >= -eps
            viol = int(np.sum(tight & ~wave) + np.sum(wave & ~relax))
            W = w.sum()
            p0 = float((w * back).sum() / W)
            pl = float((w * wave).sum() / W)
            # referee's Lipschitz constant from the actual conditional densities
            a4_1 = fam.a(u1) ** 4
            a4_2 = fam.a(u2) ** 4
            sup1, sup2 = a4_1.max() / a4_1.mean(), a4_2.max() / a4_2.mean()
            # |d/dT a^4| <= 4 sqrt(q) tan(sqrt(q)/2) a^4 du/dT,
            # with du/dT = a^2/sqrt2 <= 1/sqrt2
            dlog = 4.0 * fam.rq * np.tan(0.5 * fam.rq) / SQRT2 * hc
            rbar = sup1 * sup2
            Lc = (2.0 * rbar + dlog * rbar) / hc**2
            bound = np.sqrt(Lc / 2.0) * eps
            diff = abs(np.sqrt(pl) - np.sqrt(p0))
            # 5-sigma binomial allowance on the sqrt scale (delta method), tiny
            allow = (
                5.0
                * 0.5
                * np.sqrt(max(pl, p0, 1.0 / NS) / NS)
                / max(np.sqrt(max(pl, p0)), 1e-9)
                if max(pl, p0) > 0
                else 0.0
            )
            ok = viol == 0 and eps <= c0d and diff <= bound + allow
            if s <= 0.0:
                ok = ok and p0 == 0.0
                if pl > 0.0:
                    any_zero_ref_positive = True
            worst_ratio = max(
                worst_ratio,
                diff * hc / fam.delta / (np.sqrt(Lc * hc * hc / 2.0) * fam.c0),
            )
            rep.check(
                f"G1/G2 q={q:.0f} lam=1/{round(1 / lam)} hc={hc:.3f} {label} "
                f"s/hc={s_rel:+.2f}",
                ok,
                f"p0={p0:.3e} p_lam={pl:.3e} |dsqrt|={diff:.3e} bound={bound:.3e} "
                f"eps/delta={eps / fam.delta:.2f} (c0={fam.c0:.1f}) "
                f"sandwich violations={viol}",
            )

rep.check(
    "G3 a zero-reference-probability relation becomes possible in the wave and stays "
    "under the sqrt bound",
    any_zero_ref_positive,
)
rep.check(
    "G4 max |sqrt p_lam - sqrt p_0| hc/delta relative to sqrt(L hc^2/2) c0",
    worst_ratio <= 1.0,
    f"max ratio {worst_ratio:.4f}",
)
rep.finish()
