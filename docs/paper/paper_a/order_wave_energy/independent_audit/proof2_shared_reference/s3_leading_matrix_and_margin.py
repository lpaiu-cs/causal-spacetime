# ruff: noqa: E501, UP031
"""Referee check of Section 3 of local_null_challenge.md (geometry of the witness pairs).

All matrices come from direct quadrature of the exact primitive (audit_geometry.Wave);
nothing is taken from the author's code.

A  Integration-by-parts formula for F_lambda (eq. (1) of null_coordinate_results.md).
B  Leading matrix (5): value, sign, factor 2 i sqrt(q) f(u)/(m (M-m)), O(lambda) remainder.
C  Uniform lower bound of the leading norm for tau in [pi/2, 3pi/2].
D  Slice map: both slacks preserved exactly; v_wave = v_tilde - z^T F z / 2 + O(lambda^4);
   y = z + O(lambda^2); Jacobian 1 + O(lambda^4).
E  Equal canonical z and equal v_tilde: exact causal margin = -lambda^2 z^T L z / 2 + O(lambda^3),
   transverse difference O(lambda^2), action O(lambda^3), fixed fraction of directions.
F  Mechanism: along a null generator (same v, same y) the margin is exactly 0 while
   v_tilde changes by lambda^2 z^T L z / 2.
G  Admissible neighbourhood: transverse ball radius c' lambda^(3/2), null-time window
   kappa lambda^2: every sampled pair is causally related with decreasing v_tilde.
H  Conditional mass of violating partners scales as lambda^6 (scaled-box Monte Carlo).

Deterministic (fixed seeds).  Prints PASS/FAIL lines; exits nonzero on any FAIL.
Run:  python s3_leading_matrix_and_margin.py      (system Python: numpy only)
"""

from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_geometry import A_TIP, B_TIP, Wave  # noqa: E402

FAILS = []
U_SET = (-0.3, 0.3)
TAU_SET = (1.7, 4.5)  # compact subinterval of (pi/2, 3pi/2)


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def dens(w, u):
    return 1.0 / math.cos(w.k * u) ** 2


def big_d(w, u):
    m = float(w.m(A_TIP, u))
    return m * (w.Mtot - m)


def leading(w, u, tau):
    """Right-hand side of (5) without the O(lambda) term (complex representative)."""
    return (
        2j
        * w.k
        * dens(w, u)
        / big_d(w, u)
        * cmath.exp(1j * (u / w.lam + w.phase))
        * (cmath.exp(1j * tau) - 1)
    )


def quad(c, z):
    """z^T S(c) z for the trace-free symmetric matrix with complex representative c."""
    return c.real * (z[0] ** 2 - z[1] ** 2) + 2 * c.imag * z[0] * z[1]


def part_ab():
    for q in (1.0, 0.25):
        errs_a, errs_b, wrong_sign, wrong_conj = {}, {}, {}, {}
        for r in (256, 1024, 4096):
            for phase in (0.37, 2.9):
                w = Wave(q, 1.0 / r, phase=phase)
                ea = eb = 0.0
                ws = wc = 1e9
                for u in np.linspace(U_SET[0], U_SET[1], 25):
                    m = float(w.m(A_TIP, u))
                    t = m / w.Mtot
                    pred = (
                        2j
                        * w.k
                        / big_d(w, u)
                        * (
                            dens(w, u) * cmath.exp(1j * u * r)
                            - (1 - t) * dens(w, A_TIP) * cmath.exp(1j * A_TIP * r)
                            - t * dens(w, B_TIP) * cmath.exp(1j * B_TIP * r)
                        )
                        * cmath.exp(1j * phase)
                    )
                    ea = max(ea, abs(w.f_matrix_complex(u) * r * r - pred))
                    for tau in (1.6, 2.4, math.pi, 4.0, 4.7):
                        exact = (
                            (w.f_matrix_complex(u + tau / r) - w.f_matrix_complex(u))
                            * r
                            * r
                        )
                        lead = leading(w, u, tau)
                        eb = max(eb, abs(exact - lead))
                        ws = min(ws, abs(exact + lead) / abs(lead))
                        wc = min(wc, abs(exact - lead / 1j) / abs(lead))
                errs_a[r] = max(errs_a.get(r, 0.0), ea)
                errs_b[r] = max(errs_b.get(r, 0.0), eb)
                wrong_sign[r] = min(wrong_sign.get(r, 1e9), ws)
                wrong_conj[r] = min(wrong_conj.get(r, 1e9), wc)
        # The O(lambda) coefficient oscillates with the endpoint phases, so the test is
        # boundedness of error/lambda (a slower rate would make it grow by 4x per rung).
        scaled_a = [errs_a[r] * r for r in (256, 1024, 4096)]
        check(
            f"A  q={q}: F_lambda/lambda^2 = integration-by-parts formula + O(lambda)",
            max(scaled_a) < 2.5 * min(scaled_a) and errs_a[4096] < 0.3 * errs_a[256],
            "sup error r=256,1024,4096: "
            + ", ".join(f"{errs_a[r]:.4f}" for r in (256, 1024, 4096))
            + "; error/lambda: "
            + ", ".join(f"{x:.2f}" for x in scaled_a),
        )
        scaled_b = [errs_b[r] * r for r in (256, 1024, 4096)]
        check(
            f"B1 q={q}: eq.(5) leading matrix, remainder O(lambda) (both phases, 5 tau, 25 u)",
            max(scaled_b) < 2.5 * min(scaled_b) and errs_b[4096] < 0.3 * errs_b[256],
            "sup error r=256,1024,4096: "
            + ", ".join(f"{errs_b[r]:.4f}" for r in (256, 1024, 4096))
            + "; error/lambda: "
            + ", ".join(f"{x:.0f}" for x in scaled_b),
        )
        check(
            f"B2 q={q}: sign and factor i of (5) (opposite sign is 2 off, a missing i is sqrt 2 off)",
            wrong_sign[4096] > 1.9 and wrong_conj[4096] > 1.3,
            f"min relative distance at r=4096: opposite sign {wrong_sign[4096]:.3f}, "
            f"real prefactor {wrong_conj[4096]:.3f}",
        )


def part_c():
    for q in (1.0, 0.25):
        w = Wave(q, 1.0 / 4096)
        ratio = 1e9
        for u in np.linspace(U_SET[0], U_SET[1], 61):
            floor = 2 * w.k * dens(w, u) * math.sqrt(2.0) / big_d(w, u)
            for tau in np.linspace(math.pi / 2, 3 * math.pi / 2, 41):
                exact = (
                    abs(w.f_matrix_complex(u + tau * w.lam) - w.f_matrix_complex(u))
                    / w.lam**2
                )
                ratio = min(ratio, exact / floor)
        check(
            f"C  q={q}: |F(u')-F(u)|/lambda^2 >= 2 sqrt(2q) f/(m(M-m)) (1-o(1)) on the u,tau sets",
            ratio > 0.97,
            f"min exact/floor at r=4096: {ratio:.4f}",
        )


def part_de():
    for q in (1.0, 0.25):
        rows = {}
        for r in (64, 256, 1024, 4096):
            w = Wave(q, 1.0 / r)
            lam = w.lam
            slack = exp4 = yz = jac = merr = act = dyn = gen = 0.0
            frac = 1.0
            for u in (-0.25, -0.1, 0.05, 0.2, 0.3):
                m = float(w.m(A_TIP, u))
                jac = max(jac, abs(np.linalg.det(w.slice_matrix(u)) - 1.0))
                f_u = w.f_matrix_complex(u)
                for tau in (1.7, math.pi, 4.5):
                    lead = leading(w, u, tau)
                    good = 0
                    for ang in np.linspace(0, math.pi, 36, endpoint=False):
                        z = 0.2 * np.array([math.cos(ang), math.sin(ang)])
                        first = w.from_canonical(u, 0.5, z)
                        second = w.from_canonical(u + tau * lam, 0.5, z)
                        eta, zeta = w.slacks(first)
                        slack = max(
                            slack,
                            abs(eta - (0.5 - 0.04 / (2 * m))),
                            abs(zeta - (0.5 - 0.04 / (2 * (w.Mtot - m)))),
                            abs(w.v_tilde(first) - 0.5),
                        )
                        exp4 = max(
                            exp4, abs(first[1] - (0.5 - 0.5 * quad(f_u, z))) / lam**4
                        )
                        yz = max(yz, np.linalg.norm(first[2:] - z) / lam**2)
                        margin = w.margin(first, second)
                        pred = -0.5 * quad(lead, z) * lam**2
                        merr = max(merr, abs(margin - pred) / lam**3)
                        dy = second[2:] - first[2:]
                        dyn = max(dyn, np.linalg.norm(dy) / lam**2)
                        act = max(
                            act,
                            0.5
                            * dy
                            @ np.linalg.solve(w.mt(first[0], second[0]), dy)
                            / lam**3,
                        )
                        good += margin >= 0.2 * abs(lead) * 0.04 * lam**2
                        # F: move along the null generator instead (same v, same y)
                        shifted = first.copy()
                        shifted[0] = u + tau * lam
                        gen = max(
                            gen,
                            abs(w.margin(first, shifted)),
                            abs(
                                (w.v_tilde(shifted) - 0.5)
                                - 0.5 * quad(lead, z) * lam**2
                            )
                            / lam**3,
                        )
                    frac = min(frac, good / 36)
            rows[r] = (slack, exp4, yz, jac, merr, act, dyn, frac, gen)
        check(
            f"D1 q={q}: slice map preserves both slacks and v_tilde exactly",
            max(rows[r][0] for r in rows) < 1e-13,
            f"max deviation {max(rows[r][0] for r in rows):.2e}",
        )
        check(
            f"D2 q={q}: v_wave = v_tilde - z^T F z/2 + O(lambda^4), y = z + O(lambda^2), |z|=0.2",
            max(rows[r][1] for r in rows) < 5 and max(rows[r][2] for r in rows) < 3,
            "remainder/lambda^4: "
            + ", ".join(f"{rows[r][1]:.2f}" for r in rows)
            + "; |y-z|/lambda^2: "
            + ", ".join(f"{rows[r][2]:.2f}" for r in rows),
        )
        scaled_j = [rows[r][3] * r**4 for r in rows]
        check(
            f"D3 q={q}: slice-map Jacobian = 1 + O(lambda^4) (bounded coefficient)",
            max(scaled_j) < 5 * min(scaled_j) and max(rows[r][3] for r in rows) < 1e-4,
            "|det-1|/lambda^4 at r=64,256,1024,4096: "
            + ", ".join(f"{x:.1f}" for x in scaled_j),
        )
        check(
            f"E1 q={q}: equal z, equal v_tilde: margin = -lambda^2 z^T L z/2 + O(lambda^3)",
            max(rows[r][4] for r in rows) < 15,
            "|margin-pred|/lambda^3 at r=64,256,1024,4096: "
            + ", ".join(f"{rows[r][4]:.2f}" for r in rows),
        )
        check(
            f"E2 q={q}: transverse difference O(lambda^2), action O(lambda^3)",
            max(rows[r][6] for r in rows) < 3 and max(rows[r][5] for r in rows) < 1,
            "|dy|/lambda^2: "
            + ", ".join(f"{rows[r][6]:.2f}" for r in rows)
            + "; action/lambda^3: "
            + ", ".join(f"{rows[r][5]:.3f}" for r in rows),
        )
        check(
            f"E3 q={q}: >= 1/3 of directions have margin >= 0.2 |L| |z|^2 lambda^2 (r >= 256)",
            min(rows[r][7] for r in (256, 1024, 4096)) >= 1 / 3,
            "fractions: " + ", ".join(f"{rows[r][7]:.3f}" for r in rows),
        )
        check(
            f"F  q={q}: null generator: margin 0, v_tilde shift = lambda^2 z^T L z/2 + O(lambda^3)",
            max(rows[r][8] for r in rows) < 15,
            "max(|margin|, |shift-pred|/lambda^3): "
            + ", ".join(f"{rows[r][8]:.2f}" for r in rows),
        )


def constants(q):
    """Explicit kappa, c' for the neighbourhood (no lambda dependence)."""
    w = Wave(q, 1.0 / 64)
    lmin = min(
        2 * w.k * dens(w, u) * 2 * math.sin(tau / 2) / big_d(w, u)
        for u in np.linspace(U_SET[0], U_SET[1], 121)
        for tau in np.linspace(TAU_SET[0], TAU_SET[1], 57)
    )
    z_lo = 0.15
    kappa = 0.25 * z_lo**2 * lmin / 3.0  # leading margin >= 3 kappa lambda^2
    c_prime = math.sqrt(TAU_SET[0] * kappa)  # action <= kappa lambda^2 / 2 (f >= 1)
    return lmin, kappa, c_prime


def sample_pair(w, rng, kappa, c_prime):
    lam = w.lam
    u = rng.uniform(*U_SET)
    tau = rng.uniform(*TAU_SET)
    lead = leading(w, u, tau)
    # directions with -z^T L z >= |L||z|^2/2  (one third of all directions)
    psi = 0.5 * (cmath.phase(lead) + math.pi) + rng.uniform(-math.pi / 6, math.pi / 6)
    rad = rng.uniform(0.15, 0.25)
    z1 = rad * np.array([math.cos(psi), math.sin(psi)])
    vt1 = rng.uniform(0.4, 0.6)
    rho, ang = (
        c_prime * lam**1.5 * math.sqrt(rng.uniform()),
        rng.uniform(0, 2 * math.pi),
    )
    z2 = z1 + rho * np.array([math.cos(ang), math.sin(ang)])
    vt2 = vt1 - kappa * lam * lam * rng.uniform(0.5, 1.0)
    first = w.from_canonical(u, vt1, z1)
    second = w.from_canonical(u + tau * lam, vt2, z2)
    return first, second, vt1 - vt2


def part_g():
    for q in (1.0, 0.25):
        lmin, kappa, c_prime = constants(q)
        rng = np.random.default_rng(2026)
        stats = {}
        for r in (64, 256, 1024, 4096, 16384):
            w = Wave(q, 1.0 / r)
            margins, interior = [], 1.0
            for _ in range(1200):
                first, second, _ = sample_pair(w, rng, kappa, c_prime)
                margins.append(w.margin(first, second) * r * r / kappa)
                interior = min(interior, *w.slacks(first), *w.slacks(second))
            margins = np.array(margins)
            stats[r] = (float(np.mean(margins > 0)), float(margins.min()), interior)
        print(
            f"INFO q={q}: L_min={lmin:.3f}, kappa={kappa:.5f}, c'={c_prime:.4f}; "
            "fraction of sampled pairs with positive exact margin: "
            + ", ".join(f"r={r}: {stats[r][0]:.3f}" for r in stats)
        )
        check(
            f"G1 q={q}: every pair of the neighbourhood is causally related for r >= 1024",
            all(stats[r][0] == 1.0 for r in (1024, 4096, 16384)),
            "min margin/(kappa lambda^2): "
            + ", ".join(f"r={r}: {stats[r][1]:.3f}" for r in stats),
        )
        check(
            f"G2 q={q}: sampled margins respect the derived floor 1.5 kappa lambda^2 (r >= 1024)",
            all(stats[r][1] > 1.5 for r in (1024, 4096, 16384)),
            f"r=1024: {stats[1024][1]:.3f}, r=4096: {stats[4096][1]:.3f}, "
            f"r=16384: {stats[16384][1]:.3f}; "
            f"min slack of sampled events {min(stats[r][2] for r in stats):.3f}",
        )


def part_h():
    """Scaled-box estimate of the conditional mass of violating partners of a fixed X1."""
    q = 1.0
    _, kappa, _ = constants(q)
    box_tau, box_w, box_d = 4 * math.pi, 6.0 * kappa, 0.6
    volume = box_tau * box_w * (2 * box_d) ** 2
    out = {}
    for r in (256, 1024, 4096):
        w = Wave(q, 1.0 / r)
        lam = w.lam
        rng = np.random.default_rng(31)
        u1, vt1 = 0.1, 0.5
        lead = leading(w, u1, math.pi)
        psi = 0.5 * (cmath.phase(lead) + math.pi)
        z1 = 0.2 * np.array([math.cos(psi), math.sin(psi)])
        first = w.from_canonical(u1, vt1, z1)
        hits, n = 0, 6000
        for _ in range(n):
            tau = rng.uniform(0, box_tau)
            gap = rng.uniform(0, box_w)
            dz = rng.uniform(-box_d, box_d, 2)
            second = w.from_canonical(
                u1 + tau * lam, vt1 - gap * lam**2, z1 + dz * lam**1.5
            )
            hits += (w.margin(first, second) > 0) and gap >= 0.5 * kappa
        out[r] = hits / n * volume
    vals = [out[r] for r in (256, 1024, 4096)]
    check(
        "H  conditional mass of violating partners / lambda^6 is stable in lambda (q=1)",
        min(vals) > 0 and max(vals) / min(vals) < 1.5,
        "scaled mass (canonical coordinates, unit Jacobian) at r=256,1024,4096: "
        + ", ".join(f"{v:.5f}" for v in vals),
    )


if __name__ == "__main__":
    part_ab()
    part_c()
    part_de()
    part_g()
    part_h()
    print(f"SUMMARY s3_leading_matrix_and_margin: {len(FAILS)} FAIL")
    sys.exit(1 if FAILS else 0)
