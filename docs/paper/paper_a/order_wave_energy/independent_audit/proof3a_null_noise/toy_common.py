"""Shared machinery for the 1+1 analogue used by the proof3a referee scripts.

Model (referee-designed, NOT the causal-set instrument of the paper):

* Events are iid uniform on the unit square in null coordinates (u, v); this
  is normalised volume in a 1+1 Minkowski causal diamond.  r < x means
  u_r < u_x and v_r < v_x.  Lines of fixed v are null generators.
* Raw order moments of a probe X=(u,v), all estimated from n SHARED references:
      p = Pr(R<X) = u v
      f = Pr(X<R) = (1-u)(1-v)
      a = Pr(R<X, R'<X, R~R') = p^2/2      (order-two U-statistic, all pairs)
      b = Pr(A<R<X) = (u-uA)(v-vA)         (anchor interval)
  and of one fixed anchor A: pA, fA.
* Smooth nonlinear calibration (valid where u<v, which holds on the window):
      s = 1+p-f, d = sqrt(s^2-4p), c = 2a/p^2 - 1   (c = 0 in the population)
      u_hat = (s-d)/2 + BETA*c,   v_hat = (s+d)/2
  and the same square-root map for the anchor.
* Population-zero functional, genuinely nonlinear in the raw moments:
      H = (u_hat-uA_hat)(v_hat-vA_hat) - b + KAPPA*c
* Score density  F_omega = w(u_hat, v_hat) H exp(-i omega u_hat), with w a
  smooth bump window strictly inside the diamond and to the future of A.

Probes are uniform on the block [U_LO,U_HI]x[V_LO,V_HI] (density 1 there).
A 1+1 diamond has a u<->v reflection symmetry that order data cannot break,
so the side u<v is supplied by restricting probes to this block; the 3+1
instrument of the paper needs no such input.  This does not touch the noise
mechanism under test.

The probe integral is evaluated on a tensor midpoint grid of probe positions.
The counts at every grid node are EXACT order counts; the only approximation
is the quadrature, which is validated against exact cell integration.
"""

from __future__ import annotations

import numpy as np

UA, VA = 0.03, 0.30  # fixed anchor
UC, VC, HW = 0.25, 0.75, 0.15  # window centre and half-width
BETA, KAPPA = 0.15, 0.3
U_LO, U_HI = 0.04, 0.46  # probe block (contains the window with margin)
V_LO, V_HI = 0.54, 0.96
TINY = 1e-300
PMIN = 1e-12  # guards 0/0 where no reference precedes the probe (weight is 0 there)


def bump(x):
    out = np.zeros_like(x, dtype=float)
    inside = np.abs(x) < 1
    out[inside] = np.exp(-1.0 / (1.0 - x[inside] ** 2))
    return out


def dominance_counts(u, v):
    """D[k] = #{j : u_j < u_k and v_j < v_k}, exact, O(n log^2 n), vectorised."""
    n = u.size
    order = np.argsort(u, kind="stable")
    vr = np.empty(n, dtype=np.int64)
    vr[np.argsort(v[order], kind="stable")] = np.arange(n, dtype=np.int64)
    idx = np.arange(n, dtype=np.int64)
    cnt = np.zeros(n, dtype=np.int64)
    size = 1
    while size < n:
        blk = idx // (2 * size)
        left = ((idx // size) & 1) == 0
        o = np.argsort(blk * n + vr, kind="stable")
        lo = left[o]
        c = np.cumsum(lo)
        b = blk[o]
        first = np.searchsorted(b, b, side="left")
        within = c - (c[first] - lo[first])
        sel = ~lo
        cnt[o[sel]] += within[sel]
        size *= 2
    out = np.empty(n, dtype=np.int64)
    out[order] = cnt
    return out


def calibrate(p, f, a, b, pa, fa):
    """Return (amplitude w*H, u_hat) from raw moments (arrays) and anchor moments."""
    s = 1.0 + p - f
    d = np.sqrt(np.maximum(s * s - 4.0 * p, TINY))
    c = 2.0 * a / np.maximum(p, PMIN) ** 2 - 1.0
    uh = 0.5 * (s - d) + BETA * c
    vh = 0.5 * (s + d)
    sa = 1.0 + pa - fa
    da = np.sqrt(max(sa * sa - 4.0 * pa, TINY))
    ua, va = 0.5 * (sa - da), 0.5 * (sa + da)
    h = (uh - ua) * (vh - va) - b + KAPPA * c
    w = bump((uh - UC) / HW) * bump((vh - VC) / HW)
    return w * h, uh


def grid(nu, nv):
    hu, hv = (U_HI - U_LO) / nu, (V_HI - V_LO) / nv
    ug = U_LO + hu * (np.arange(nu) + 0.5)
    vg = V_LO + hv * (np.arange(nv) + 0.5)
    return ug, vg, hu, hv


def population_moments(ug, vg):
    uu, vv = np.meshgrid(ug, vg, indexing="ij")
    p = uu * vv
    f = (1 - uu) * (1 - vv)
    a = 0.5 * p * p
    b = np.clip(uu - UA, 0, None) * np.clip(vv - VA, 0, None)
    return p, f, a, b, UA * VA, (1 - UA) * (1 - VA)


def empirical_moments(u, v, ug, vg):
    """Exact empirical raw moments at every grid node from n shared references."""
    n, nu, nv = u.size, ug.size, vg.size
    dom = dominance_counts(u, v)
    iu = np.searchsorted(ug, u, side="right")  # node a has ug[a] > u_r iff a >= iu
    iv = np.searchsorted(vg, v, side="right")
    ok = (iu < nu) & (iv < nv)
    flat = iu[ok] * nv + iv[ok]
    size = nu * nv

    def past(weights):
        h2 = np.bincount(flat, weights=weights, minlength=size).reshape(nu, nv)
        return np.cumsum(np.cumsum(h2, axis=0), axis=1)

    ones = np.ones(int(ok.sum()))
    pc = past(ones)
    bc = past(((u > UA) & (v > VA))[ok].astype(float))
    qc = past(dom[ok].astype(float))
    ju = np.searchsorted(ug, u, side="left") - 1  # node a has ug[a] < u_r iff a <= ju
    jv = np.searchsorted(vg, v, side="left") - 1
    okf = (ju >= 0) & (jv >= 0)
    h2 = np.bincount(ju[okf] * nv + jv[okf], minlength=size).reshape(nu, nv)
    fc = np.cumsum(np.cumsum(h2[::-1, ::-1], axis=0), axis=1)[::-1, ::-1]
    pa = np.count_nonzero((u < UA) & (v < VA)) / n
    fa = np.count_nonzero((u > UA) & (v > VA)) / n
    return pc / n, fc / n, qc / (0.5 * n * (n - 1)), bc / n, pa, fa


def fourier_sums(amp, uh, omegas, cell):
    """cell * sum_nodes amp * exp(-i omega u_hat) for every omega."""
    mask = amp != 0.0
    if not np.any(mask):
        return np.zeros(len(omegas), dtype=complex)
    vals, phase = amp[mask], uh[mask]
    return cell * (np.exp(-1j * np.outer(omegas, phase)) @ vals)


def scaled_scores(rng, n, ug, vg, hu, hv, omegas, ts, pop=None):
    """G[t, omega] = probe integral of F(mu_0 + t*(mu_hat - mu_0)) for one replicate.

    t=1 is the actual statistic.  Other t values rescale the SAME empirical
    error and are used only to separate Taylor orders.
    Also returns the probe second moment  int |w H|^2 dX  at t=1 (or nan).
    """
    if pop is None:
        pop = population_moments(ug, vg)
    u, v = rng.random(n), rng.random(n)
    emp = empirical_moments(u, v, ug, vg)
    out = np.empty((len(ts), len(omegas)), dtype=complex)
    second = np.nan
    for k, t in enumerate(ts):
        mixed = [m0 + t * (m1 - m0) for m0, m1 in zip(pop, emp, strict=True)]
        amp, uh = calibrate(*mixed)
        out[k] = fourier_sums(amp, uh, omegas, hu * hv)
        if t == 1:
            second = hu * hv * float(np.sum(amp * amp))
    return out, second


def exact_cell_scores(u, v, omegas):
    """Exact probe integral over the probe block: the integrand is constant on
    the cells cut out by the reference coordinates, so no quadrature is used.

    Returns (scores per omega, int |w H|^2 dX).  Memory is O(n^2) int32.
    """
    n = u.size
    order = np.argsort(u)
    us, vsort = u[order], np.sort(v)
    vrank = np.empty(n, dtype=np.int64)
    vrank[np.argsort(v[order])] = np.arange(1, n + 1)
    rows = np.arange(1, n + 1)
    e = np.zeros((n + 1, n + 1), dtype=np.int32)
    e[rows, vrank] = 1
    pc = np.cumsum(np.cumsum(e, axis=0, dtype=np.int32), axis=1, dtype=np.int32)
    dom = pc[rows - 1, vrank - 1]  # references strictly preceding each reference
    e[rows, vrank] = dom
    qc = np.cumsum(np.cumsum(e, axis=0, dtype=np.int64), axis=1, dtype=np.int64)
    del e
    # Probes are restricted to the block: the 1+1 toy has a u<->v reflection
    # symmetry, so the half u<v must be supplied (see module docstring).
    du = np.diff(np.clip(np.r_[0.0, us, 1.0], U_LO, U_HI))
    dv = np.diff(np.clip(np.r_[0.0, vsort, 1.0], V_LO, V_HI))
    ri, ci = np.flatnonzero(du > 0), np.flatnonzero(dv > 0)
    i, j = ri[:, None], ci[None, :]
    ia = int(np.count_nonzero(u < UA))
    ja = int(np.count_nonzero(v < VA))
    pblk = pc[np.ix_(ri, ci)].astype(float)
    qblk = qc[np.ix_(ri, ci)].astype(float)
    fblk = n - i - j + pblk
    bblk = pblk - pc[ia, ci][None, :] - pc[ri, ja][:, None] + pc[ia, ja]
    bblk = np.where((i >= ia) & (j >= ja), bblk, 0.0)
    pa = pc[ia, ja] / n
    fa = (n - ia - ja + pc[ia, ja]) / n
    amp, uh = calibrate(
        pblk / n, fblk / n, qblk / (0.5 * n * (n - 1)), bblk / n, pa, fa
    )
    area = du[ri][:, None] * dv[ci][None, :]
    mask = amp != 0.0
    vals, phase = (amp * area)[mask], uh[mask]
    scores = np.array([np.sum(vals * np.exp(-1j * w * phase)) for w in omegas])
    return scores, float(np.sum(amp * amp * area))


def first_order_variance(omega, nu=4200, nv=160):
    """Deterministic value of sigma_1(omega)^2 = n * E|first Taylor term|^2.

    The first-order term is (1/n) sum_j psi(R_j) with the influence function
        psi(r) = int_block w(X) DH(mu_0(X))[kernel(r; X)] exp(-i omega u) dX
    (centred).  Its variance over one uniform reference r is computed by
    quadrature: cumulative sums give psi on a tensor grid of r, with the
    probe cells' edges as r-nodes so that the jump of each indicator kernel
    sits exactly on a quadrature boundary.  No random numbers are used.
    """
    from numpy.polynomial.legendre import leggauss

    ug, vg, hu, hv = grid(nu, nv)
    u, v = np.meshgrid(ug, vg, indexing="ij")
    s, d, p = u + v, v - u, u * v
    du_dp = 0.5 * (1 - (s - 2) / d) - 2 * BETA / p
    du_df = 0.5 * (-1 + s / d)
    dv_dp = 0.5 * (1 + (s - 2) / d)
    dv_df = 0.5 * (-1 - s / d)
    hp = (v - VA) * du_dp + (u - UA) * dv_dp - 2 * KAPPA / p
    hf = (v - VA) * du_df + (u - UA) * dv_df
    ha = (v - VA) * 2 * BETA / p**2 + 2 * KAPPA / p**2
    sa, da = UA + VA, VA - UA
    hpa = -(v - VA) * 0.5 * (1 - (sa - 2) / da) - (u - UA) * 0.5 * (1 + (sa - 2) / da)
    hfa = -(v - VA) * 0.5 * (-1 + sa / da) - (u - UA) * 0.5 * (-1 - sa / da)
    w = bump((u - UC) / HW) * bump((v - VC) / HW)
    base = hu * hv * w * np.exp(-1j * omega * u)

    def upper(g):  # sum over probe cells with index >= edge index
        out = np.zeros((nu + 1, nv + 1), dtype=complex)
        out[:nu, :nv] = np.cumsum(np.cumsum(g[::-1, ::-1], axis=0), axis=1)[::-1, ::-1]
        return out

    def lower(g):  # sum over probe cells with index < edge index
        out = np.zeros((nu + 1, nv + 1), dtype=complex)
        out[1:, 1:] = np.cumsum(np.cumsum(g, axis=0), axis=1)
        return out

    sp_, sb = upper(base * hp), upper(-base)  # dH/db = -1
    lf = lower(base * hf)
    ga = base * ha
    a0, au, av, auv = upper(ga), upper(ga * u), upper(ga * v), upper(ga * u * v)
    gpa, gfa = np.sum(base * hpa), np.sum(base * hfa)
    x4, w4 = leggauss(4)

    def nodes(lo, hi, anchor, inner):
        xs, ws, idx = [], [], []
        for a_, b_ in ((0.0, anchor), (anchor, lo)):
            xs.append(0.5 * (a_ + b_) + 0.5 * (b_ - a_) * x4)
            ws.append(0.5 * (b_ - a_) * w4)
            idx.append(np.zeros(4, dtype=int))
        wt = np.full(inner + 1, (hi - lo) / inner)
        wt[0] *= 0.5
        wt[-1] *= 0.5
        xs.append(lo + (hi - lo) * np.arange(inner + 1) / inner)
        ws.append(wt)
        idx.append(np.arange(inner + 1))
        xs.append(0.5 * (hi + 1) + 0.5 * (1 - hi) * x4)
        ws.append(0.5 * (1 - hi) * w4)
        idx.append(np.full(4, inner))
        return np.concatenate(xs), np.concatenate(ws), np.concatenate(idx)

    ur, wu, iu = nodes(U_LO, U_HI, UA, nu)
    vr, wv, iv = nodes(V_LO, V_HI, VA, nv)
    ii, jj = np.ix_(iu, iv)
    urr, vrr = ur[:, None], vr[None, :]
    after = (urr > UA) & (vrr > VA)
    before = (urr < UA) & (vrr < VA)
    # ancestor-pair h_1 kernel: 2*1[r<X]*(u_r v_r + (u-u_r)(v-v_r))
    pair = 2 * (
        2 * urr * vrr * a0[ii, jj] - urr * av[ii, jj] - vrr * au[ii, jj] + auv[ii, jj]
    )
    psi = (
        sp_[ii, jj]
        + lf[ii, jj]
        + after * sb[ii, jj]
        + pair
        + before * gpa
        + after * gfa
    )
    wr = wu[:, None] * wv[None, :]
    mean = np.sum(wr * psi)
    return float(np.sum(wr * np.abs(psi) ** 2) - abs(mean) ** 2)
