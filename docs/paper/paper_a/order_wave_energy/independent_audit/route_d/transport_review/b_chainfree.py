# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""Chain-free checks of the transport draft in 3+1 (no chain is run).

Sampling, event sets, bins, linear projection of the 7-parameter group fit,
induced correlation of the readout under iid/heteroscedastic noise, and the
C2 Brownian field (expected values by projection).
"""
import math, sys, pathlib, time
import numpy as np

AUD = pathlib.Path(r"E:/lab/causal-spacetime/.claude/worktrees/upbeat-newton-3e92fe/docs/paper/paper_a/order_wave_energy/independent_audit")
sys.path.insert(0, str(AUD / "route_d"))
sys.path.insert(0, str(AUD / "route_a"))
import flip_probe as fp  # noqa
import drift_probe as dp  # noqa

EDGES = (0.0, 0.25, 0.5, 0.75, np.inf)
O = (dp.TIP_P + dp.TIP_Q) / 2


def scale_of(x):
    """Smallest sub-diamond scale s (about the midpoint) containing x."""
    half = (dp.TIP_Q - dp.TIP_P) / 2
    lo_s, hi_s = np.zeros(len(x)), np.ones(len(x)) * 1.0001
    for _ in range(40):
        mid = (lo_s + hi_s) / 2
        low = O - mid[:, None] * half
        high = O + mid[:, None] * half
        d1, d2 = x - low, high - x
        ok = (d1[:, 0] > np.linalg.norm(d1[:, 1:], axis=1)) & (d2[:, 0] > np.linalg.norm(d2[:, 1:], axis=1))
        hi_s = np.where(ok, mid, hi_s)
        lo_s = np.where(ok, lo_s, mid)
    return hi_s


def tangent(x):
    h = 1e-6
    cols = []
    for i in range(7):
        p = np.zeros(7); p[i] = h
        m = np.zeros(7); m[i] = -h
        cols.append(((dp.group_map(p, x) - dp.group_map(m, x)) / (2 * h)).ravel())
    return np.array(cols).T  # (4n, 7)


def proj_complement(x):
    J = tangent(x)
    Q, _ = np.linalg.qr(J)
    return np.eye(J.shape[0]) - Q @ Q.T


def bin_stats(cov, x, wave=False):
    """Expected S, I per bin from the 4n x 4n residual covariance."""
    n = len(x)
    C = cov.reshape(n, 4, n, 4)
    if wave:
        w = np.array([1, 1, 0, 0]) / math.sqrt(2)
        Cp = np.einsum("a,iajb,b->ij", w, C, w)
    else:
        Cp = np.einsum("iaja->ij", C)
    var = np.diag(Cp)
    S = var[:, None] + var[None, :] - 2 * Cp
    I = var[:, None] + var[None, :]
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    up = np.triu(np.ones((n, n), bool), 1)
    out = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        pick = up & (d >= lo) & (d < hi)
        if pick.any():
            out.append((S[pick].mean(), I[pick].mean(), pick.sum(), d[pick].mean()))
        else:
            out.append((np.nan, np.nan, 0, np.nan))
    return out


def brownian_cov(x, kappa, n_events):
    r = np.linalg.norm(x - O, axis=1)
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    k = (kappa ** 2 / (8 * n_events)) * (r[:, None] + r[None, :] - d)
    return np.kron(k, np.eye(4)), k


def main():
    rng = np.random.default_rng(20261126)
    for N, configs in ((128, 12), (256, 12), (512, 8)):
        t0 = time.time()
        acc = {k: [] for k in ("nc", "no", "pairs", "dbar", "far_frac_outer", "s_int", "s90",
                               "Y_iid", "Y_iid_c", "Y_het", "Y_het_w", "Xn_iid", "Xf_iid",
                               "Y_k1", "Y_k03", "absorb_far_k1", "fieldfrac_far_k1", "fieldfrac_far_k03",
                               "beta_c", "beta_cm", "beta_int_far", "mineig", "Y_k1_w",
                               "Y_k2_w", "Xn_k1", "Xf_k1", "Sfar_k1", "Ifar_k1", "Snear_k1", "Inear_k1")}
        for c in range(configs):
            xa = fp.sample(fp.BACKGROUND, N, rng)
            x = xa[:N]
            rel = fp.order_matrix(xa)[:N, :N]
            central = dp.central_mask(x)
            outer = ~central & rel.any(axis=0) & rel.any(axis=1)
            inter = central | outer
            xi = x[inter]
            n = len(xi)
            acc["nc"].append(central.sum()); acc["no"].append(outer.sum())
            s = scale_of(xi)
            acc["s_int"].append(s.mean()); acc["s90"].append(np.mean(s > 0.9))
            d = np.linalg.norm(xi[:, None] - xi[None], axis=2)
            up = np.triu(np.ones((n, n), bool), 1)
            cnt = [np.sum(up & (d >= lo) & (d < hi)) for lo, hi in zip(EDGES[:-1], EDGES[1:])]
            acc["pairs"].append(cnt)
            far = up & (d >= 0.75)
            acc["dbar"].append(d[far].mean() if far.any() else np.nan)
            isout = outer[inter]
            acc["far_frac_outer"].append(np.mean((isout[:, None] | isout[None, :])[far]))
            M = proj_complement(xi)
            # iid isotropic noise on both states: residual cov = 2 sigma^2 (I-P)
            cov = M @ M.T
            b = bin_stats(cov, xi)
            X = [s_ / i_ for s_, i_, _, _ in b]
            acc["Y_iid"].append(X[3] / X[0]); acc["Xn_iid"].append(X[0]); acc["Xf_iid"].append(X[3])
            # central-only fit, central pairs, iid
            # heteroscedastic: outer variance 3x central, events with scale>0.9 6x
            var = np.where(isout, 3.0, 1.0) * np.where(s > 0.9, 2.0, 1.0)
            D = np.kron(np.diag(var), np.eye(4))
            covh = M @ D @ M.T
            b = bin_stats(covh, xi); X = [s_ / i_ for s_, i_, _, _ in b]
            acc["Y_het"].append(X[3] / X[0])
            b = bin_stats(covh, xi, wave=True); X = [s_ / i_ for s_, i_, _, _ in b]
            acc["Y_het_w"].append(X[3] / X[0])
            # Brownian field + noise 10/N rms per event per state (sigma^2 = 25/N^2 per coord, both states)
            noise = 2 * (10.0 / N / 2) ** 2
            for kappa, key in ((1.0, "k1"), (0.3, "k03")):
                K4, k1 = brownian_cov(xi, kappa, N)
                if kappa == 1.0:
                    acc["mineig"].append(np.linalg.eigvalsh(k1).min() / np.linalg.eigvalsh(k1).max())
                tot = M @ (K4 + noise * np.eye(4 * n)) @ M.T
                b = bin_stats(tot, xi); X = [s_ / i_ for s_, i_, _, _ in b]
                acc["Y_" + key].append(X[3] / X[0])
                bf = bin_stats(M @ K4 @ M.T, xi)
                bn = bin_stats(M @ (noise * np.eye(4 * n)) @ M.T, xi)
                acc["fieldfrac_far_" + key].append(bf[3][0] / b[3][0])
                if kappa == 1.0:
                    acc["Xn_k1"].append(X[0]); acc["Xf_k1"].append(X[3])
                    acc["Sfar_k1"].append(b[3][0] * N * N); acc["Ifar_k1"].append(b[3][1] * N * N)
                    acc["Snear_k1"].append(b[0][0] * N * N); acc["Inear_k1"].append(b[0][1] * N * N)
                    pre = bin_stats(K4, xi)
                    acc["absorb_far_k1"].append(1 - bf[3][0] / pre[3][0])
                    bw = bin_stats(tot, xi, wave=True); Xw = [s_ / i_ for s_, i_, _, _ in bw]
                    acc["Y_k1_w"].append(Xw[3] / Xw[0])
                    K4b, _ = brownian_cov(xi, 2.0, N)
                    totb = M @ (K4b + noise * np.eye(4 * n)) @ M.T
                    bw = bin_stats(totb, xi, wave=True); Xw = [s_ / i_ for s_, i_, _, _ in bw]
                    acc["Y_k2_w"].append(Xw[3] / Xw[0])
            # beta: central-only fit with field on central events (conversion of existing c bound)
            xc = x[central]
            Mc = proj_complement(xc)
            Kc, _ = brownian_cov(xc, 1.0, N)
            post = Mc @ Kc @ Mc.T
            nc = len(xc)
            acc["beta_c"].append(np.trace(post) / nc * N * 4 * nc / (4 * nc - 7))
            # wave-only field: tau = tau_m * e_m, with E|tau_m(x)-tau_m(y)|^2 = |x-y|/N
            r = np.linalg.norm(xc - O, axis=1)
            dd = np.linalg.norm(xc[:, None] - xc[None], axis=2)
            km = (1.0 / (2 * N)) * (r[:, None] + r[None, :] - dd)
            e = np.array([1, 1, 0, 0]) / math.sqrt(2)
            Km = np.kron(km, np.outer(e, e))
            postm = Mc @ Km @ Mc.T
            C = postm.reshape(nc, 4, nc, 4)
            wv = np.einsum("a,iaia,a->i", e, C, e)
            acc["beta_cm"].append(wv.mean() * N * 4 * nc / (4 * nc - 7))
        print(f"N={N} ({time.time()-t0:.0f}s, {configs} configs)")
        for k, v in acc.items():
            if len(v) == 0:
                continue
            v = np.array(v, dtype=float)
            if v.ndim == 2:
                print(f"  {k}: mean {np.round(v.mean(axis=0),1)}")
            else:
                print(f"  {k}: mean {np.nanmean(v):.4f} sd {np.nanstd(v, ddof=1):.4f} min {np.nanmin(v):.4f} max {np.nanmax(v):.4f}")


if __name__ == "__main__":
    main()
