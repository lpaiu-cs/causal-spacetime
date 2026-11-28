# ruff: noqa: E501, UP031
"""Direct test of the order-by-order bound (4) of gram_fourier_noise.md,
    || k-th Taylor term of the probe-integrated score ||_2 <= C_k omega^(k-2) n^(-k/2),
in the referee's 1+1 analogue (toy_common.py).

Method.  For one set of n shared references let eps(X) be the centred
empirical error of the raw moments.  G(t) = int F(mu_0(X) + t*eps(X)) dX is
smooth in t, G(1) is the statistic, and its Taylor coefficients a_k at t=0
are exactly the probe-integrated order-k Taylor terms of the proof.  They are
extracted from G(+-tau), G(+-2 tau) by Richardson combination (tau = 1/4):
    a_1 = [8 O(tau) - O(2tau)]/(6 tau),   a_3 = [O(2tau) - 2 O(tau)]/(6 tau^3),
    a_2 = [16 E(tau) - E(2tau)]/(12 tau^2),
with O, E the odd and even parts.  Truncation errors are -4 tau^4 a_5,
5 tau^2 a_5 and -4 tau^4 a_6.

Independent cross-check of k=1: n*E|a_1|^2 must equal the variance of the
influence function, which toy_common.first_order_variance computes by
deterministic quadrature with no random numbers.

Acceptance thresholds were fixed after exploratory runs with other seeds
(see report.md); they are tolerances on exponents, not fitted constants.

Deterministic; prints PASS/FAIL; exits nonzero on any FAIL.
"""

from __future__ import annotations

import math
import sys

import numpy as np
import toy_common as tc

FAILS = []
TAU = 0.25
OMEGAS = np.array([16.0, 32.0, 64.0, 128.0, 256.0])
RUNS = ((4000, 240), (16000, 240), (64000, 240))
SIGMA_U = 0.6  # sd(u_hat - u) ~ SIGMA_U / sqrt(n) on the window (referee estimate)
VALID = 1.5  # use points with omega * sd(u_hat) <= VALID, where truncation is < few %


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL") + " " + name + (" :: " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def collect():
    ug, vg, hu, hv = tc.grid(4200, 24)
    pop = tc.population_moments(ug, vg)
    ts = [1.0, TAU, -TAU, 2 * TAU, -2 * TAU]
    data = {}
    for n, reps in RUNS:
        rng = np.random.default_rng(4000 + n)
        rows = np.empty((reps, 5, OMEGAS.size))
        for r in range(reps):
            g, _ = tc.scaled_scores(rng, n, ug, vg, hu, hv, OMEGAS, ts, pop)
            o1, o2 = 0.5 * (g[1] - g[2]), 0.5 * (g[3] - g[4])
            e1, e2 = 0.5 * (g[1] + g[2]), 0.5 * (g[3] + g[4])
            a1 = (8 * o1 - o2) / (6 * TAU)
            a3 = (o2 - 2 * o1) / (6 * TAU**3)
            a2 = (16 * e1 - e2) / (12 * TAU**2)
            rows[r] = np.abs(np.array([g[0], a1, a2, a3, g[0] - a1 - a2 - a3])) ** 2
        data[n] = rows
    return data


def fit(means, rel_se, mask):
    """Weighted LS of log E|a|^2 on (1, log omega, log n) over the masked points."""
    ns = np.array([n for n, _ in RUNS], dtype=float)
    lw, ln = np.meshgrid(np.log(OMEGAS), np.log(ns), indexing="xy")
    x = np.stack([np.ones(mask.sum()), lw[mask], ln[mask]], axis=1)
    wts = 1 / rel_se[mask] ** 2
    y = np.log(means[mask])
    beta = np.linalg.solve(x.T @ (wts[:, None] * x), x.T @ (wts * y))
    return beta[1], beta[2]


def fit_fixed_effects(means, rel_se, mask):
    """Each exponent with the other variable absorbed by free intercepts:
    log y = c_omega + e_n log n   and   log y = c_n + e_omega log omega.
    This cannot leak a non-power-law omega dependence into the n exponent."""
    ns = np.array([n for n, _ in RUNS], dtype=float)
    out = []
    for axis in (
        0,
        1,
    ):  # axis 0: slope in log n per omega column; axis 1: slope in log omega per n row
        num = den = 0.0
        for idx in range(mask.shape[1 - axis]):
            sel = mask[:, idx] if axis == 0 else mask[idx, :]
            if sel.sum() < 2:
                continue
            xs = (np.log(ns) if axis == 0 else np.log(OMEGAS))[sel]
            ys = np.log(means[:, idx] if axis == 0 else means[idx, :])[sel]
            wts = 1 / (rel_se[:, idx] if axis == 0 else rel_se[idx, :])[sel] ** 2
            xbar = np.sum(wts * xs) / np.sum(wts)
            num += np.sum(wts * (xs - xbar) * ys)
            den += np.sum(wts * (xs - xbar) ** 2)
        out.append(num / den)
    return out[1], out[0]  # (omega exponent, n exponent) of the mean square


def exponents_fixed(data, k, mask, draws=400):
    reps = [data[n][:, k, :] for n, _ in RUNS]
    means = np.array([r.mean(axis=0) for r in reps])
    rel_se = (
        np.array([r.std(axis=0, ddof=1) / math.sqrt(r.shape[0]) for r in reps]) / means
    )
    point = fit_fixed_effects(means, rel_se, mask)
    rng = np.random.default_rng(1999 + k)
    boot = []
    for _ in range(draws):
        res = np.array(
            [r[rng.integers(0, r.shape[0], r.shape[0])].mean(axis=0) for r in reps]
        )
        boot.append(fit_fixed_effects(res, rel_se, mask))
    se = np.std(np.array(boot), axis=0, ddof=1)
    return point[0] / 2, se[0] / 2, point[1] / 2, se[1] / 2


def exponents(data, k, mask, draws=400):
    reps = [data[n][:, k, :] for n, _ in RUNS]
    means = np.array([r.mean(axis=0) for r in reps])
    rel_se = (
        np.array([r.std(axis=0, ddof=1) / math.sqrt(r.shape[0]) for r in reps]) / means
    )
    point = fit(means, rel_se, mask)
    rng = np.random.default_rng(999 + k)
    boot = []
    for _ in range(draws):
        res = np.array(
            [r[rng.integers(0, r.shape[0], r.shape[0])].mean(axis=0) for r in reps]
        )
        boot.append(fit(res, rel_se, mask))
    se = np.std(np.array(boot), axis=0, ddof=1)
    # exponents of the L2 norm are half those of the mean square
    return point[0] / 2, se[0] / 2, point[1] / 2, se[1] / 2, means, rel_se


def main():
    data = collect()
    ns = np.array([n for n, _ in RUNS], dtype=float)
    valid = (OMEGAS[None, :] * SIGMA_U / np.sqrt(ns)[:, None]) <= VALID
    print(
        f"valid mask (rows n={[int(n) for n in ns]}, cols omega={[int(w) for w in OMEGAS]}):"
    )
    print(valid.astype(int))

    # --- deterministic first-order variance and its asymptote -------------
    sigma1 = np.array([tc.first_order_variance(float(w)) for w in OMEGAS])
    tail = np.array(
        [tc.first_order_variance(float(w)) * w**2 for w in (512.0, 1024.0, 2048.0)]
    )
    print(
        "deterministic sigma_1(omega)*omega at omega=16..256 : "
        + " ".join(
            f"{math.sqrt(s * w * w):.5f}" for s, w in zip(sigma1, OMEGAS, strict=False)
        )
    )
    print(
        "deterministic sigma_1(omega)*omega at 512,1024,2048 : "
        + " ".join(f"{math.sqrt(t):.5f}" for t in tail)
    )
    check(
        "T0 first-order term is exactly C/omega at large omega (deterministic)",
        abs(tail[2] / tail[1] - 1) < 0.01 and abs(tail[1] / tail[0] - 1) < 0.01,
        f"omega^2 sigma_1^2 changes by {tail[1] / tail[0] - 1:.2e}, {tail[2] / tail[1] - 1:.2e} between octaves; sigma_inf = {math.sqrt(tail[2]):.5f}",
    )

    names = {1: "a_1", 2: "a_2", 3: "a_3"}
    norm = {
        1: lambda n: OMEGAS * math.sqrt(n),
        2: lambda n: np.full(OMEGAS.size, float(n)),
        3: lambda n: n**1.5 / OMEGAS,
    }
    label = {1: "omega*sqrt(n)", 2: "n", 3: "n^1.5/omega"}
    results = {}
    for k in (1, 2, 3):
        ew, sew, en, sen, means, rel_se = exponents(data, k, valid)
        results[k] = (ew, sew, en, sen)
        print(
            f"rms |{names[k]}| * {label[k]}  (rows n, cols omega; rel SE of rms in brackets)"
        )
        for i, n in enumerate(ns):
            print(
                "  n=%6d " % n
                + " ".join(
                    "{:.5f}[{:.2f}]{}".format(
                        math.sqrt(means[i, j]) * norm[k](n)[j],
                        rel_se[i, j] / 2,
                        "" if valid[i, j] else "*",
                    )
                    for j in range(OMEGAS.size)
                )
            )
        print(
            "  fitted L2 exponents for %s: omega^(%.3f +- %.3f)  n^(%.3f +- %.3f)   [claim (4): omega^(<=%d) n^(%.1f)]"
            % (names[k], ew, sew, en, sen, k - 2, -k / 2)
        )

    # simulation vs deterministic first-order variance
    means1 = np.array([data[n][:, 1, :].mean(axis=0) for n, _ in RUNS])
    se1 = np.array(
        [
            data[n][:, 1, :].std(axis=0, ddof=1) / math.sqrt(data[n].shape[0])
            for n, _ in RUNS
        ]
    )
    ratio = means1 * ns[:, None] / sigma1[None, :]
    zmax = float(np.max(np.abs(ratio - 1) / (se1 * ns[:, None] / sigma1[None, :])))
    check(
        "T1 simulated n*E|a_1|^2 equals the deterministic influence-function variance",
        zmax < 4.0 and float(np.max(np.abs(ratio - 1))) < 0.35,
        "ratio range %.3f..%.3f, worst |z| = %.2f over %d points"
        % (ratio.min(), ratio.max(), zmax, ratio.size),
    )

    ew, sew, en, sen = results[1]
    check(
        "T2 k=1 n-exponent = -1/2",
        abs(en + 0.5) < max(3 * sen, 0.05) and sen < 0.05,
        f"{en:.3f} +- {sen:.3f}",
    )
    check(
        "T2 k=1 omega-exponent <= -1 (bound (3))",
        ew < -1 + max(3 * sew, 0.05) and sew < 0.05,
        f"{ew:.3f} +- {sew:.3f}",
    )
    ew, sew, en, sen = results[2]
    check(
        "T3 k=2 n-exponent = -1",
        abs(en + 1.0) < max(3 * sen, 0.1) and sen < 0.1,
        f"{en:.3f} +- {sen:.3f}",
    )
    check(
        "T3 k=2 omega-exponent <= 0",
        ew < 0 + max(3 * sew, 0.1) and sew < 0.1,
        f"{ew:.3f} +- {sew:.3f}",
    )
    ew, sew, en, sen = results[3]
    check(
        "T4 k=3 n-exponent = -3/2",
        abs(en + 1.5) < max(3 * sen, 0.1) and sen < 0.1,
        f"{en:.3f} +- {sen:.3f}",
    )
    check(
        "T4 k=3 omega-exponent <= +1",
        ew < 1 + max(3 * sew, 0.1) and sew < 0.1,
        f"{ew:.3f} +- {sew:.3f}",
    )

    # consistency of the extraction: G(1) = a_1 + a_2 + a_3 + higher orders
    res = np.array([data[n][:, 4, :].mean(axis=0) for n, _ in RUNS])
    tot = np.array([data[n][:, 0, :].mean(axis=0) for n, _ in RUNS])
    rel = np.sqrt(res / tot)
    print("rms|G(1)-a_1-a_2-a_3| / rms|G(1)| (rows n, cols omega):")
    for i, n in enumerate(ns):
        print("  n=%6d " % n + " ".join(f"{rel[i, j]:.4f}" for j in range(OMEGAS.size)))
    small = (OMEGAS[None, :] * SIGMA_U / np.sqrt(ns)[:, None]) <= 0.2
    check(
        "T5 statistic = a_1+a_2+a_3 up to <2%% where omega*sd(u_hat) <= 0.2 (%d points)"
        % small.sum(),
        small.sum() >= 3 and float(rel[small].max()) < 0.02,
        f"max relative residual {float(rel[small].max()):.4f}",
    )
    # remainder after K=3 orders: the proof bounds it by C omega^3 n^(-2)
    low = (OMEGAS[None, :] * SIGMA_U / np.sqrt(ns)[:, None]) <= 0.5
    ew, sew, en, sen, _, _ = exponents(data, 4, low)
    print(
        "  fitted L2 exponents for the K=3 remainder G(1)-a_1-a_2-a_3 (%d points with omega*sd(u_hat)<=0.5):"
        " omega^(%.3f +- %.3f) n^(%.3f +- %.3f)   [claim: omega^(<=3) n^(-2)]"
        % (low.sum(), ew, sew, en, sen)
    )
    check(
        "T6 K=3 remainder n-exponent = -(K+1)/2 = -2",
        abs(en + 2.0) < max(3 * sen, 0.1) and sen < 0.1,
        f"{en:.3f} +- {sen:.3f}",
    )
    check(
        "T6 K=3 remainder omega-exponent <= K = 3",
        ew < 3 + max(3 * sew, 0.1) and sew < 0.1,
        f"{ew:.3f} +- {sew:.3f}",
    )
    # secondary, criterion-free report: exponents with free intercepts (added after
    # the first run, because a single power law in omega is only asymptotic here)
    print(
        "INFO fixed-effects fits (n-slope with one intercept per omega; omega-slope with one intercept per n):"
    )
    for kk, name, m in (
        (1, "a_1", valid),
        (2, "a_2", valid),
        (3, "a_3", valid),
        (4, "K=3 remainder", low),
        (0, "full statistic G(1)", small),
    ):
        fw, sfw, fn, sfn = exponents_fixed(data, kk, m)
        print(
            "INFO   %-20s omega^(%.3f +- %.3f)  n^(%.3f +- %.3f)"
            % (name, fw, sfw, fn, sfn)
        )
    print("SUMMARY failures=%d" % len(FAILS))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
