# ruff: noqa: E501, UP031
"""The full probe-integrated statistic Sbar(omega) (all Taylor orders, t=1)
in the referee's 1+1 analogue, against the bound (5) of gram_fourier_noise.md
and against the three-term bound used on the scan range.

Part A (tensor-grid quadrature, n up to 128000).  On the scan range
    1 <= omega <= n^(K/(2K+2)),  K = 1, 2, 3,
the draft claims ||S||_2 <= C [1/(omega sqrt n) + 1/n + omega^K n^(-(K+1)/2)].
Here the leading constant is known: the first-order term has L2 norm
sigma_1(omega)/sqrt(n), sigma_1 from deterministic quadrature.  The test is
that the simulated rms equals sigma_1(omega)/sqrt(n) within tolerance on the
whole scan range, for every n, and that the ratio to the three-term bracket
never exceeds max_omega omega*sigma_1(omega) by more than the tolerance.

Part B (exact cell integration, no quadrature, n <= 2000).  Frequencies far
beyond sqrt(n), where the Taylor expansion of the proof is useless.  Reported:
the noise there, and whether the remainder term of (5) is attained.

Sbar is the probe INTEGRAL; the finite-probe term is tested in s3.
Deterministic; prints PASS/FAIL; exits nonzero on any FAIL.
"""

from __future__ import annotations

import math
import sys

import numpy as np
import toy_common as tc

FAILS = []
OMEGAS_A = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0])
RUNS_A = ((2000, 200), (8000, 160), (32000, 120), (128000, 60))
OMEGAS_B = np.array([8.0, 32.0, 128.0, 512.0, 2048.0])
RUNS_B = ((500, 300), (1000, 150), (2000, 60))
SIGMA_U = 0.6  # sd(u_hat - u) ~ SIGMA_U / sqrt(n) on the window


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL") + " " + name + (" :: " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def bracket(omega, n, k):
    return 1 / (omega * math.sqrt(n)) + 1 / n + omega**k / n ** ((k + 1) / 2)


def full_bound(omega, n, k):
    return (
        sum(omega ** (j - 2) * n ** (-j / 2) for j in range(1, k + 1))
        + 1 / n
        + omega**k / n ** ((k + 1) / 2)
    )


def part_a():
    ug, vg, hu, hv = tc.grid(4200, 40)
    pop = tc.population_moments(ug, vg)
    sigma1 = np.sqrt(np.array([tc.first_order_variance(float(w)) for w in OMEGAS_A]))
    peak = float(np.max(OMEGAS_A * sigma1))
    print(
        "A: deterministic omega*sigma_1(omega): "
        + " ".join(f"{x:.5f}" for x in OMEGAS_A * sigma1)
    )
    stats = {}
    for n, reps in RUNS_A:
        rng = np.random.default_rng(5000 + n)
        sq = np.empty((reps, OMEGAS_A.size))
        for r in range(reps):
            g, _ = tc.scaled_scores(rng, n, ug, vg, hu, hv, OMEGAS_A, [1.0], pop)
            sq[r] = np.abs(g[0]) ** 2
        stats[n] = sq
        rms = np.sqrt(sq.mean(axis=0))
        rel = 0.5 * sq.std(axis=0, ddof=1) / math.sqrt(reps) / sq.mean(axis=0)
        print(
            "A: n=%6d rms*sqrt(n)/sigma_1 : " % n
            + " ".join(
                f"{a:.3f}[{b:.2f}]"
                for a, b in zip(rms * math.sqrt(n) / sigma1, rel, strict=False)
            )
        )
    print(
        "A: columns are omega = "
        + " ".join(f"{w:g}" for w in OMEGAS_A)
        + " ; [rel SE of rms]"
    )

    for k in (1, 2, 3):
        worst_dev, worst_ratio, count, ok = 0.0, 0.0, 0, True
        for n, reps in RUNS_A:
            cap = n ** (k / (2 * k + 2))
            sq = stats[n]
            rms = np.sqrt(sq.mean(axis=0))
            rel = 0.5 * sq.std(axis=0, ddof=1) / math.sqrt(reps) / sq.mean(axis=0)
            for j, w in enumerate(OMEGAS_A):
                if w > cap:
                    continue
                count += 1
                dev = rms[j] * math.sqrt(n) / sigma1[j] - 1
                worst_dev = max(worst_dev, abs(dev))
                ok &= abs(dev) < max(4 * rel[j], 0.12)
                ratio = rms[j] / bracket(w, n, k)
                worst_ratio = max(worst_ratio, ratio)
                ok &= ratio < 1.12 * peak or ratio - 4 * rel[j] * ratio < 1.12 * peak
        check(
            "A%d scan range K=%d (%d points): rms = sigma_1(omega)/sqrt(n) and rms/bracket <= 1.12*max(omega sigma_1)"
            % (k, k, count),
            bool(ok),
            f"worst |rms*sqrt(n)/sigma_1 - 1| = {worst_dev:.3f}; sup rms/bracket = {worst_ratio:.5f} vs max(omega sigma_1) = {peak:.5f}",
        )

    # power-law fit where the first-order term dominates and 1/omega has set in
    pts, means, rels, reps_data = [], [], [], []
    for n, reps in RUNS_A:
        for j, w in enumerate(OMEGAS_A):
            if w >= 32 and w * SIGMA_U / math.sqrt(n) <= 0.5:
                sq = stats[n][:, j]
                pts.append((math.log(w), math.log(n)))
                means.append(sq.mean())
                rels.append(sq.std(ddof=1) / math.sqrt(reps) / sq.mean())
                reps_data.append((n, j))
    pts, means, rels = np.array(pts), np.array(means), np.array(rels)

    def fit(values):
        x = np.column_stack([np.ones(len(pts)), pts])
        wts = 1 / rels**2
        beta = np.linalg.solve(x.T @ (wts[:, None] * x), x.T @ (wts * np.log(values)))
        return beta[1] / 2, beta[2] / 2

    point = fit(means)
    rng = np.random.default_rng(77)
    boot = []
    for _ in range(400):
        idx = {n: rng.integers(0, reps, reps) for n, reps in RUNS_A}
        boot.append(fit(np.array([stats[n][idx[n], j].mean() for n, j in reps_data])))
    se = np.std(np.array(boot), axis=0, ddof=1)
    check(
        "A4 full statistic, first-order regime (%d points, 32<=omega, omega*sd(u_hat)<=0.5): n-exponent -1/2"
        % len(pts),
        abs(point[1] + 0.5) < max(3 * se[1], 0.05) and se[1] < 0.05,
        f"L2 norm ~ omega^({point[0]:.3f} +- {se[0]:.3f}) n^({point[1]:.3f} +- {se[1]:.3f})",
    )
    check(
        "A5 full statistic, first-order regime: omega-exponent <= -1",
        point[0] < -1 + max(3 * se[0], 0.05) and se[0] < 0.05,
        f"{point[0]:.3f} +- {se[0]:.3f}",
    )
    return peak


def part_b(peak):
    sigma_inf = math.sqrt(tc.first_order_variance(2048.0)) * 2048.0
    worst_two_term, attained = 0.0, 0.0
    for n, reps in RUNS_B:
        rng = np.random.default_rng(6000 + n)
        sq = np.empty((reps, OMEGAS_B.size))
        for r in range(reps):
            scores, _ = tc.exact_cell_scores(rng.random(n), rng.random(n), OMEGAS_B)
            sq[r] = np.abs(scores) ** 2
        rms = np.sqrt(sq.mean(axis=0))
        rel = 0.5 * sq.std(axis=0, ddof=1) / math.sqrt(reps) / sq.mean(axis=0)
        print(
            "B: n=%5d omega/sqrt(n)     : " % n
            + " ".join("%8.2f" % (w / math.sqrt(n)) for w in OMEGAS_B)
        )
        print(
            "B: n=%5d rms*omega*sqrt(n) : " % n
            + " ".join(f"{x:8.4f}" for x in rms * OMEGAS_B * math.sqrt(n))
            + "   rel SE "
            + " ".join(f"{x:.2f}" for x in rel)
        )
        print(
            "B: n=%5d rms*n             : " % n + " ".join(f"{x:8.4f}" for x in rms * n)
        )
        for kk in (1, 2, 3):
            print(
                "B: n=%5d rms/(5)-bound K=%d : " % (n, kk)
                + " ".join(
                    "%8.1e" % (rms[j] / full_bound(w, n, kk))
                    for j, w in enumerate(OMEGAS_B)
                )
            )
        for j, w in enumerate(OMEGAS_B):
            if math.sqrt(n) < w <= n / 4:
                worst_two_term = max(
                    worst_two_term, rms[j] / (1 / (w * math.sqrt(n)) + 1 / n)
                )
                attained = max(attained, rms[j] / (w**2 / n**1.5))
    print(f"B: sigma_inf = lim omega*sigma_1 = {sigma_inf:.5f}")
    check(
        "B1 beyond the Taylor range (sqrt(n) < omega <= n/4, exact cells): rms <= 2.5*max(omega sigma_1)*[1/(omega sqrt n)+1/n]",
        worst_two_term < 2.5 * peak,
        f"sup rms/[1/(omega sqrt n)+1/n] = {worst_two_term:.5f} vs 2.5*max(omega sigma_1) = {2.5 * peak:.5f}",
    )
    check(
        "B2 the K=2 remainder term omega^2 n^(-3/2) of (5) is NOT attained beyond the Taylor range (ratio < 1e-2)",
        attained < 1e-2,
        f"sup rms/(omega^2 n^-1.5) = {attained:.2e}",
    )


if __name__ == "__main__":
    part_b(part_a())
    print("SUMMARY failures=%d" % len(FAILS))
    sys.exit(1 if FAILS else 0)
