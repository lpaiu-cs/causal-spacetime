# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written in this folder by the fresh-context AI referee of local_null_challenge_r2.md.
"""Referee check r2: Lemmas 2.1-2.3 and the counting in Lemma 4.3.

checks.py N4 tests Lemma 2.1 for means of i.i.d. Bernoulli variables and N5
tests the averaging identity and the MGF comparison on one kernel. This
script tests the statement actually used for the a-moment and its increment:
the data-driven radius r_D of (3), with kk = floor(n/2), applied to a
COMPLETE order-two U-statistic, by exact enumeration of its law.

R takes values in a 4-point set with probabilities pi; U depends only on the
occupation counts, so its exact law is a sum over compositions of n.

  U1  Pr(|U - d| > r_H) <= 2 e^-B (Hoeffding part, kk = floor(n/2)).
  U2  Pr(|U - d| > sqrt(2 sigma^2 B/kk) + 2B/(3kk)) <= 2 e^-B with the TRUE
      kernel variance sigma^2 (Bernstein part).
  U3  Pr(|U - d| > r_D(U)) <= 4 e^-B (Lemma 2.1 conclusion via Lemma 2.2).
  U4  Negative control: the same with kk = C(n, 2) (pairs treated as
      independent) must exceed 4 e^-B somewhere. Information: kk = n.
  K1  Union-bound count: 6P + 6 C(P, 2) <= 6N + 3N^2 and
      (6N + 3N^2) 2^(1-B) <= (3/4) N^-2 with the least B, N = 4..20000.
  K2  Lemma 4.3: floor(floor(N/2)/2) >= N/5 for every N >= 20 (checked to
      10^5); the N below 20 where it fails are printed.
  K3  Bernstein root of Lemma 2.1: exponent = B at r*, r* <= stated bound.
Run: python r2_confidence_bounds.py
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

OK = []


def rec(name, ok, msg):
    OK.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def compositions(n, parts):
    for cuts in itertools.combinations(range(n + parts - 1), parts - 1):
        prev, out = -1, []
        for c in cuts:
            out.append(c - prev - 1)
            prev = c
        out.append(n + parts - 1 - prev - 1)
        yield out


def exact_law(n, pi, H):
    """Support and probabilities of U = C(n,2)^-1 sum_{j<l} H[R_j, R_l]."""
    comps = np.array(list(compositions(n, len(pi))), dtype=float)
    lg = math.lgamma(n + 1) - np.sum([np.vectorize(math.lgamma)(c + 1) for c in comps.T], axis=0)
    lg = lg + comps @ np.log(pi)
    prob = np.exp(lg)
    pairs = np.zeros(len(comps))
    S = len(pi)
    for i in range(S):
        pairs += comps[:, i] * (comps[:, i] - 1) / 2 * H[i, i]
        for j in range(i + 1, S):
            pairs += comps[:, i] * comps[:, j] * H[i, j]
    U = pairs / math.comb(n, 2)
    return U, prob


def main():
    # kernels: symmetric {0,1}; the first is a nested-increment kernel of a toy
    # poset (states: 0 = below X, 1 = below Y only, 2 = incomparable, 3 = other)
    kernels = {
        "increment": np.array([[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]),
        "ancestor-pair": np.array([[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]),
        "dense": np.array([[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1]]),
        "star": np.array([[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]]),
    }
    pis = [
        np.array([0.05, 0.01, 0.5, 0.44]),
        np.array([0.3, 0.05, 0.3, 0.35]),
        np.array([0.6, 0.2, 0.1, 0.1]),
        np.array([0.15, 0.003, 0.4, 0.447]),
    ]
    ns = (6, 9, 14, 21, 40, 80)
    Bs = (1.0, 2.0, 3.0, 5.0)
    w1 = w2 = w3 = 0.0
    neg_pairs = neg_n = 0.0
    for (name, H), pi, n in itertools.product(kernels.items(), pis, ns):
        U, prob = exact_law(n, pi, H)
        d = float(pi @ H @ pi)
        sig2 = d - d * d
        for B in Bs:
            kk = n // 2
            rH = math.sqrt(B / (2 * kk))
            dU = np.minimum(1.0, U + rH)
            rD = np.minimum(rH, np.sqrt(2 * dU * B / kk) + 2 * B / (3 * kk))
            dev = np.abs(U - d)
            f1 = prob[dev > rH + 1e-15].sum() / (2 * math.exp(-B))
            rb = math.sqrt(2 * sig2 * B / kk) + 2 * B / (3 * kk)
            f2 = prob[dev > rb + 1e-15].sum() / (2 * math.exp(-B))
            f3 = prob[dev > rD + 1e-15].sum() / (4 * math.exp(-B))
            w1, w2, w3 = max(w1, f1), max(w2, f2), max(w3, f3)
            for kk_bad, tag in ((math.comb(n, 2), "pairs"), (n, "n")):
                rHb = math.sqrt(B / (2 * kk_bad))
                dUb = np.minimum(1.0, U + rHb)
                rDb = np.minimum(rHb, np.sqrt(2 * dUb * B / kk_bad) + 2 * B / (3 * kk_bad))
                fb = prob[dev > rDb + 1e-15].sum() / (4 * math.exp(-B))
                if tag == "pairs":
                    neg_pairs = max(neg_pairs, fb)
                else:
                    neg_n = max(neg_n, fb)
    n_cases = len(kernels) * len(pis) * len(ns) * len(Bs)
    rec("U1 Hoeffding for U, kk = floor(n/2)", w1 <= 1, f"{n_cases} cases; worst / 2e^-B = {w1:.3f}")
    rec("U2 Bernstein for U, true variance", w2 <= 1, f"worst / 2e^-B = {w2:.3f}")
    rec("U3 radius r_D of (3) for U", w3 <= 1, f"worst / 4e^-B = {w3:.3f}")
    rec("U4 negative control kk = C(n,2) violates", neg_pairs > 1, f"worst / 4e^-B = {neg_pairs:.2f}")
    print(f"[INFO] U4 kk = n (information): worst / 4e^-B = {neg_n:.3f}", flush=True)

    # K1 union bound counting
    ok_cnt = ok_union = True
    for N in range(4, 20001):
        P = N - N // 2
        ok_cnt &= 6 * P + 6 * math.comb(P, 2) <= 6 * N + 3 * N * N
        B = (32 * N**4 - 1).bit_length()
        ok_union &= (6 * N + 3 * N * N) * 2.0 ** (1 - B) <= 0.75 / N**2 * (1 + 1e-12)
    rec("K1 count 6P + 6C(P,2) <= 6N + 3N^2", ok_cnt, "N = 4..20000")
    rec("K1 union bound <= (3/4) N^-2 with least B", ok_union, "N = 4..20000")

    # K2 kk >= N/5
    bad = [N for N in range(4, 100001) if (N // 2) // 2 < N / 5]
    rec(
        "K2 floor(floor(N/2)/2) >= N/5 for N >= 20",
        all(N < 20 for N in bad),
        f"fails only at N = {bad}",
    )

    # K3 Bernstein root algebra
    worst_root, worst_exp = 0.0, 0.0
    for kk, B, s2 in itertools.product((1, 3, 10, 100, 10**4), (0.5, 3, 30), (0, 1e-6, 1e-3, 0.1, 0.25)):
        r = B / (3 * kk) + math.sqrt(B * B / (9 * kk * kk) + 2 * s2 * B / kk)
        expo = kk * r * r / (2 * s2 + 2 * r / 3)
        worst_exp = max(worst_exp, abs(expo - B) / B)
        worst_root = max(worst_root, r / (math.sqrt(2 * s2 * B / kk) + 2 * B / (3 * kk)))
    rec(
        "K3 Bernstein root",
        worst_exp < 1e-12 and worst_root <= 1 + 1e-12,
        f"max |exponent - B|/B = {worst_exp:.1e}; max r*/bound = {worst_root:.4f}",
    )
    n_fail = OK.count(False)
    print(f"\n{OK.count(True)} PASS, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
