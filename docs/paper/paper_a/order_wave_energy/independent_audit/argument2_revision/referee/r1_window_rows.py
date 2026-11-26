# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written in this folder by the fresh-context AI referee of local_null_challenge_r2.md.
"""Referee check r1: the calibration window and the 12 frozen witnesses.

Independent of checks.py: p0 is computed here from the background interval
kernel of [GFL] Remark K1 by Gauss-Legendre quadrature,

    p0(q l^2) = Lambda_2,0(0, l) / (3 Lambda_0(0, l)^2),
    Lambda_0(s, t)   = int_s^t K(s, w, t) dw,
    Lambda_2,0(s, t) = int_s^t Lambda_0(w, t) K(s, w, t) dw,
    K(s, w, t) = a(w)^2 m(s, w) m(w, t) / m(s, t),

with a = cos(k u), m = int a^-2, not from the 0F1 series that checks.py and
the implementation use. Lambda_0(w, t) is itself computed by quadrature, so
the closed form alpha(l) is checked, not assumed.

Claims tested (local_null_challenge_r2.md):
  A  p0(1) = 0.095919, p0(0.25) = 0.099031; widths 4.1% and 1.0% (Remark 2.6).
  B  Lambda_0 = alpha(l) and the flat limit p0(0) = 1/10.
  C  p0 decreasing on a 201-point grid of [0, 1] (not a certificate).
  D  Section 6 table: limit ranges [a_X/p_Y^2, a_Y/p_X^2] to 1e-6, the
     skipped / evaluated verdicts, in particular all six q = 0.25 rows
     (the audit record's item 4 lists only r = 1024, 2048 there).
  E  Which rows have a positive population decrease; how many of those are
     evaluable ("at most two").
Run: python r1_window_rows.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROWS = HERE.parents[2] / "local_null_evidence" / "rows.jsonl"
X, W = np.polynomial.legendre.leggauss(40)
OK = []


def rec(name, ok, msg):
    OK.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def gl(a, b, n_panels=8):
    e = np.linspace(a, b, n_panels + 1)
    lo, hi = e[:-1, None], e[1:, None]
    x = (0.5 * (hi - lo) * X[None, :] + 0.5 * (hi + lo)).ravel()
    w = (0.5 * (hi - lo) * W[None, :]).ravel()
    return x, w


def kernel(k, s, w, t):
    def m(a, b):
        return (np.tan(k * b) - np.tan(k * a)) / k

    return np.cos(k * w) ** 2 * m(s, w) * m(w, t) / m(s, t)


def lam0(k, s, t):
    x, w = gl(s, t)
    return float(np.sum(w * kernel(k, s, x, t)))


def p0_quad(z):
    """p0(z) with l = 1 and q = z (k = sqrt z); flat case handled by k -> tiny."""
    k = max(math.sqrt(z), 1e-4)
    s, t = -0.5, 0.5  # any unit segment; K is translation invariant
    x, w = gl(s, t)
    inner = np.array([lam0(k, xi, t) for xi in x])
    l2 = float(np.sum(w * inner * kernel(k, s, x, t)))
    l0 = lam0(k, s, t)
    return l2 / (3 * l0**2)


def alpha(k, ell):
    q = k * k
    x = k * ell
    return (1 - x / math.tan(x)) / (2 * q)


def main():
    # B: Lambda_0 = alpha, translation invariance
    errs = []
    for q in (0.25, 1.0):
        k = math.sqrt(q)
        for s, t in ((-0.5, 0.5), (-0.5, -0.1), (0.05, 0.45), (-0.3, 0.3)):
            errs.append(abs(lam0(k, s, t) - alpha(k, t - s)) / alpha(k, t - s))
    rec("B Lambda_0 = alpha(l)", max(errs) < 1e-12, f"max rel err {max(errs):.1e}")
    p_flat = p0_quad(1e-8)
    rec("B flat limit p0(0) = 1/10", abs(p_flat - 0.1) < 1e-8, f"p0(1e-8) = {p_flat:.12f}")

    # A: window
    p1, p25 = p0_quad(1.0), p0_quad(0.25)
    w1, w25 = (0.1 - p1) / 0.1, (0.1 - p25) / 0.1
    rec(
        "A p0(1), p0(0.25)",
        abs(p1 - 0.095919) < 1e-6 and abs(p25 - 0.099031) < 1e-6,
        f"p0(1) = {p1:.8f}, p0(0.25) = {p25:.8f}; widths {100 * w1:.3f}%, {100 * w25:.3f}%",
    )
    # C: monotone on a grid
    zs = np.linspace(0, 1, 201)
    vals = np.array([p0_quad(max(z, 1e-8)) for z in zs])
    d = np.diff(vals)
    dl = np.diff(np.log(vals)) / np.diff(zs)
    rec(
        "C p0 decreasing on grid",
        bool(np.all(d < 0)),
        f"max step diff {d.max():.2e}; max d log p0/dz (finite diff) {dl.max():.4f}",
    )

    # D, E: rows
    rows = [json.loads(s) for s in ROWS.read_text().splitlines() if s.strip()]
    paper = {
        (0.25, 64): (0.074560, 0.132767), (0.25, 128): (0.085922, 0.115252),
        (0.25, 256): (0.092414, 0.107175), (0.25, 512): (0.095890, 0.103299),
        (0.25, 1024): (0.097688, 0.101401), (0.25, 2048): (0.098604, 0.100462),
        (1.0, 64): (0.072984, 0.131380), (1.0, 128): (0.084392, 0.113799),
        (1.0, 256): (0.090916, 0.105712), (1.0, 512): (0.094410, 0.101836),
        (1.0, 1024): (0.096219, 0.099940), (1.0, 2048): (0.097139, 0.099002),
    }
    win = {0.25: p25, 1.0: p1}
    worst, evaluable, positive, nested = 0.0, [], [], True
    for row in rows:
        (pX, fX, cX), (pY, fY, cY) = row["moments"]
        aX, aY = pX * pX * cX, pY * pY * cY
        nested &= pY > pX and fY < fX and aY > aX
        lo, hi = aX / pY**2, aY / pX**2
        key = (row["q"], row["r"])
        worst = max(worst, abs(lo - paper[key][0]), abs(hi - paper[key][1]))
        inside = win[row["q"]] < lo and hi < 0.1
        # margin by which the limit range misses / fits the window
        miss = min(lo - win[row["q"]], 0.1 - hi)
        dec = row["population_v_decrease"]
        if inside:
            evaluable.append(key)
        if dec > 0:
            positive.append(key)
        print(
            f"[INFO] row q={row['q']}, r={row['r']}: [{lo:.6f}, {hi:.6f}] "
            f"window ({win[row['q']]:.6f}, 0.1) -> {'evaluated' if inside else 'SKIPPED'}"
            f" (signed margin {miss:+.2e}); decrease {dec:+.3e}",
            flush=True,
        )
    rec("D nested signs, strict", nested, "p up, f down, a up in all 12 rows")
    rec("D table of Section 6 to 1e-6", worst < 1e-6, f"max |diff| = {worst:.2e}")
    rec(
        "D only q=1, r=1024 and 2048 evaluable",
        sorted(evaluable) == [(1.0, 1024), (1.0, 2048)],
        f"evaluable: {sorted(evaluable)}",
    )
    q25_small = [(0.25, r) for r in (64, 128, 256, 512)]
    rec(
        "D q=0.25 rows with r <= 512 are skipped (missing from AUDIT item 4)",
        all(k not in evaluable for k in q25_small),
        "all four skipped",
    )
    rec(
        "E ten positive decreases, two of them evaluable",
        len(positive) == 10 and set(evaluable) <= set(positive),
        f"positive: {len(positive)}; non-positive: "
        f"{sorted(set(paper) - set(positive))}",
    )
    n_fail = OK.count(False)
    print(f"\n{OK.count(True)} PASS, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
