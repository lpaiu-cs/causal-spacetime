# ruff: noqa: E501, UP031
"""Cross-check of the AUTHOR's certified decision routine against the referee's
independent population moments and derivatives.

This is the only script that imports the author's source
(experiments/positive_control/order_wave_local_null.py, read-only).  It asks whether
the explicit test of Section 2, as implemented, agrees with the referee's own
calibration F and with the referee's own exact wave moments.

X1  author's certified gradient/Hessian boxes contain the referee's finite differences.
X2  author's difference_interval, fed the referee's exact moments of the q=1, r=2048
    witness with a synthetic reference count 10^30, has a negative upper endpoint and
    reproduces the archived post-run interval.
X3  background pairs (true difference >= 0): the certified interval never excludes
    the truth (no false rejection on exact background moments).
X4  the rectangular domain R forces a skip of the frozen q=0.25 witnesses (r <= 2048)
    at every reference count; the q=1, r=2048 witness is evaluated.
INFO  smallest synthetic count (decades) at which the witness is certified.

Synthetic counts are population moments rounded to a denominator; no order was
sampled and nothing here is a power result.

Deterministic.  Prints PASS/FAIL lines; exits nonzero on any FAIL.
Run (needs scipy because the author's module imports it):
  E:/lab/causal-spacetime/.venv/Scripts/python.exe s6_author_routine_crosscheck.py
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments" / "positive_control"))

import order_wave_local_null as author  # noqa: E402
from audit_geometry import Wave, alpha0, calibration_f, p0, sym  # noqa: E402
from order_wave_null_coordinate import Box  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def rounded_moments(theta, count):
    dens = (count, count, math.comb(count, 2))
    return tuple(
        F(round(F(float(v)) * d), d) for v, d in zip(theta, dens, strict=False)
    )


def certified(theta_x, theta_y, count):
    exponent = (32 * (2 * count) ** 4).bit_length()
    return author.difference_interval(
        rounded_moments(theta_x, count),
        rounded_moments(theta_y, count),
        count,
        F(1),
        exponent,
        128,
    )


def witness(q=1.0, r=2048):
    w = Wave(q, 1.0 / r)
    u1, u2 = 0.2, 0.2 + math.pi / r
    _, vec = np.linalg.eigh(sym(w.f_matrix_complex(u1) - w.f_matrix_complex(u2)))
    z = 0.2 * vec[:, -1]
    first = w.from_canonical(u1, 0.5, z)
    second = w.from_canonical(u2, 0.5, z)
    second = w.from_canonical(u2, 0.5 - w.margin(first, second) / 2, z)
    gamma = w.gamma()
    return w.moments(first, gamma=gamma), w.moments(second, gamma=gamma)


def part_x1(theta):
    q = 1.0

    def fun(t):
        return calibration_f(t[0], t[1], t[2], q)

    theta = np.array(theta)
    step = 2e-4 * theta
    grad = np.zeros(3)
    hess = np.zeros((3, 3))
    for i in range(3):
        e = np.zeros(3)
        e[i] = step[i]
        grad[i] = (fun(theta + e) - fun(theta - e)) / (2 * step[i])
        for j in range(3):
            g = np.zeros(3)
            g[j] = step[j]
            hess[i, j] = (
                fun(theta + e + g)
                - fun(theta + e - g)
                - fun(theta - e + g)
                + fun(theta - e - g)
            ) / (4 * step[i] * step[j])
    exact = [F(float(v)) for v in theta]
    radii = [F(1, 10**9), F(1, 10**9), F(1, 10**12)]
    boxes = [Box(v - s, v + s) for v, s in zip(exact, radii, strict=False)]
    gbox, hbox = author.derivative_boxes(boxes, F(1), 128)

    def inside(box, value):
        slack = 2e-3 * abs(value) + 1e-9
        return float(box.lo) - slack <= value <= float(box.hi) + slack

    ok_g = all(inside(gbox[i], grad[i]) for i in range(3))
    ok_h = all(inside(hbox[i][j], hess[i, j]) for i in range(3) for j in range(3))
    check(
        "X1 author's certified derivative boxes contain the referee's finite differences",
        ok_g and ok_h,
        "referee grad "
        + np.array2string(grad, precision=4)
        + f"; largest |Hessian| entry {np.abs(hess).max():.4e} "
        + f"(author box [{float(hbox[2][2].lo):.4e}, {float(hbox[2][2].hi):.4e}])",
    )


def part_x2(theta_x, theta_y):
    bound = certified(theta_x, theta_y, 10**30)
    archived = (-3.87906123997e-8, -3.84986837602e-8)
    ok = (
        bound is not None
        and bound.hi < 0
        and abs(float(bound.lo) / archived[0] - 1) < 1e-4
        and abs(float(bound.hi) / archived[1] - 1) < 1e-4
    )
    check(
        "X2 author's routine certifies the witness from the referee's moments (count 10^30)",
        ok,
        "interval "
        + (
            "None"
            if bound is None
            else f"[{float(bound.lo):.9e}, {float(bound.hi):.9e}]"
        )
        + f"; archived post-run [{archived[0]:.9e}, {archived[1]:.9e}]",
    )
    reached = None
    report = []
    for expo in (24, 25, 26):
        b = certified(theta_x, theta_y, 10**expo)
        negative = b is not None and b.hi < 0
        report.append(
            f"10^{expo}: upper {'None' if b is None else f'{float(b.hi):+.2e}'}"
        )
        if negative and reached is None:
            reached = expo
    print(
        "INFO synthetic reference count needed by the implemented test for the r=2048 "
        "witness: " + "; ".join(report) + f" (lambda^-3 = {2048.0**3:.2e})"
    )


def background_moments(q, u, v, rad2):
    k = math.sqrt(q)
    big = 2 * math.tan(k / 2) / k
    s = u + 0.5
    m = (math.tan(k * u) + math.tan(k / 2)) / k
    eta, zeta = v - rad2 / (2 * m), 1 - v - rad2 / (2 * (big - m))
    gamma = alpha0(1.0, q)
    p = alpha0(s, q) * eta**2 / gamma
    f = alpha0(1 - s, q) * zeta**2 / gamma
    return np.array([p, f, p * p * p0(q * s * s)])


def part_x3():
    q = 1.0
    cases = [
        # (u_X, u_Y, v_Y - v_X): same transverse position, so X < Y iff v_Y >= v_X
        (0.2, 0.2 + math.pi / 2048, 0.0),
        (0.2, 0.2 + math.pi / 2048, 2.0e-8),
        (0.0, 5.0e-4, 0.0),
        (-0.1, -0.1 + 2.0e-4, 1.0e-6),
        # A wider pair: the rectangular domain R leaves the calibrated chain range
        # (p0(q), 1/10) and the routine must SKIP, never reject.
        (0.0, 0.01, 0.0),
    ]
    ok, report, evaluated, skipped = True, [], 0, 0
    for u_x, u_y, gap in cases:
        th_x = background_moments(q, u_x, 0.5, 0.04)
        th_y = background_moments(q, u_y, 0.5 + gap, 0.04)
        for expo in (20, 30):
            b = certified(th_x, th_y, 10**expo)
            if b is None:
                skipped += 1
                report.append(f"du {u_y - u_x:.1e}, gap {gap:g}, 10^{expo}: skipped")
                continue
            evaluated += 1
            covers = float(b.lo) <= gap <= float(b.hi)
            ok &= covers and b.hi >= 0
            report.append(
                f"du {u_y - u_x:.1e}, gap {gap:g}, 10^{expo}: [{float(b.lo):+.2e}, {float(b.hi):+.2e}]"
            )
    check(
        "X3 exact background moments: certified interval covers the true v difference >= 0",
        ok and evaluated == 8 and skipped == 2,
        "; ".join(report),
    )


def part_x4(theta_one):
    """The rectangular domain R of Section 2 pairs a from one endpoint with p from the
    other, so the chain coordinate a/p^2 sweeps [a_X/p_Y^2, a_Y/p_X^2].  If that range
    leaves (p0(q), 1/10) the specified test skips the pair for EVERY reference count."""
    rows, ok = [], True
    cases = [(F(1, 4), 0.25, 1024), (F(1, 4), 0.25, 2048), (F(1), 1.0, 2048)]
    for q_exact, q, r in cases:
        th_x, th_y = theta_one if (q, r) == (1.0, 2048) else witness(q, r)
        corner = (th_x[2] / th_y[0] ** 2, th_y[2] / th_x[0] ** 2)
        predicted_skip = corner[0] <= p0(q) or corner[1] >= 0.1
        count = 10**40
        bound = author.difference_interval(
            rounded_moments(th_x, count),
            rounded_moments(th_y, count),
            count,
            q_exact,
            (32 * (2 * count) ** 4).bit_length(),
            128,
        )
        ok &= (bound is None) == predicted_skip
        rows.append(
            f"q={q}, r={r}: corner chain range [{corner[0]:.6f}, {corner[1]:.6f}] vs "
            f"({p0(q):.6f}, 0.1) -> {'SKIPPED' if bound is None else 'evaluated'}"
        )
    check(
        "X4 skip behaviour of the rectangular domain matches the corner-range prediction",
        ok and "SKIPPED" in rows[0] and "SKIPPED" in rows[1] and "evaluated" in rows[2],
        "; ".join(rows) + " (synthetic count 10^40)",
    )


if __name__ == "__main__":
    tx, ty = witness()
    part_x1(tx)
    part_x2(tx, ty)
    part_x3()
    part_x4((tx, ty))
    print(f"SUMMARY s6_author_routine_crosscheck: {len(FAILS)} FAIL")
    sys.exit(1 if FAILS else 0)
