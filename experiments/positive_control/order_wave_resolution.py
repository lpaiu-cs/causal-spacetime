"""Fixed-q circular-wave chain asymptotics; no finite-poset accuracy claim.

--self-check compares independent kernels and the weak-energy control.
--output DIR runs only the grid frozen in full_order_chain_asymptotics.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import scipy
from order_wave_energy import (
    background_probability,
    relation_probability,
    rosen_kernel,
    rule,
    transfer,
)
from scipy.integrate import solve_ivp
from scipy.special import hyp0f1

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/full_order_chain_asymptotics.md"


def convolution(q, r):
    return float(hyp0f1(r + 0.5, -q / 4) / math.factorial(2 * r - 1))


def coefficient(q):
    return (
        4
        * q
        * convolution(q, 1)
        / convolution(q, 2)
        * (1 - 4 * background_probability(q))
    )


def circular_correction(q, lam, *, tighter=False):
    """Return K0 and (K_lambda/c-K0)/K0 without subtractive cancellation."""
    if not (0 <= q < math.pi**2 and lam > 0):
        raise ValueError("q must be in [0,pi^2); lambda must be positive")
    k = math.sqrt(q)
    sol = solve_ivp(
        lambda u, _: np.array([np.exp(1j * u / lam) / np.cos(k * u) ** 2]),
        (-0.5, 0.5),
        np.array([0j]),
        method="DOP853",
        dense_output=True,
        rtol=2e-13 if tighter else 2e-11,
        atol=2e-15 if tighter else 2e-13,
        max_step=lam / (16 if tighter else 8),
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    def kernel(s, t):
        s, t = np.broadcast_arrays(s, t)
        ell = t - s
        if np.any((ell < 0) | (s < -0.5) | (t > 0.5)):
            raise ValueError("require -1/2 <= s <= t <= 1/2")
        base = ell * np.sinc(k * ell / np.pi)
        m = base / (np.cos(k * s) * np.cos(k * t))
        integral = (sol.sol(t.ravel())[0] - sol.sol(s.ravel())[0]).reshape(s.shape)
        ratio = np.divide(abs(integral), m, out=np.ones_like(m), where=ell > 0)
        z = 4 * q * lam**2 / (1 + 4 * q * lam**2) * ratio**2
        if np.any((z < 0) | (z >= 1)):
            raise RuntimeError("loss of positive segment determinant")
        relative = -z / (1 + np.sqrt(1 - z))
        return base, relative

    return kernel, sol.nfev


def pair_change(kernel, q, n):
    """Exact hierarchy evaluated as changes, including all product cross terms."""
    x, w = rule(n)
    v = x - 0.5
    left, dl = kernel(-0.5, v)
    right, dr = kernel(v, 0.5)
    dc2 = np.sum(w * left * right * (dl + dr + dl * dr))
    u = -0.5 + x[:, None] * x[None, :]
    first, d1 = kernel(-0.5, u)
    middle, d2 = kernel(u, v[None, :])
    d3 = dr[None, :]
    change = d1 + d2 + d3 + d1 * d2 + d1 * d3 + d2 * d3 + d1 * d2 * d3
    dc3 = np.sum(
        w[:, None] * w[None, :] * x[None, :] * first * middle * right[None, :] * change
    )
    _, dk = kernel(-0.5, 0.5)
    e2, e3 = dc2 / convolution(q, 2), dc3 / convolution(q, 3)
    dp = (
        background_probability(q) * (dk + e3 + dk * e3 - 2 * e2 - e2**2) / (1 + e2) ** 2
    )
    return {"p2_change": float(dp), "c2_change": float(dc2), "c3_change": float(dc3)}


def reject_from_pair_count(comparable_pairs, n, q, alpha=0.025):
    """Order-only test; requires fixed known q, no wavelength or coordinates."""
    if n < 2 or not 0 <= comparable_pairs <= n * (n - 1) // 2:
        raise ValueError("invalid event/pair counts")
    if not (0 < q < math.pi**2 and 0 < alpha < 1):
        raise ValueError("require 0 < q < pi^2 and 0 < alpha < 1")
    observed = comparable_pairs / (n * (n - 1) / 2)
    threshold = background_probability(q) - math.sqrt(
        math.log(1 / alpha) / (2 * (n // 2))
    )
    return observed < threshold


def self_check():
    assert abs(coefficient(1e-6) / 1e-6 - 72 / 5) < 2e-6
    assert all(coefficient(q) > 0 for q in (0.25, 1, 4, 9.0))
    kernel, _ = circular_correction(1, 1 / 8)
    independent, _ = rosen_kernel(1, 1 / 8, "circular", tighter=True)
    s = np.array([-0.5, -0.31, -0.12, 0.17])
    t = np.array([0.5, -0.30, 0.25, 0.499])
    base, relative = kernel(s, t)
    c = math.sqrt(1 + 4 / 8**2)
    assert np.max(abs(c * base * (1 + relative) - independent(s, t))) < 2e-11
    got = pair_change(kernel, 1, 128)["p2_change"]
    reference = relation_probability(independent, 128) - background_probability(1)
    assert abs(got - reference) < 2e-12
    weak_q, r = 1e-5, 16
    kernel, _ = circular_correction(weak_q, 1 / r)
    got = pair_change(kernel, weak_q, 128)["p2_change"]
    predicted = -2 * weak_q / 525 * (transfer(r) - 1)
    assert abs(got / predicted - 1) < 3e-5
    assert reject_from_pair_count(0, 10_000, 1)
    assert not reject_from_pair_count(4_999_500, 10_000, 1)
    assert not reject_from_pair_count(1, 2, 1)
    print(
        "Controls passed: positive coefficient, independent Rosen kernel/P2, "
        "weak limit, test."
    )


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit design and source; run from a clean worktree")
    stamp = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    started = time.time()
    sources = [Path(__file__), Path(__file__).with_name("order_wave_energy.py"), DESIGN]
    manifest = {
        "status": "in_progress",
        "code_version": stamp,
        "files": {
            path.relative_to(ROOT).as_posix(): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in sources
        },
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "utc_started": time.strftime("%FT%TZ", time.gmtime()),
        "scope": "deterministic fixed-q asymptotics; no boundary/power campaign",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:
        for q in (0.25, 1.0, 4.0):
            for r in (32, 64, 128, 256, 512, 1024):
                n = max(128, math.ceil(r / 2) + 64)
                kernel, evaluations = circular_correction(q, 1 / r)
                row = {
                    "q": q,
                    "r": r,
                    "quadrature_order": n,
                    "ode_evaluations": evaluations,
                    **pair_change(kernel, q, n),
                    "predicted_A": coefficient(q),
                }
                row["scaled_coefficient"] = -row["p2_change"] * r**4 / math.log(r)
                row["bounded_scale_remainder"] = row["p2_change"] * r**4 + row[
                    "predicted_A"
                ] * math.log(r)
                row["c2_log_coefficient"] = -row["c2_change"] * r**4 / math.log(r)
                row["predicted_c2_coefficient"] = 8 * q * convolution(q, 1)
                row["c3_log_coefficient"] = -row["c3_change"] * r**4 / math.log(r)
                row["predicted_c3_coefficient"] = 12 * q * convolution(q, 2)
                if q == 1 and r in (128, 1024):
                    refined, evaluations = circular_correction(q, 1 / r, tighter=True)
                    row["refinement"] = {
                        "quadrature_order": 2 * n,
                        "ode_evaluations": evaluations,
                        **pair_change(refined, q, 2 * n),
                    }
                stream.write(json.dumps(row) + "\n")
                stream.flush()
                print(json.dumps(row), flush=True)
    manifest.update(
        status="completed",
        elapsed_seconds=time.time() - started,
        rows_sha256=hashlib.sha256((output / "rows.jsonl").read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
