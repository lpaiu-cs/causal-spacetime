"""Coordinate-oracle barrier: deterministic geometry checks, not an intrinsic test.

--output DIR runs only the frozen grid in coordinate_oracle_barrier.md.
The finite-lambda integral is exact loss on reference-coordinate pairs;
twice that loss is NOT the exact full disagreement probability.
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
from order_wave_energy import rule
from order_wave_resolution import ROOT, circular_primitive, convolution
from scipy.special import ellipe

DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/coordinate_oracle_barrier.md"
GRID = [(q, r) for q in (0.25, 1.0, 4.0) for r in (32, 64, 128, 256, 512, 1024)]


def parts(q, s, t):
    if not (0 < q < math.pi**2):
        raise ValueError("require fixed 0 < q < pi^2")
    k = math.sqrt(q)
    a, b = (np.tan(k * s) + np.tan(k / 2)) / k, (np.tan(k / 2) - np.tan(k * t)) / k
    total = 2 * math.tan(k / 2) / k
    middle = (t - s) * np.sinc(k * (t - s) / np.pi) / (np.cos(k * s) * np.cos(k * t))
    density = np.cos(k * s) ** 2 * np.cos(k * t) ** 2
    volume = math.pi * convolution(q, 2) / convolution(q, 1)
    return a, middle, b, total, density, volume


def oracle_coefficient(q, n=128):
    x, w = rule(n)
    s = -0.5 + (1 - x[:, None]) * x[None, :]
    t = s + x[:, None]
    a, _, b, total, density, volume = parts(q, s, t)
    k = math.sqrt(q)
    f, g = 1 / np.cos(k * s) ** 2, 1 / np.cos(k * t) ** 2
    parameter = np.minimum(1.0, 4 * f * g / (f + g) ** 2)
    phase_average = 2 * (f + g) / math.pi * ellipe(parameter)
    integral = np.sum(
        w[:, None]
        * w[None, :]
        * (1 - x[:, None])
        * density
        * a
        * b
        * (a + b)
        * phase_average
    )
    return float(4 * math.pi * math.sqrt(q) / (3 * volume**2 * total**2) * integral)


def loss_fraction(chi, outer_fraction, n=24):
    """Positive angular sector, avoiding subtraction of two O(1) values."""
    chi, outer_fraction = np.broadcast_arrays(chi, outer_fraction)
    if np.any((chi < 0) | (chi >= 1) | (outer_fraction < 0) | (outer_fraction > 1)):
        raise ValueError("invalid ellipticity or segment fraction")
    x, w = rule(n)
    start = 0.5 * np.arccos(chi)
    length = math.pi / 2 - start
    theta = start[..., None] + length[..., None] * x
    d = (outer_fraction * chi / (1 - chi**2))[..., None] * (
        chi[..., None] - np.cos(2 * theta)
    )
    return 2 / math.pi * length * np.sum(w * d / (1 + d), axis=-1)


def reference_loss(q, r, *, refine=False):
    """Exact background-pair loss integral; the wave sampling law differs."""
    lam = 1 / r
    sol = circular_primitive(q, lam, tighter=refine)
    multiplier = 2 if refine else 1
    x, w = rule(32 * multiplier)
    mid, mw = rule(64 * multiplier)
    boundaries = np.append(np.arange(0, 1, 2 * math.pi / r), 1.0)
    length = np.diff(boundaries)
    separations = (boundaries[:-1, None] + length[:, None] * x).ravel()
    weights = (length[:, None] * w).ravel()
    loss, linear = 0.0, 0.0
    # ponytail: bounded quadrature blocks; no parallel job framework for this CPU probe.
    for first in range(0, len(separations), 128):
        ell = separations[first : first + 128, None]
        s = -0.5 + (1 - ell) * mid[None, :]
        t = s + ell
        a, b, last, total, density, volume = parts(q, s, t)
        integral = (sol.sol(t.ravel())[0] - sol.sol(s.ravel())[0]).reshape(s.shape)
        e = 2 * lam * math.sqrt(q) / math.sqrt(1 + 4 * q * lam**2) * abs(integral)
        chi, outer = e / b, (a + last) / total
        chain = math.pi**2 / 6 * a * b * last / total
        factor = (
            2
            / volume**2
            * weights[first : first + 128, None]
            * mw[None, :]
            * (1 - ell)
            * density
            * chain
        )
        loss += np.sum(factor * loss_fraction(chi, outer, 24 * multiplier))
        linear += np.sum(factor * outer * chi / math.pi)
    return {
        "reference_loss_probability": float(loss),
        "twice_reference_loss_scaled": float(2 * loss / lam**2),
        "linear_disagreement_scaled": float(2 * linear / lam**2),
        "ode_evaluations": sol.nfev,
        "separation_nodes": len(separations),
        "midpoint_nodes": len(mid),
        "angular_nodes": 24 * multiplier,
    }


def self_check():
    expected = 32 / (5 * math.pi**2)
    assert abs(oracle_coefficient(1e-8) / 1e-4 - expected) < 2e-9
    a, b, last, eps = 0.23, 0.41, 0.36, 1e-5

    def chain(middle):
        return math.pi**2 / 6 * a * middle * last / (a + middle + last)

    derivative = math.pi**2 / 6 * a * last * (a + last) / (a + b + last) ** 2
    assert abs((chain(b + eps) - chain(b - eps)) / (2 * eps) / derivative - 1) < 2e-9
    for chi in (0.001, 0.03, 0.3):
        outer = 0.63
        start = math.acos(chi) / 2
        b1 = 1 - outer * chi / (1 + chi)
        b2 = 1 + outer * chi / (1 - chi)
        antiderivative = (
            2
            / math.pi
            * (
                math.pi / 2
                - start
                - (math.pi / 2 - math.atan(math.sqrt(b2 / b1) * math.tan(start)))
                / math.sqrt(b1 * b2)
            )
        )
        assert abs(float(loss_fraction(chi, outer)) - antiderivative) < 2e-14
    sol = circular_primitive(1, 1 / 16)
    x, w = rule(128)
    s, t = -0.31, 0.27
    u = s + (t - s) * x
    direct = (t - s) * np.sum(w * np.exp(16j * u) / np.cos(u) ** 2)
    assert abs((sol.sol(t)[0] - sol.sol(s)[0]) - direct) < 2e-12
    assert loss_fraction(0, 0.5) == 0
    print(
        "Controls passed: small-q coefficient, slice derivative, "
        "angular integral, primitive."
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
    files = [
        Path(__file__),
        Path(__file__).with_name("order_wave_resolution.py"),
        Path(__file__).with_name("order_wave_energy.py"),
        DESIGN,
    ]
    manifest = {
        "status": "in_progress",
        "code_version": stamp,
        "utc_started": time.strftime("%FT%TZ", time.gmtime()),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "files": {
            path.relative_to(ROOT).as_posix(): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in files
        },
        "scope": "exact reference-coordinate loss integral; not order-only power",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    started = time.time()
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:
        for q, r in GRID:
            row = {
                "q": q,
                "r": r,
                "kappa": oracle_coefficient(q),
                **reference_loss(q, r),
            }
            row["ratio_to_limit"] = row["twice_reference_loss_scaled"] / row["kappa"]
            if (q, r) == (1, 1024):
                row["refinement"] = {
                    "kappa": oracle_coefficient(q, 256),
                    **reference_loss(q, r, refine=True),
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


def summarize(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    raw = (folder / "rows.jsonl").read_bytes()
    assert manifest["status"] == "completed"
    assert hashlib.sha256(raw).hexdigest() == manifest["rows_sha256"]
    stamp = manifest["code_version"]
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=ROOT, check=True
    )
    for path, digest in manifest["files"].items():
        saved = subprocess.check_output(["git", "show", f"{stamp}:{path}"], cwd=ROOT)
        assert hashlib.sha256(saved).hexdigest() == digest
    rows = [json.loads(line) for line in raw.splitlines()]
    assert [(row["q"], row["r"]) for row in rows] == GRID
    refined = next(row for row in rows if "refinement" in row)
    change = abs(
        refined["refinement"]["twice_reference_loss_scaled"]
        / refined["twice_reference_loss_scaled"]
        - 1
    )
    result = {
        "code_version": stamp,
        "raw_sha256": manifest["rows_sha256"],
        "elapsed_seconds": manifest["elapsed_seconds"],
        "positive_coefficients_and_losses": all(
            row["kappa"] > 0 and row["reference_loss_probability"] > 0 for row in rows
        ),
        "ratios_to_limit": [
            {key: row[key] for key in ("q", "r", "kappa", "ratio_to_limit")}
            for row in rows
        ],
        "refinement_relative_change": change,
        "refinement_within_frozen_tolerance": change < 1e-5,
        "coefficient_refinement_relative_change": abs(
            refined["refinement"]["kappa"] / refined["kappa"] - 1
        ),
    }
    (folder / "summary.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    assert result["positive_coefficients_and_losses"]
    assert result["refinement_within_frozen_tolerance"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summarize", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
    if args.summarize:
        summarize(args.summarize)
