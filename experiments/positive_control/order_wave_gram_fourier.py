"""Population challenge for an intrinsic anchor Gram residual, not a power run."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from order_wave_energy import background_probability, rule
from order_wave_null_coordinate import Box, WaveFixture, population_v, series
from scipy.special import hyp0f1

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/gram_fourier_challenge.md"


def alpha(length, q):
    z = -q * np.asarray(length) ** 2 / 4
    return np.asarray(length) ** 2 * hyp0f1(2.5, z) / (6 * hyp0f1(1.5, z))


def m_coordinate(u, q):
    k = math.sqrt(q)
    return (np.tan(k * u) + math.tan(k / 2)) / k


def state(p, f, chain, q):
    """Recover (u,m,v,radius squared) from order probabilities, chain scalar."""
    v, u = population_v(p, f, chain, q)
    m, total = m_coordinate(u, q), 2 * m_coordinate(0, q)
    eta = np.sqrt(p * alpha(1, q) / alpha(u + 0.5, q))
    zeta = np.sqrt(f * alpha(1, q) / alpha(0.5 - u, q))
    radius = 2 * m * (total - m) / total * (1 - eta - zeta)
    return np.stack(np.broadcast_arrays(u, m, v, radius), axis=-1)


def gram(first, second, interval_probability, q):
    """Background inner product expressed in calibrated states and order counts."""
    gap = second[..., 0] - first[..., 0]
    if np.any(gap <= 0) or np.any(interval_probability < 0):
        raise ValueError(
            "positive calibrated time gap and interval probability required"
        )
    slack = np.sqrt(interval_probability * alpha(1, q) / alpha(gap, q))
    return (first[..., 3] + second[..., 3]) / 2 - (second[..., 1] - first[..., 1]) * (
        second[..., 2] - first[..., 2] - slack
    )


def bump(value):
    result = np.zeros_like(value, dtype=float)
    inside = np.abs(value) < 1
    result[inside] = np.exp(-1 / (1 - value[inside] ** 2))
    return result


def weight(states):
    return (
        bump((states[..., 0] - 0.2) / 0.1)
        * bump((states[..., 2] - 0.7) / 0.05)
        * bump((states[..., 3] - 0.02125) / 0.01875)
    )


def anchors(theta):
    rotation = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )
    return (
        np.array([-0.35, -0.25, -0.15]),
        np.array([0.15, 0.25, 0.35]),
        0.08 * np.array([[1, 0], [0, 1], [1, 1]]) @ rotation.T,
    )


def predicted(q, theta):
    x, w = rule(512)
    u = 0.1 + 0.2 * x
    w = w * bump((u - 0.2) / 0.1)
    ma = m_coordinate(anchors(theta)[0], q)
    gaps = m_coordinate(u, q)[:, None] - ma
    coefficient = 0.08**2 * ((1 / gaps[:, 0] - 1 / gaps[:, 1]) / 2 + 1j / gaps[:, 2])
    return (
        -1j
        * math.sqrt(q)
        * np.exp(1j * (0.37 - 2 * theta))
        * np.sum(w * coefficient)
        / np.sum(w * np.cos(math.sqrt(q) * u) ** 2)
    )


def population(wave, theta, factor=1):
    q, r = wave.q, wave.r
    ua, va, za = anchors(theta)
    points = np.array(
        [wave.canonical_point(u, z, v) for u, v, z in zip(ua, va, za, strict=True)]
    )
    recovered = np.array([state(*wave.moments(p, factor), q) for p in points])
    gamma = wave.segment(-0.5, 0.5, factor)[0]
    anchor_gram = np.diag(recovered[:, 3])
    minimum_slack = min(
        wave.margin(points[i], points[j]) for i, j in ((0, 1), (0, 2), (1, 2))
    )
    for i, j in ((0, 1), (0, 2), (1, 2)):
        slack = wave.margin(points[i], points[j])
        if slack <= 0:
            raise ValueError("anchor pair is not timelike")
        b = wave.segment(ua[i], ua[j], factor)[0] * slack**2 / gamma
        anchor_gram[i, j] = anchor_gram[j, i] = gram(recovered[i], recovered[j], b, q)
    block = anchor_gram[:2, :2]
    if np.linalg.eigvalsh(block)[0] <= 0:
        raise ValueError("anchor Gram block is not positive definite")
    ell = np.r_[-np.linalg.solve(block, anchor_gram[:2, 2]), 1.0]
    x, wx = rule(factor * max(96, math.ceil(r / 4) + 48))
    vnodes, vw = rule(factor * 16)
    tnodes, tw = rule(factor * 16)
    v, t = np.meshgrid(0.60 + 0.20 * vnodes, 0.0001 + 0.0599 * tnodes, indexing="ij")
    area_weight = vw[:, None] * tw[None, :] * 0.20 * 0.0599
    angles = np.arange(4) * np.pi / 2
    z = np.sqrt(t)[..., None, None] * np.stack(
        (np.cos(angles), np.sin(angles)), axis=-1
    )
    total, mass = 0j, 0.0
    max_state_error = 0.0
    for u, wu in zip(0.05 + 0.30 * x, 0.30 * wx, strict=True):
        m = m_coordinate(u, q)
        eta, zeta = v - t / (2 * m), 1 - v - t / (2 * (wave.M - m))
        if min(np.min(eta), np.min(zeta)) <= 0:
            raise ValueError("integration point left the diamond")
        prefix, chain = wave.segment(-0.5, u, factor, pairs=True)
        suffix = wave.segment(u, 0.5, factor)[0]
        states = state(prefix * eta**2 / gamma, suffix * zeta**2 / gamma, chain, q)
        weights = weight(states)
        if not np.any(weights):
            continue
        max_state_error = max(
            max_state_error, float(np.max(np.abs(states[..., 0] - u)))
        )
        left, right = wave.matrix(-0.5, u), wave.matrix(u, 0.5)
        values, vectors = np.linalg.eigh(
            m * (wave.M - m) / wave.M * (np.linalg.inv(left) + np.linalg.inv(right))
        )
        transform = (vectors * values**-0.5) @ vectors.T
        y = z @ transform.T
        actual_v = eta[..., None] + 0.5 * np.einsum(
            "...i,ij,...j->...", y, np.linalg.inv(left), y
        )
        h = np.zeros_like(actual_v)
        for i in range(3):
            dy = y - points[i, 2:]
            slack = (
                actual_v
                - points[i, 1]
                - 0.5
                * np.einsum(
                    "...i,ij,...j->...", dy, np.linalg.inv(wave.matrix(ua[i], u)), dy
                )
            )
            active_slack = slack[weights > 0]
            minimum_slack = min(minimum_slack, float(np.min(active_slack)))
            if np.min(active_slack) <= 0:
                raise ValueError("nonzero-weight probe is not in all anchor futures")
            b = wave.segment(ua[i], u, factor)[0] * np.maximum(slack, 0) ** 2 / gamma
            h += ell[i] * gram(recovered[i], states[..., None, :], b, q)
        density = math.cos(math.sqrt(q) * u) ** 2 * np.linalg.det(transform) / gamma
        measure = wu * density * area_weight * weights
        mass += float(np.sum(measure))
        total += np.sum(
            measure * np.mean(h, axis=-1) * np.exp(-1j * r * states[..., 0])
        )
    scaled = total / mass * r**2
    prediction = predicted(q, theta)
    return {
        "q": q,
        "r": r,
        "theta": theta,
        "phase": 0.37,
        "scaled_complex_signal": [float(scaled.real), float(scaled.imag)],
        "predicted_complex_limit": [float(prediction.real), float(prediction.imag)],
        "relative_prediction_error": float(abs(scaled - prediction) / abs(prediction)),
        "weight_mass": mass,
        "anchor_gram": anchor_gram.tolist(),
        "ell": ell.tolist(),
        "anchor_block_condition": float(np.linalg.cond(block)),
        "minimum_active_causal_slack": minimum_slack,
        "maximum_active_u_calibration_error": max_state_error,
        "quadrature_factor": factor,
        "primitive_evaluations": wave.sol.nfev,
    }


def self_check():
    q, total = 1.0, 2 * m_coordinate(0, 1)
    ua, va, za = anchors(0.2)
    u = np.r_[ua, 0.2]
    v = np.r_[va, 0.7]
    z = np.vstack((za, [0.1, 0.12]))
    m = m_coordinate(u, q)
    radius = np.sum(z**2, axis=-1)
    eta, zeta = v - radius / (2 * m), 1 - v - radius / (2 * (total - m))
    recovered = np.array(
        [
            state(
                alpha(s + 0.5, q) * e**2 / alpha(1, q),
                alpha(0.5 - s, q) * f**2 / alpha(1, q),
                background_probability(q * (s + 0.5) ** 2),
                q,
            )
            for s, e, f in zip(u, eta, zeta, strict=True)
        ]
    )
    assert np.max(np.abs(recovered - np.stack((u, m, v, radius), axis=-1))) < 2e-12
    g = np.diag(recovered[:, 3])
    for i in range(4):
        for j in range(i + 1, 4):
            slack = v[j] - v[i] - np.sum((z[j] - z[i]) ** 2) / (2 * (m[j] - m[i]))
            assert slack > 0
            b = alpha(u[j] - u[i], q) * slack**2 / alpha(1, q)
            g[i, j] = g[j, i] = gram(recovered[i], recovered[j], b, q)
    assert np.max(np.abs(g - z @ z.T)) < 2e-12
    ell = np.r_[-np.linalg.solve(g[:2, :2], g[:2, 2]), 1.0]
    assert abs(ell @ g[:3, 3]) < 2e-12
    means = []
    for count in (4, 4096):
        theta = 2 * np.pi * np.arange(count) / count
        y = 0.2 * np.stack((np.cos(theta), np.sin(theta)), axis=-1)
        means.append(
            np.mean(
                np.einsum("...i,ij,...j->...", y, [[1.3, 0.4], [0.4, -0.7]], y)
                + y @ [0.3, -0.2]
                + 0.1
            )
        )
    assert abs(means[0] - means[1]) < 1e-14
    print(
        json.dumps(
            {
                "controls": "passed",
                "background_gram_max_error": float(np.max(np.abs(g - z @ z.T))),
            }
        )
    )


def guard_check(output):
    """Exact-rational outer-anchor bounds for the subsequent proof draft."""
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("commit the guard proof and source before checking")
    if output.exists():
        raise FileExistsError(output)
    small = Box(F(638, 100000), F(642, 100000))
    large = Box(F(1278, 100000), F(1282, 100000))
    off = Box(-F(2, 100000), F(2, 100000))
    x1 = small.sqrt(128)
    x2 = off / x1
    y2 = (small - x2 * x2).sqrt(128)
    x3 = small / x1
    y3 = (small - x2 * x3) / y2
    ell2 = -(small - off * small / small) / (small - off * off / small)
    ell_bound = max(abs(ell2.lo), abs(ell2.hi))
    cosine = series(F(1, 4), 0, 128)
    derivative_bound = (Box.point(1) / (cosine * cosine)).hi
    assert derivative_bound < F(13, 10)
    assert min(x3.lo, y3.lo) > F(79, 1000)
    assert ell_bound < F(102, 100)
    determinant = small * small - off * off
    assert determinant.lo > 0
    imaginary_lower = x3.lo * y3.lo / (
        F(13, 10) * F(452, 1000)
    ) - ell_bound * off.hi / F(348, 1000)
    absolute_upper = (
        ell_bound * small.hi / F(448, 1000)
        + ell_bound * small.hi / F(348, 1000)
        + large.hi / F(248, 1000)
    ) / 2
    assert imaginary_lower > F(1, 100) and absolute_upper < F(1, 20)
    # A frequency-grid error <=1/8 rotates phases by <=0.1/8 about u=0.2.
    grid_lower = imaginary_lower - absolute_upper / 80
    assert grid_lower > F(9, 1000)
    radii, cross = [small, small, large], {(0, 1): off, (0, 2): small, (1, 2): small}
    ua = [F(-35, 100), F(-25, 100), F(-15, 100)]
    va = [F(15, 100), F(25, 100), F(35, 100)]
    pair_slacks, probe_slacks = [], []
    for (i, j), g in cross.items():
        gap = ua[j] - ua[i] - F(4, 1000)
        distance = (radii[i] + radii[j] - 2 * g).hi
        pair_slacks.append(va[j] - va[i] - F(4, 1000) - distance / (2 * gap))
    for i in range(3):
        distance = (F(1, 5) + radii[i].sqrt(128).hi) ** 2
        gap = F(1, 10) - ua[i] - F(2, 1000)
        probe_slacks.append(F(65, 100) - va[i] - F(2, 1000) - distance / (2 * gap))
    assert min(pair_slacks + probe_slacks) > 0
    files = [Path(__file__), DESIGN.with_name("gram_fourier_noise.md")]
    result = {
        "code_version": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
        .decode()
        .strip(),
        "scope": "rational anchor bounds; not finite-order noise or power",
        "files": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
        "passed": True,
        "m_derivative_upper": float(derivative_bound),
        "third_anchor_coordinate_lower": float(min(x3.lo, y3.lo)),
        "ell_magnitude_upper": float(ell_bound),
        "D_imaginary_lower": float(imaginary_lower),
        "D_absolute_upper": float(absolute_upper),
        "grid_coefficient_lower": float(grid_lower),
        "anchor_pair_slack_lower": float(min(pair_slacks)),
        "active_probe_slack_lower": float(min(probe_slacks)),
    }
    output.write_text(json.dumps(result, indent=2))
    print(json.dumps(result))


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("commit source and design before running")
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), DESIGN] + [
        Path(__file__).with_name(name)
        for name in (
            "order_wave_null_coordinate.py",
            "order_wave_resolution.py",
            "order_wave_energy.py",
        )
    ]
    manifest = {
        "status": "in_progress",
        "code_version": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
        .decode()
        .strip(),
        "utc_started": time.strftime("%FT%TZ", time.gmtime()),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "scope": "deterministic population Fourier coefficient, not finite-order power",
        "files": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    started = time.time()
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:
        for r in (64, 128, 256, 512):
            wave = WaveFixture(1, r)
            for theta in (0.0, np.pi / 8):
                try:
                    row = population(wave, theta)
                    if r == 512 and theta == 0:
                        refined = population(WaveFixture(1, r, tighter=True), theta, 2)
                        main_signal = complex(*row["scaled_complex_signal"])
                        refined_signal = complex(*refined["scaled_complex_signal"])
                        refined["relative_signal_change"] = abs(
                            refined_signal - main_signal
                        ) / abs(main_signal)
                        refined["numerical_consistency_pass"] = (
                            refined["relative_signal_change"] < 0.01
                        )
                        row["refinement"] = refined
                except (ValueError, RuntimeError) as error:
                    row = {"q": 1, "r": r, "theta": theta, "failure": str(error)}
                stream.write(json.dumps(row) + "\n")
                stream.flush()
                print(json.dumps(row), flush=True)
    manifest.update(
        status="completed",
        elapsed_seconds=time.time() - started,
        rows_sha256=hashlib.sha256((output / "rows.jsonl").read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--guard-output", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
    if args.guard_output:
        guard_check(args.guard_output)
