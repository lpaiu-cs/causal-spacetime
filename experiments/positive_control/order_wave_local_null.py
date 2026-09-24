"""Shared-reference null-coordinate differences; no coordinate/phase test input."""

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
from order_wave_null_coordinate import (
    ROOT,
    Box,
    WaveFixture,
    alpha_box,
    cross_check,
    inverse_z,
    monotonicity_certificate,
    null_v_box,
    p0_box,
    population_v,
    series,
    t_box,
)

DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/local_null_challenge.md"
GRID = [(q, r) for q in (0.25, 1.0) for r in (64, 128, 256, 512, 1024, 2048)]


def rounded(value, bits):
    scale = 1 << bits
    return Box(
        F((value.lo * scale).__floor__(), scale),
        F((value.hi * scale).__ceil__(), scale),
    )


def square(value):
    return Box(
        F(0) if value.lo <= 0 <= value.hi else min(value.lo**2, value.hi**2),
        max(value.lo**2, value.hi**2),
    )


def magnitude(value):
    return max(abs(value.lo), abs(value.hi))


def cr(z, r, bits):
    """C_r is positive and decreasing on [0,1]."""
    return rounded(
        Box(series(z.hi, r, bits + 16).lo, series(z.lo, r, bits + 16).hi), bits
    )


def scaled_tangent(u, q, bits):
    def at(value):
        return (
            value
            * series(q * value**2, 1, bits + 16)
            / series(q * value**2, 0, bits + 16)
        )

    return rounded(Box(at(u.lo).lo, at(u.hi).hi), bits)


def shape_derivatives(s, q, bits):
    """Interval values of t,t',t'' and A,A',A'' on the whole s interval."""
    z = q * square(s)
    c1, c2, c3, c4 = (cr(z, r, bits) for r in (1, 2, 3, 4))
    u0 = rounded(c2 / c1, bits)
    u1 = rounded(-2 * c3 / c1 + square(c2 / c1), bits)
    u2 = rounded(
        6 * c4 / c1 - 6 * c2 * c3 / square(c1) + 2 * c2 * c2 * c2 / (c1 * c1 * c1), bits
    )
    alpha = rounded(square(s) * u0, bits)
    alpha1 = rounded(2 * s * u0 + 2 * q * s * square(s) * u1, bits)
    alpha2 = rounded(
        2 * u0 + 10 * q * square(s) * u1 + 4 * q * q * square(square(s)) * u2, bits
    )
    gamma = alpha_box(F(1), q, bits + 16)
    tangent = scaled_tangent(s - F(1, 2), q, bits)
    half = scaled_tangent(Box.point(F(1, 2)), q, bits)
    total = 2 * half
    t0 = rounded((tangent + half) / total, bits)
    t1 = rounded((1 + q * square(tangent)) / total, bits)
    t2 = rounded(2 * q * tangent * (1 + q * square(tangent)) / total, bits)
    b0 = rounded((gamma / alpha).sqrt(bits), bits)
    b1 = rounded(-b0 * alpha1 / (2 * alpha), bits)
    b2 = rounded(b0 * (F(3, 4) * square(alpha1 / alpha) - alpha2 / (2 * alpha)), bits)
    return (t0, t1, t2), (
        rounded(t0 * b0, bits),
        rounded(t1 * b0 + t0 * b1, bits),
        rounded(t2 * b0 + 2 * t1 * b1 + t0 * b2, bits),
    )


def derivative_boxes(theta, q, bits=96):
    """Enclose gradient and Hessian of F(p,f,a) over a Cartesian moment box."""
    q = F(q)
    if not 0 < q <= 1:
        raise ValueError("require exact known 0<q<=1")
    p, f, a = theta
    if p.lo <= 0 or f.lo <= 0:
        return None
    c = rounded(a / square(p), bits)
    if c.lo <= p0_box(q, bits + 16).hi or c.hi >= F(1, 10):
        return None
    # Split the inverse coordinate to limit interval dependency inflation.
    pieces = [
        derivatives_on_chain_band(
            p,
            f,
            Box(c.lo + (c.hi - c.lo) * F(i, 8), c.lo + (c.hi - c.lo) * F(i + 1, 8)),
            q,
            bits,
        )
        for i in range(8)
    ]
    if any(piece is None for piece in pieces):
        return None
    gradient = [
        Box(
            min(piece[0][i].lo for piece in pieces),
            max(piece[0][i].hi for piece in pieces),
        )
        for i in range(3)
    ]
    hessian = [
        [
            Box(
                min(piece[1][i][j].lo for piece in pieces),
                max(piece[1][i][j].hi for piece in pieces),
            )
            for j in range(3)
        ]
        for i in range(3)
    ]
    return gradient, hessian


def derivatives_on_chain_band(p, f, c, q, bits):
    """Enclose derivatives where a/p^2 lies in this band; a remains a coordinate."""
    z = Box(inverse_z(c.hi, q, bits).lo, inverse_z(c.lo, q, bits).hi)
    s = rounded((z / q).sqrt(bits), bits)
    if s.lo <= 0 or s.hi >= 1:
        return None
    c1, c2, c3, c4, c5 = (cr(z, r, bits) for r in range(1, 6))
    prob = rounded(c1 * c3 / (3 * square(c2)), bits)
    log1 = rounded(-c2 / c1 - 3 * c4 / c3 + 4 * c3 / c2, bits)
    log2 = rounded(
        2 * c3 / c1
        - square(c2 / c1)
        + 12 * c5 / c3
        - 9 * square(c4 / c3)
        - 12 * c4 / c2
        + 8 * square(c3 / c2),
        bits,
    )
    prob1 = rounded(prob * log1, bits)
    # The existing whole-[0,1] certificate proves p0'<-0.09/50.
    prob1 = Box(prob1.lo, min(prob1.hi, -F(9, 5000)))
    prob2 = rounded(prob * (square(log1) + log2), bits)
    z1 = rounded(Box.point(1) / prob1, bits)
    z2 = rounded(-prob2 / (prob1 * prob1 * prob1), bits)
    s1 = rounded(z1 / (2 * q * s), bits)
    s2 = rounded(z2 / (2 * q * s) - square(z1) / (4 * q * q * s * square(s)), bits)
    zero = Box.point(0)
    cg = [
        rounded(-2 * c / p, bits),
        zero,
        rounded(Box.point(1) / square(p), bits),
    ]
    ch = [
        [
            rounded(6 * c / square(p), bits),
            zero,
            rounded(Box.point(-2) / (p * square(p)), bits),
        ],
        [zero, zero, zero],
        [rounded(Box.point(-2) / (p * square(p)), bits), zero, zero],
    ]
    sg = [rounded(s1 * item, bits) for item in cg]
    sh = [
        [rounded(s2 * cg[i] * cg[j] + s1 * ch[i][j], bits) for j in range(3)]
        for i in range(3)
    ]
    t, left = shape_derivatives(s, q, bits)
    _, right = shape_derivatives(1 - s, q, bits)
    rp, rf = p.sqrt(bits), f.sqrt(bits)
    gs = rounded(-t[1] + left[1] * rp + right[1] * rf, bits)
    gss = rounded(-t[2] + left[2] * rp - right[2] * rf, bits)
    gp, gf = rounded(left[0] / (2 * rp), bits), rounded(-right[0] / (2 * rf), bits)
    gsp, gsf = rounded(left[1] / (2 * rp), bits), rounded(right[1] / (2 * rf), bits)
    gpp, gff = (
        rounded(-left[0] / (4 * p * rp), bits),
        rounded(right[0] / (4 * f * rf), bits),
    )
    gradient = [
        rounded(gs * sg[i] + (gp if i == 0 else gf if i == 1 else zero), bits)
        for i in range(3)
    ]
    hessian = []
    for i in range(3):
        row = []
        for j in range(3):
            value = gss * sg[i] * sg[j] + gs * sh[i][j]
            value += gsp * ((sg[j] if i == 0 else zero) + (sg[i] if j == 0 else zero))
            value += gsf * ((sg[j] if i == 1 else zero) + (sg[i] if j == 1 else zero))
            if i == j == 0:
                value += gpp
            if i == j == 1:
                value += gff
            row.append(rounded(value, bits))
        hessian.append(row)
    return gradient, hessian


def calibrated(theta, q, bits):
    p, f, a = theta
    return null_v_box(Box.point(p), Box.point(f), Box.point(a / p**2), q, bits)


def mean_box(mean, effective_n, exponent, bits):
    if not 0 <= mean <= 1 or effective_n < 1:
        raise ValueError("invalid bounded mean or effective sample size")
    loose = Box.point(F(exponent, 2 * effective_n)).sqrt(bits).hi
    upper = min(F(1), mean + loose)
    sharp = Box.point(2 * upper * exponent / effective_n).sqrt(bits).hi + F(
        2 * exponent, 3 * effective_n
    )
    radius = min(loose, sharp)
    return Box(max(F(0), mean - radius), min(F(1), mean + radius))


def difference_interval(first, second, references, q, exponent, bits=96):
    """Certified F(second)-F(first) for an observed first<second probe pair."""
    sizes, signs = (references, references, references // 2), (1, -1, 1)
    base = [
        mean_box(value, size, exponent, bits)
        for value, size in zip(first, sizes, strict=True)
    ]
    observed = [second[i] - first[i] for i in range(3)]
    increments = [
        sign * mean_box(sign * value, size, exponent, bits)
        for value, size, sign in zip(observed, sizes, signs, strict=True)
    ]
    domain = [
        box + Box(min(F(0), change.lo), max(F(0), change.hi))
        for box, change in zip(base, increments, strict=True)
    ]
    derivatives = derivative_boxes(domain, q, bits)
    if derivatives is None:
        return None
    gradient, hessian = derivatives
    base_error = [
        max(value - box.lo, box.hi - value)
        for value, box in zip(first, base, strict=True)
    ]
    change_error = [
        max(value - box.lo, box.hi - value)
        for value, box in zip(observed, increments, strict=True)
    ]
    radius = sum(magnitude(gradient[i]) * change_error[i] for i in range(3))
    radius += sum(
        magnitude(hessian[i][j]) * base_error[i] * magnitude(increments[j])
        for i in range(3)
        for j in range(3)
    )
    left, right = calibrated(first, q, bits), calibrated(second, q, bits)
    if left is None or right is None:
        return None
    return right - left + Box(-radius, radius)


def order_test(order, q, rng):
    """Order-only, wavelength-free detection with analytical level at most N^-2."""
    order = np.asarray(order)
    if order.dtype != bool or order.ndim != 2 or order.shape[0] != order.shape[1]:
        raise ValueError("square boolean transitive-order matrix required")
    n = len(order)
    if n < 4 or np.any(np.diag(order)) or np.any(order & order.T):
        raise ValueError("at least four events in a strict order required")
    q = F(q)
    if not 0 < q <= 1:
        raise ValueError("require known 0<q<=1")
    rows = [int.from_bytes(np.packbits(row).tobytes(), "big") for row in order]
    for i in range(n):
        if any(rows[j] & ~rows[i] for j in np.flatnonzero(order[i])):
            raise ValueError("input is not transitively closed")
    shuffled = rng.permutation(n)
    reference, probes = shuffled[: n // 2], shuffled[n // 2 :]
    m = len(reference)
    exponent, bits = (32 * n**4).bit_length(), 40 + n.bit_length()
    moments = {}
    # ponytail: cubic counts and quadratic rational bounds; batch before a power run.
    for i in probes:
        ancestors = reference[order[reference, i]]
        count = int(np.count_nonzero(order[np.ix_(ancestors, ancestors)]))
        moments[int(i)] = (
            F(len(ancestors), m),
            F(int(np.count_nonzero(order[i, reference])), m),
            F(count, math.comb(m, 2)),
        )
    for i in probes:
        for j in probes[order[i, probes]]:
            result = difference_interval(
                moments[int(i)], moments[int(j)], m, q, exponent, bits
            )
            if result is not None and result.hi < 0:
                return True
    return False


def background_moments(s, q=F(1), bits=128):
    """Control only: canonical radius 1/5 and null time 1/2."""
    total = 2 * scaled_tangent(Box.point(F(1, 2)), q, bits)
    m = t_box(s, q, bits) * total
    eta = F(1, 2) - Box.point(F(1, 25)) / (2 * m)
    zeta = F(1, 2) - Box.point(F(1, 25)) / (2 * (total - m))
    gamma = alpha_box(F(1), q, bits)
    p = alpha_box(s, q, bits) * square(eta) / gamma
    f = alpha_box(1 - s, q, bits) * square(zeta) / gamma
    a = square(p) * p0_box(q * s * s, bits)
    return tuple((item.lo + item.hi) / 2 for item in (p, f, a))


def self_check():
    monotonicity_certificate()
    centre = background_moments(F(7, 10))
    steps = (F(1, 10**7), F(1, 10**7), F(1, 10**10))
    domain = [
        Box(value - step, value + step)
        for value, step in zip(centre, steps, strict=True)
    ]
    gradient, hessian = derivative_boxes(domain, F(1), 128)
    origin = calibrated(centre, F(1), 128)
    for i, step in enumerate(steps):
        plus, minus = list(centre), list(centre)
        plus[i] += step
        minus[i] -= step
        vp, vm = calibrated(plus, F(1), 128), calibrated(minus, F(1), 128)
        finite_first = (vp - vm) / (2 * step)
        finite_second = (vp - 2 * origin + vm) / (step * step)
        assert gradient[i].lo <= finite_first.lo <= finite_first.hi <= gradient[i].hi
        assert (
            hessian[i][i].lo <= finite_second.lo <= finite_second.hi <= hessian[i][i].hi
        )
    for i in range(3):
        for j in range(i):
            mixed = Box.point(0)
            for si in (-1, 1):
                for sj in (-1, 1):
                    corner = list(centre)
                    corner[i] += si * steps[i]
                    corner[j] += sj * steps[j]
                    mixed += si * sj * calibrated(corner, F(1), 128)
            mixed /= 4 * steps[i] * steps[j]
            assert hessian[i][j].lo <= mixed.lo <= mixed.hi <= hessian[i][j].hi
    other = background_moments(F(70001, 100000))
    interval = difference_interval(centre, other, 10**18, F(1), 256, 128)
    assert interval is not None and interval.lo <= 0 <= interval.hi
    assert not order_test(np.zeros((8, 8), dtype=bool), 1, np.random.default_rng(7))
    invalid = np.zeros((4, 4), dtype=bool)
    invalid[0, 1] = invalid[1, 2] = True
    try:
        order_test(invalid, 1, np.random.default_rng(7))
    except ValueError:
        pass
    else:
        raise AssertionError("nontransitive input accepted")
    print(
        json.dumps(
            {
                "controls": "derivatives, background difference, strict order",
                "background_difference_interval": [
                    float(interval.lo),
                    float(interval.hi),
                ],
                "gradient_absolute_bounds": [
                    float(magnitude(value)) for value in gradient
                ],
                "hessian_absolute_bounds": [
                    [float(magnitude(value)) for value in row] for row in hessian
                ],
            },
            indent=2,
        )
    )


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit source and design; run from a clean worktree")
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    stamp = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )
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
        "scope": "short-separation population witnesses; not finite-N power",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    started = time.time()
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:
        for q, r in GRID:
            geometry = WaveFixture(q, r)
            u1, u2 = 0.2, 0.2 + math.pi / r
            _, vectors = np.linalg.eigh(geometry.F_matrix(u1) - geometry.F_matrix(u2))
            z = 0.2 * vectors[:, -1]
            first = geometry.canonical_point(u1, z, 0.5)
            second = geometry.canonical_point(u2, z, 0.5)
            initial_margin = geometry.margin(first, second)
            second[1] -= max(initial_margin, 0) / 2
            moments = [geometry.moments(point) for point in (first, second)]
            recovered = [population_v(*row, q)[0] for row in moments]
            theta = [np.array([p, f, p * p * c]) for p, f, c in moments]
            row = {
                "q": q,
                "r": r,
                "phase": 0.37,
                "points": [first.tolist(), second.tolist()],
                "moments": moments,
                "initial_wave_margin": initial_margin,
                "wave_margin": geometry.margin(first, second),
                "population_v_decrease": recovered[0] - recovered[1],
                "scaled_population_decrease": (recovered[0] - recovered[1]) * r**2,
                "scaled_nonnegative_increments": (
                    (theta[1] - theta[0]) * [1, -1, 1] * r
                ).tolist(),
                "ode_evaluations": geometry.sol.nfev,
            }
            if (q, r) == (1, 2048):
                synthetic_m = 10**30
                denominators = (synthetic_m, synthetic_m, math.comb(synthetic_m, 2))
                synthetic = [
                    tuple(
                        F(round(F(float(value)) * count), count)
                        for value, count in zip(item, denominators, strict=True)
                    )
                    for item in theta
                ]
                bound = difference_interval(
                    *synthetic,
                    synthetic_m,
                    F(1),
                    (32 * (2 * synthetic_m) ** 4).bit_length(),
                    128,
                )
                row["synthetic_count_control"] = {
                    "scope": "population moments; no sampled order or measured power",
                    "reference_count": synthetic_m,
                    "difference_interval": None
                    if bound is None
                    else [float(bound.lo), float(bound.hi)],
                }
                refined = WaveFixture(q, r, tighter=True)
                refined_moments = [
                    refined.moments(point, 2) for point in (first, second)
                ]
                values = [population_v(*item, q)[0] for item in refined_moments]
                row["refinement"] = {
                    "moments": refined_moments,
                    "population_v_decrease": values[0] - values[1],
                    "wave_margin": refined.margin(first, second),
                }
            stream.write(json.dumps(row) + "\n")
            stream.flush()
            print(
                json.dumps(
                    {
                        key: value
                        for key, value in row.items()
                        if key not in ("points", "moments")
                    }
                ),
                flush=True,
            )
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
        assert (
            hashlib.sha256(
                subprocess.check_output(["git", "show", f"{stamp}:{path}"], cwd=ROOT)
            ).hexdigest()
            == digest
        )
    rows = [json.loads(line) for line in raw.splitlines()]
    assert [(row["q"], row["r"]) for row in rows] == GRID
    result = {
        "code_version": stamp,
        "raw_sha256": manifest["rows_sha256"],
        "elapsed_seconds": manifest["elapsed_seconds"],
        "positive_population_witnesses": sum(
            row["population_v_decrease"] > 0 and row["wave_margin"] > 0 for row in rows
        ),
        "rows": [
            {
                key: row[key]
                for key in (
                    "q",
                    "r",
                    "scaled_population_decrease",
                    "scaled_nonnegative_increments",
                )
            }
            for row in rows
        ],
        "refinement": [
            {
                "q": row["q"],
                "r": row["r"],
                "population_decrease_change": row["refinement"]["population_v_decrease"]
                - row["population_v_decrease"],
            }
            for row in rows
            if "refinement" in row
        ],
        "synthetic_count_controls": [
            row["synthetic_count_control"]
            for row in rows
            if "synthetic_count_control" in row
        ],
    }
    (folder / "summary.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


def recheck_enclosure(folder):
    """Post-run check on unchanged synthetic moments; preserve the primary failure."""
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit enclosure revision and plan before rechecking")
    path = folder / "post_run_enclosure_check.json"
    if path.exists():
        raise FileExistsError(path)
    raw = (folder / "rows.jsonl").read_bytes()
    manifest = json.loads((folder / "manifest.json").read_text())
    assert hashlib.sha256(raw).hexdigest() == manifest["rows_sha256"]
    row = next(
        item
        for item in map(json.loads, raw.splitlines())
        if (item["q"], item["r"]) == (1, 2048)
    )
    count = row["synthetic_count_control"]["reference_count"]
    assert count == 10**30
    denominators = (count, count, math.comb(count, 2))
    moments = [
        tuple(
            F(round(F(float(value)) * den), den)
            for value, den in zip((p, f, p * p * c), denominators, strict=True)
        )
        for p, f, c in row["moments"]
    ]
    bound = difference_interval(
        *moments, count, F(1), (32 * (2 * count) ** 4).bit_length(), 128
    )
    plan = DESIGN.with_name("local_null_results.md")
    result = {
        "scope": "post-run eight-band enclosure; synthetic moments, not sampled power",
        "code_version": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
        .decode()
        .strip(),
        "files": {
            item.relative_to(ROOT).as_posix(): hashlib.sha256(
                item.read_bytes()
            ).hexdigest()
            for item in (Path(__file__), plan)
        },
        "raw_sha256": manifest["rows_sha256"],
        "reference_count": count,
        "bands": 8,
        "original_difference_interval": row["synthetic_count_control"][
            "difference_interval"
        ],
        "difference_interval": None
        if bound is None
        else [float(bound.lo), float(bound.hi)],
        "negative_upper_endpoint": bound is not None and bound.hi < 0,
    }
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summarize", type=Path)
    parser.add_argument("--cross-check", type=Path)
    parser.add_argument("--recheck-enclosure", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
    if args.summarize:
        summarize(args.summarize)
    if args.recheck_enclosure:
        recheck_enclosure(args.recheck_enclosure)
    if args.cross_check:
        if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
            raise RuntimeError("Run the independent formulation from a clean worktree")
        cross_check(args.cross_check)
        path = args.cross_check / "independent_formulation_check.json"
        result = json.loads(path.read_text())
        result["code_version"] = (
            subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
            .decode()
            .strip()
        )
        result["caller_source_sha256"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        path.write_text(json.dumps(result, indent=2))
