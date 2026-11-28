"""Deterministic development diagnostics, NOT an order-only recovery campaign.

The evaluator uses coordinates/profiles to generate exact interval-volume
tables. Only those tables enter the reconstruction module. Finite sampling and
landmark selection are separate, still-unvalidated requirements.
"""

import hashlib
import itertools
import json
import math
import platform
import sys
from pathlib import Path
from unittest.mock import patch

import intrinsic_weyl
import numpy as np
from intrinsic_weyl import ETA, WEYL, fit_from_interval_counts, residual_design
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs/paper/paper_a/window_response"))
from check_plane_wave_chains import transfer  # noqa: E402

DIRECTIONS = [
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
    *itertools.product([-0.6, 0.6], repeat=3),
]
TANGENTS = np.array(
    [
        [0, 0, 0, 0],
        [1, 0, 0, 0],
        [2, 0.5, 0, 0],
        [3, 0, 0.5, 0],
        [4, 0, 0, 0.5],
        *[[t, *v] for t in (6, 9) for v in DIRECTIONS],
    ],
    dtype=float,
)


def physical_profile(name):
    def profile(u):
        if name == "flat":
            return np.zeros((2, 2))
        theta = 0 if name == "constant" else 0.31 + 1.7 * u
        amplitude = 1 if name == "constant" else 1 + 0.3 * math.sin(3 * u)
        c, s = math.cos(2 * theta), math.sin(2 * theta)
        return amplitude * np.array([[c, s], [s, -c]])

    return profile


def exact_interval_tables(profile, scale, order=48):
    """Oracle: exponentiate specified tangent vectors, then evaluate all volumes."""
    tangent = scale * TANGENTS
    times = (tangent[:, 0] + tangent[:, 3]) / math.sqrt(2)
    nodes, weights = leggauss(order)
    pairs = [
        (i, j) if times[i] < times[j] else (j, i)
        for i, j in itertools.combinations(range(len(times)), 2)
        if abs(times[i] - times[j]) > 1e-14
    ]
    intervals = [times[i] + (times[j] - times[i]) * (nodes + 1) / 2 for i, j in pairs]
    evaluation = np.unique(np.r_[times, np.concatenate(intervals)])
    fs = transfer(profile, 0, max(times), evaluation)
    inv = np.linalg.inv(fs)
    endpoints = np.searchsorted(evaluation, times)
    points = np.zeros((len(times), 4))
    points[:, 0] = times
    for i in range(1, len(times)):
        f = fs[endpoints[i]]
        transverse_velocity = tangent[i, 1:3] / times[i]
        x = f[:2, 2:] @ transverse_velocity
        xp = f[2:, 2:] @ transverse_velocity
        norm = tangent[i] @ ETA @ tangent[i]
        points[i, 1] = (times[i] * (x @ xp) - norm) / (2 * times[i])
        points[i, 2:] = x
    volumes, proper_squared = (np.zeros((len(times), len(times))) for _ in range(2))
    for (i, j), quadrature in zip(pairs, intervals, strict=True):
        f = fs[endpoints[j]] @ inv[endpoints[i]]
        a, b, d = f[:2, :2], f[:2, 2:], f[2:, 2:]
        bi = np.linalg.inv(b)
        x, y = points[i, 2:], points[j, 2:]
        action = (x @ bi @ a @ x - 2 * x @ bi @ y + y @ d @ bi @ y) / 2
        slack = points[j, 1] - points[i, 1] - action
        if slack <= 0:
            continue
        length = times[j] - times[i]
        indices = np.searchsorted(evaluation, quadrature)
        left = fs[indices] @ inv[endpoints[i]]
        right = fs[endpoints[j]] @ inv[indices]
        detleft = np.linalg.det(left[:, :2, 2:])
        detright = np.linalg.det(right[:, :2, 2:])
        assert np.min(detleft) > 0 and np.min(detright) > 0
        c2 = length / 2 * np.dot(weights, np.sqrt(detleft * detright))
        v = math.pi * slack**2 * c2 / math.sqrt(np.linalg.det(b))
        volumes[i, j] = volumes[j, i] = v
        proper_squared[i, j] = proper_squared[j, i] = 2 * length * slack
    return volumes, proper_squared


def target_tensor(profile):
    """Ground-truth Weyl tensor at the origin, in the generator's t,x,y,z frame."""
    tensor = np.zeros((4, 4, 4, 4))
    for i in range(2):
        for j in range(2):
            value = -profile(0)[i, j]
            tensor[0, i + 2, 0, j + 2] = value
            tensor[i + 2, 0, j + 2, 0] = value
            tensor[i + 2, 0, 0, j + 2] = -value
            tensor[0, i + 2, j + 2, 0] = -value
    jac = np.array(
        [[1, 0, 0, 1], [1, 0, 0, -1], [0, math.sqrt(2), 0, 0], [0, 0, math.sqrt(2), 0]]
    ) / math.sqrt(2)
    return np.einsum("ai,bj,ck,dl,abcd->ijkl", jac, jac, jac, jac, tensor)


def run():
    assert np.max(abs(np.einsum("ac,kabcd->kbd", ETA, WEYL))) == 0
    assert (
        np.max(
            abs(WEYL + WEYL.transpose(0, 1, 3, 4, 2) + WEYL.transpose(0, 1, 4, 2, 3))
        )
        == 0
    )
    for i in range(5):
        for j in range(i + 1, len(TANGENTS)):
            delta = TANGENTS[j] - TANGENTS[i]
            assert delta @ ETA @ delta < 0
    refs = TANGENTS[1:5]
    coefficients = TANGENTS[5:] @ np.linalg.inv(refs)
    design, pairs = residual_design(refs, coefficients)
    available = [
        i == j
        or (TANGENTS[j + 5] - TANGENTS[i + 5])
        @ ETA
        @ (TANGENTS[j + 5] - TANGENTS[i + 5])
        < 0
        for i, j in pairs
    ]
    design = design[available]
    singular = np.linalg.svd(design, compute_uv=False)
    assert np.linalg.matrix_rank(design) == 10
    # Verify the order-only boundary separately from the continuum oracle.
    causal = np.triu(np.ones((21, 21), dtype=bool), 1)
    ids = np.arange(0, 21, 2)
    expected = np.maximum(np.abs(ids[:, None] - ids[None, :]) - 1, 0)
    permutation = np.arange(21)[::-1]
    with patch.object(
        intrinsic_weyl, "fit_from_interval_counts", side_effect=lambda a: a
    ):
        extracted = intrinsic_weyl.fit_from_order(causal, ids)
        relabelled = intrinsic_weyl.fit_from_order(
            causal[permutation][:, permutation], 20 - ids
        )
    assert np.array_equal(extracted, expected) and np.array_equal(relabelled, expected)
    cases = []
    for name in ("flat", "constant", "rotating"):
        profile = physical_profile(name)
        rows = []
        for scale in (0.04, 0.02, 0.01, 0.005):
            volumes, tau2 = exact_interval_tables(profile, scale)
            fitted = fit_from_interval_counts(volumes)
            exact_clocks = fit_from_interval_counts(tau2**2)
            if scale == 0.02:
                rescaled = fit_from_interval_counts(17 * volumes)
                assert (
                    np.max(
                        abs(np.array(rescaled["coefficients"]) - fitted["coefficients"])
                    )
                    < 2e-9
                )
                permutation = np.r_[np.arange(5), np.arange(5, len(volumes))[::-1]]
                reordered = fit_from_interval_counts(
                    volumes[permutation][:, permutation]
                )
                assert (
                    np.max(
                        abs(
                            np.array(reordered["coefficients"]) - fitted["coefficients"]
                        )
                    )
                    < 2e-9
                )
            target = target_tensor(profile) * tau2[0, 4]
            estimated = np.einsum("k,kabcd->abcd", fitted["coefficients"], WEYL)
            tensor_error = np.linalg.norm(estimated - target)
            if name != "flat":
                tensor_error /= np.linalg.norm(target)
            row = {
                "scale": scale,
                "fit_using_continuum_volumes": fitted,
                "fit_using_exact_proper_times": exact_clocks,
                "true_electric_norm_squared": float(np.sum(target[0, 1:, 0, 1:] ** 2)),
                "relative_tensor_error_or_flat_absolute": float(tensor_error),
            }
            rows.append(row)
        cases.append({"profile": name, "rows": rows})
        errors = [r["relative_tensor_error_or_flat_absolute"] for r in rows]
        if name == "flat":
            assert max(errors) < 1e-9
        else:
            assert all(b < a for a, b in zip(errors[:-1], errors[1:], strict=True))
    return {
        "kind": (
            "Deterministic continuum-oracle feasibility diagnostic; "
            "not sampled order recovery"
        ),
        "date": "2026-09-22",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "design_rank": 10,
        "order_count_boundary_and_relabelling": (
            "passed; combinatorial test, not sampled recovery"
        ),
        "design_singular_values": singular.tolist(),
        "cases": cases,
    }


if __name__ == "__main__":
    result = run()
    result["source_sha256_lf"] = {
        p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        for p in (
            Path(__file__),
            Path(__file__).with_name("intrinsic_weyl.py"),
            ROOT / "docs/paper/paper_a/window_response/check_plane_wave_chains.py",
        )
    }
    destination = ROOT / "docs/paper/paper_a/intrinsic_weyl/continuum_diagnostic.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "written": str(destination),
                "cases": [
                    {
                        "profile": c["profile"],
                        "tensor_errors": [
                            r["relative_tensor_error_or_flat_absolute"]
                            for r in c["rows"]
                        ],
                    }
                    for c in result["cases"]
                ],
            },
            indent=2,
        )
    )
