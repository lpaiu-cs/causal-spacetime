"""Prospective finite-order development pilot, not an accuracy certification.

Only the generator reads amplitude/coordinates. fit_order_archive reads an
archive containing causal relations and random-selection IDs, and nothing else.
"""

import argparse
import hashlib
import itertools
import json
import math
import subprocess
import time
from pathlib import Path

import numpy as np
from intrinsic_weyl_joint import (
    ETA,
    fit_joint,
    fit_order_archive,
    model,
    select_landmarks,
    selection,
)
from p14_plane_wave import PlaneWaveGeometry, Slab, causal_relation
from window_response import action, order_matrix, sample_points

ROOT = Path(__file__).resolve().parents[2]
MEASUREMENT_N = 262144
SELECTION_N = 1024
CHUNK = 4096


def relation_queries(landmarks, points, amplitude):
    """Generator-side order queries; ambiguous decisions use the existing guard."""
    du = points[None, :, 0] - landmarks[:, None, 0]
    dv = points[None, :, 1] - landmarks[:, None, 1]
    cost, error = action(
        np.maximum(abs(du), 1e-300),
        landmarks[:, None, 2],
        points[None, :, 2],
        landmarks[:, None, 3],
        points[None, :, 3],
        amplitude,
    )
    margin = -(np.sign(du) * dv + cost)
    guard = error + 256 * np.finfo(float).eps * (1 + abs(dv))
    near = (abs(margin) <= guard) | (du == 0)
    related = margin > 0
    geometry = PlaneWaveGeometry(Slab(1, 8, 16, 16), math.sqrt(amplitude))
    for i, j in np.argwhere(near):
        p, q = landmarks[i], points[j]
        if p[0] > q[0]:
            p, q = q, p
        verdict = causal_relation(geometry, p, q)
        if verdict.related is None:
            raise ArithmeticError("An order query remained ambiguous")
        related[i, j] = verdict.related
    return related & (du > 0), related & (du < 0), int(near.sum())


def generate_archive(path, amplitude, seed):
    selection_rng = np.random.default_rng([seed, 0])
    points = np.vstack(
        (
            [0, 0.5, 0, 0],
            sample_points("D", amplitude, SELECTION_N, selection_rng),
            [1, -0.5, 0, 0],
        )
    )
    # order_matrix returns the u-sorted matrix. Randomize interior labels before
    # handing it to the selector, so sorting cannot serve as hidden coordinates.
    ordered = points[np.lexsort((-points[:, 1], points[:, 0]))]
    order, escalations = order_matrix(points, amplitude)
    permutation = np.r_[
        0, selection_rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
    ]
    ordered, order = ordered[permutation], order[permutation][:, permutation]
    ids, _ = selection(order)
    landmarks = ordered[ids]
    after = np.zeros((len(ids), MEASUREMENT_N // 8), dtype=np.uint8)
    before = np.zeros_like(after)
    measurement_rng = np.random.default_rng([seed, 1])
    for start in range(0, MEASUREMENT_N, CHUNK):
        events = sample_points("D", amplitude, CHUNK, measurement_rng)
        future, past, near = relation_queries(landmarks, events, amplitude)
        after[:, start // 8 : (start + CHUNK) // 8] = np.packbits(future, axis=1)
        before[:, start // 8 : (start + CHUNK) // 8] = np.packbits(past, axis=1)
        escalations += near
    np.savez_compressed(
        path,
        selection_order=order,
        selected_ids=ids,
        after=after,
        before=before,
        measurement_n=MEASUREMENT_N,
    )
    return {"escalations": escalations, "landmarks": len(ids) - 2}


def selfcheck():
    x = np.sin(np.arange(32).reshape(8, 4)) * 0.1
    pairs = np.array(list(itertools.combinations(range(8), 2)))
    c = np.linspace(-0.2, 0.3, 10)
    _, dl, dr, _ = model(x, c, pairs)
    for point, dim in [(1, 0), (4, 2), (6, 3)]:
        delta = np.zeros_like(x)
        delta[point, dim] = 1e-6
        numeric = (model(x + delta, c, pairs)[0] - model(x - delta, c, pairs)[0]) / 2e-6
        analytic = (pairs[:, 0] == point) * dl[:, dim] + (pairs[:, 1] == point) * dr[
            :, dim
        ]
        assert np.max(abs(numeric - analytic)) < 1e-9
    # Exact flat clocks, independent of the model's curvature implementation.
    k = np.arange(48)
    z = 1 - 2 * (k + 0.5) / len(k)
    angle = k * np.pi * (3 - np.sqrt(5))
    r = 0.12 + 0.08 * np.sin(k * 2.1) ** 2
    directions = np.column_stack(
        (np.sqrt(1 - z * z) * np.cos(angle), np.sqrt(1 - z * z) * np.sin(angle), z)
    )
    x = np.zeros((50, 4))
    x[0, 0], x[-1, 0] = -0.5, 0.5
    x[1:-1, 0] = np.where(k % 2 == 0, -0.27, 0.27) + 0.025 * np.sin(k)
    x[1:-1, 1:] = r[:, None] * directions
    d = x[:, None] - x[None, :]
    tau2 = -np.einsum("ija,ab,ijb->ij", d, ETA, d)
    counts = 1e10 * np.maximum(tau2, 0) ** 2
    ids = select_landmarks(counts)
    result = fit_joint(counts[ids[:, None], ids], starts=2)
    assert np.linalg.norm(result["coefficients"]) < 1e-7
    assert result["projected_curvature_rank"] == 10
    # Compare a rectangular query to the previously validated complete predicate.
    points = np.array(
        [
            [0, 0.5, 0, 0],
            [0.2, 0.15, 0.05, 0.02],
            [0.4, 0.01, -0.03, 0.1],
            [0.7, -0.1, 0.04, -0.02],
            [1, -0.5, 0, 0],
        ]
    )
    for amplitude in (0.0, 0.25, 0.75):
        full, _ = order_matrix(points, amplitude)
        future, past, _ = relation_queries(points[[0, 2, 4]], points[[1, 3]], amplitude)
        assert np.array_equal(future, full[[0, 2, 4]][:, [1, 3]])
        assert np.array_equal(past, full.T[[0, 2, 4]][:, [1, 3]])
    return {
        "derivatives": "passed",
        "flat_clock_recovery": "passed",
        "order_queries": "matched existing predicate",
    }


def run(destination):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit protocol and code before stochastic data")
    if destination.exists():
        raise FileExistsError("Use a new directory; preserve earlier pilot data")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    check = selfcheck()
    destination.mkdir(parents=True)
    rows = []
    for case, (amplitude, repeat) in enumerate(
        itertools.product((0.0, 0.25, 0.75), range(2))
    ):
        seed = 2026092220 + case
        path = destination / f"case_{case}.npz"
        started = time.perf_counter()
        row = {
            "case": case,
            "generator_amplitude": amplitude,
            "repeat": repeat,
            "seed": seed,
            "expected_electric_eigenvalues": [-amplitude, 0, amplitude],
        }
        try:
            row["generation"] = generate_archive(path, amplitude, seed)
            # This function reloads the order-only archive and cannot read row.
            row["fit"] = fit_order_archive(path)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            row["failure"] = f"{type(error).__name__}: {error}"
        row["seconds"] = time.perf_counter() - started
        rows.append(row)
        (destination / f"case_{case}_result.json").write_text(
            json.dumps(row, indent=2) + "\n"
        )
        print(
            json.dumps(
                {k: v for k, v in row.items() if k != "fit"}
                | {"fit_eigenvalues": row.get("fit", {}).get("electric_eigenvalues")}
            ),
            flush=True,
        )
    sources = [
        Path(__file__),
        Path(__file__).with_name("intrinsic_weyl_joint.py"),
        Path(__file__).with_name("intrinsic_weyl.py"),
        Path(__file__).with_name("window_response.py"),
        Path(__file__).with_name("p14_plane_wave.py"),
        ROOT / "docs/paper/paper_a/intrinsic_weyl/joint_pilot_protocol.md",
    ]
    result = {
        "kind": "Finite-order development pilot; no accuracy or coverage verdict",
        "code_version": stamp,
        "selfcheck": check,
        "cases": rows,
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }
    (destination / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Written {destination / 'summary.json'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--fit", type=Path)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/intrinsic_weyl_joint_pilot_v1"
    )
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    elif args.fit:
        print(json.dumps(fit_order_archive(args.fit), indent=2))
    else:
        run(args.output)
