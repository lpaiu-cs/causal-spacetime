"""Exact first-order magnetic obstruction and three-slice rank certificate."""

import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/positive_control"))
from intrinsic_weyl import PAIRS, WEYL  # noqa: E402

BASIS = [
    sp.Matrix([[int(c[a, b, d, e]) for d, e in PAIRS] for a, b in PAIRS]) for c in WEYL
]


def modular_rank(matrix, prime=1000003):
    """Exact finite-field lower bound on rational rank; no floating tolerance."""
    assert sp.isprime(prime)
    work = np.array(
        [
            [int(v.p) % prime * pow(int(v.q) % prime, -1, prime) % prime for v in row]
            for row in matrix.tolist()
        ],
        dtype=np.int64,
    )
    rank = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        work[rank] = work[rank] * pow(int(work[rank, column]), -1, prime) % prime
        factors = work[rank + 1 :, column].copy()
        work[rank + 1 :] = (work[rank + 1 :] - factors[:, None] * work[rank]) % prime
        rank += 1
        if rank == len(work):
            break
    return rank


def polynomial_identity():
    t, s = sp.symbols("t s")
    x, y = sp.Matrix(sp.symbols("x1:4")), sp.Matrix(sp.symbols("y1:4"))
    b = sp.symbols("b0:5")
    magnetic = sp.Matrix(
        [[b[0], b[2], b[3]], [b[2], -b[0] + b[1], b[4]], [b[3], b[4], -b[1]]]
    )
    xx, yy = sp.Matrix([t, *x]), sp.Matrix([s, *y])
    wedge = sp.Matrix([xx[a] * yy[c] - xx[c] * yy[a] for a, c in PAIRS])
    q = sum(b[k] * (wedge.T * BASIS[k + 5] * wedge)[0] for k in range(5))
    fx, fy = x.cross(magnetic * x), y.cross(magnetic * y)
    expected = 2 * (t * x.dot(fy) + s * y.dot(fx))
    assert sp.expand(q - expected) == 0
    assert sp.expand(x.dot(fx)) == 0 and sp.expand(y.dot(fy)) == 0
    deformation = 2 * (x - y).dot(s * fx / 3 - t * fy / 3)
    assert sp.expand(deformation + q / 3) == 0
    return {
        "magnetic_quartic_identity": "verified exactly",
        "two_slice_coordinate_absorption": "verified exactly for all five B components",
    }


def rank_certificate(times):
    tenth = sp.Rational(1, 10)
    points = [sp.Matrix([-sp.Rational(1, 2), 0, 0, 0])]
    for time in times:
        for k, signs in enumerate(itertools.product((-1, 1), repeat=3)):
            # Construct rationals from integer formulas, not binary floats.
            spatial = [
                sp.Rational(signs[0] * (7 + k), 100),
                sp.Rational(signs[1] * (40 - k), 500),
                sp.Rational(signs[2] * (60 + 3 * k), 1000),
            ]
            points.append(sp.Matrix([time * tenth, *spatial]))
    points.append(sp.Matrix([sp.Rational(1, 2), 0, 0, 0]))
    eta = sp.diag(-1, 1, 1, 1)
    rows = []
    m = len(points) - 2
    for i, j in itertools.combinations(range(len(points)), 2):
        if i == 0 and j == len(points) - 1:
            continue
        x, y = points[i], points[j]
        delta = x - y
        if (delta.T * eta * delta)[0] >= 0:
            continue
        geometry = sp.zeros(1, 4 * m)
        for index, sign in ((i, 1), (j, -1)):
            if 0 < index < len(points) - 1:
                geometry[0, 4 * (index - 1) : 4 * index] = (2 * sign * eta * delta).T
        wedge = sp.Matrix([x[a] * y[b] - x[b] * y[a] for a, b in PAIRS])
        curvature = sp.Matrix([[-(wedge.T * c * wedge)[0] / 3 for c in BASIS]])
        rows.append(geometry.row_join(curvature))
    joint = sp.Matrix.vstack(*rows)
    geometry_rank = modular_rank(joint[:, : 4 * m])
    joint_rank = modular_rank(joint)
    assert sp.Matrix([list(p[1:]) for p in points[1:-1]]).rank() == 3
    # Three independent rotations give the geometric upper bound. In two
    # slices the symbolic identity adds five independent magnetic null modes.
    geometry_upper = 4 * m - 3
    joint_upper = geometry_upper + (5 if len(times) == 2 else 10)
    assert geometry_rank == geometry_upper
    assert joint_rank == joint_upper
    return {
        "time_slices": [str(t * tenth) for t in times],
        "interior_landmarks": m,
        "timelike_clocks": len(rows),
        "geometry_rank": geometry_rank,
        "joint_rank": joint_rank,
        "projected_weyl_rank": joint_rank - geometry_rank,
        "exact_rank_certificate": {
            "prime": 1000003,
            "geometry_modular_lower_bound": geometry_rank,
            "geometry_rotation_upper_bound": geometry_upper,
            "joint_modular_lower_bound": joint_rank,
            "joint_analytic_upper_bound": joint_upper,
        },
        "points": [[str(v) for v in point] for point in points],
    }


if __name__ == "__main__":
    result = polynomial_identity()
    result["two_slices"] = rank_certificate((-3, 3))
    result["three_slices"] = rank_certificate((-3, 0, 3))
    assert result["two_slices"]["projected_weyl_rank"] == 5
    assert result["three_slices"]["projected_weyl_rank"] == 10
    result["scope"] = (
        "First-order local-vacuum interval-clock model with unknown positions; "
        "not a full-order impossibility theorem or finite-sample precision guarantee"
    )
    result["source_sha256_lf"] = {
        p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        for p in (
            Path(__file__),
            ROOT / "experiments/positive_control/intrinsic_weyl.py",
        )
    }
    target = Path(__file__).with_name("two_slice_obstruction.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "two_slice_rank": result["two_slices"]["projected_weyl_rank"],
                "three_slice_rank": result["three_slices"]["projected_weyl_rank"],
                "written": str(target),
            }
        )
    )
