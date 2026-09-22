"""Exact rational rank certificate; use the existing SymPy environment.

This checks an algebraic identifiability claim, not finite-sample accuracy.
"""

import hashlib
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/positive_control"))
from intrinsic_weyl import PAIRS, WEYL  # noqa: E402


def certificate():
    half, quarter = sp.Rational(1, 2), sp.Rational(1, 4)
    x = sp.Matrix(
        [
            [1, 0, 0, 0],
            [2, half, 0, 0],
            [3, 0, half, 0],
            [4, 0, 0, half],
            [5, quarter, quarter, quarter],
            [6, -quarter, quarter, 0],
            [7, quarter, -quarter, 0],
            [8, 0, quarter, -quarter],
            [9, -quarter, 0, quarter],
        ]
    )
    eta = sp.diag(-1, 1, 1, 1)
    all_points = sp.zeros(1, 4).col_join(x)
    for i, j in itertools.combinations(range(len(all_points.tolist())), 2):
        delta = all_points[j, :] - all_points[i, :]
        assert (delta * eta * delta.T)[0] < 0
    basis = [
        sp.Matrix([[int(t[i, j, k, ell]) for k, ell in PAIRS] for i, j in PAIRS])
        for t in WEYL
    ]

    def quartic(a, b):
        wedge = sp.Matrix([a[i] * b[j] - a[j] * b[i] for i, j in PAIRS])
        return sp.Matrix([(wedge.T * tensor * wedge)[0] for tensor in basis])

    references, probes = x[:4, :], x[4:, :]
    c = probes * references.inv()
    rr = [
        [quartic(references[i, :], references[j, :]) for j in range(4)]
        for i in range(4)
    ]
    rp = [[quartic(references[i, :], probes[j, :]) for i in range(4)] for j in range(5)]
    rows = []
    for i, j in itertools.combinations_with_replacement(range(5), 2):
        value = quartic(probes[i, :], probes[j, :])
        for a in range(4):
            value -= c[j, a] * rp[i][a] + c[i, a] * rp[j][a]
            for b in range(4):
                value += c[i, a] * c[j, b] * rr[a][b]
        rows.append((value / 6).T)
    design = sp.Matrix.vstack(*rows)
    pivot_rows = list(design.T.rref()[1])
    assert len(pivot_rows) == 10
    determinant = sp.factor(design[pivot_rows, :].det())
    assert determinant != 0
    # Exact obstruction when only clocks to the five references are retained.
    # This symbolic point ranges freely: adding more probes cannot add channels.
    symbols = sp.symbols("x0:4")
    point = sp.Matrix([symbols])
    cp = point * references.inv()
    star = sp.zeros(10, 1)
    for a in range(4):
        star -= 2 * cp[a] * quartic(references[a, :], point)
        for b in range(4):
            star += cp[a] * cp[b] * rr[a][b]
    polynomials = [sp.Poly(sp.expand(v / 6), *symbols) for v in star]
    monomials = sorted(set().union(*(set(p.monoms()) for p in polynomials)))
    star_map = sp.Matrix(
        [[p.coeff_monomial(m) for p in polynomials] for m in monomials]
    )
    assert star_map.rank() == 8
    return {
        "kind": "Exact rational certificate for the projected ten-component Weyl map",
        "shape": list(design.shape),
        "rank": 10,
        "pivot_rows_zero_based": pivot_rows,
        "minor_determinant": str(determinant),
        "reference_only_clock_polynomial_rank": 8,
        "reference_only_nullspace": [
            [str(v) for v in vector] for vector in star_map.nullspace()
        ],
        "configuration": [[str(v) for v in row] for row in all_points.tolist()],
        "all_pairs_timelike": True,
        "scope": (
            "Existence/generic local identifiability, "
            "not precision or landmark-selection success"
        ),
    }


if __name__ == "__main__":
    result = certificate()
    result["sympy"] = sp.__version__
    result["source_sha256_lf"] = {
        p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        for p in (
            Path(__file__),
            ROOT / "experiments/positive_control/intrinsic_weyl.py",
        )
    }
    target = Path(__file__).with_name("identifiability_certificate.json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "rank": result["rank"],
                "minor_determinant": result["minor_determinant"],
                "written": str(target),
            },
            indent=2,
        )
    )
