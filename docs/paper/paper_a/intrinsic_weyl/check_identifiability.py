"""Exact rational rank certificate; use the existing SymPy environment.

This checks an algebraic identifiability claim, not finite-sample accuracy.
"""

import hashlib
import itertools
import json
import subprocess
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/positive_control"))
from intrinsic_weyl import PAIRS, WEYL  # noqa: E402


def reference_only_map(references, quartic, gradient=None):
    """Polynomial coefficients of the diagonal Schur residual, not noisy data."""
    symbols = sp.symbols("x0:4")
    point = sp.Matrix([symbols])
    cp = point * references.inv()

    def response(a, b):
        value = quartic(a, b)
        if gradient is None:
            return value
        derivative = sp.Matrix(
            [(a[e] + b[e]) * value[k] / 4 for e in range(4) for k in range(10)]
        )
        return value.col_join(gradient.T * derivative)

    star = sp.zeros(10 if gradient is None else 34, 1)
    for a in range(4):
        star -= 2 * cp[a] * response(references[a, :], point)
        for b in range(4):
            star += cp[a] * cp[b] * response(references[a, :], references[b, :])
    polynomials = [sp.Poly(sp.expand(v / 6), *symbols) for v in star]
    monomials = sorted(set().union(*(set(p.monoms()) for p in polynomials)))
    matrix = sp.Matrix([[p.coeff_monomial(m) for p in polynomials] for m in monomials])
    return matrix, monomials


def skew_reference_certificate(quartic):
    """A different admissible reference geometry retains all 34 vacuum jet channels."""
    half, quarter, eighth = sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 8)
    references = sp.Matrix(
        [
            [1, 0, 0, 0],
            [2, half, 0, 0],
            [3, quarter, half, 0],
            [4, eighth, quarter, half],
        ]
    )
    points = sp.zeros(1, 4).col_join(references)
    eta = sp.diag(-1, 1, 1, 1)
    for i, j in itertools.combinations(range(5), 2):
        delta = points[j, :] - points[i, :]
        assert (delta * eta * delta.T)[0] < 0

    rows = []
    for a, b, c in itertools.combinations(range(4), 3):
        for d, e in PAIRS:
            row = [0] * 40
            for k in range(10):
                row[10 * a + k] += int(WEYL[k, b, c, d, e])
                row[10 * b + k] += int(WEYL[k, c, a, d, e])
                row[10 * c + k] += int(WEYL[k, a, b, d, e])
            rows.append(row)
    constraint = sp.Matrix(rows)
    gradient = sp.Matrix.hstack(*constraint.nullspace())
    assert gradient.shape == (40, 24)
    assert constraint * gradient == sp.zeros(len(rows), 24)
    matrix, monomials = reference_only_map(references, quartic, gradient)
    prime = 1000003
    assert sp.isprime(prime)
    residues = matrix.applyfunc(
        lambda v: int(sp.numer(v)) * pow(int(sp.denom(v)), -1, prime) % prime
    )
    domain = DomainMatrix.from_Matrix(residues).convert_to(sp.GF(prime))
    pivots = list(domain.transpose().rref()[1])
    rank = len(pivots)
    assert rank == 34
    determinant = int(domain.extract(pivots, list(range(34))).det())
    assert determinant == 438280
    # Full column rank already implies independence of the 24 nuisance columns.
    curvature, _ = reference_only_map(references, quartic)
    curvature_pivots = list(curvature.T.rref()[1])
    curvature_determinant = sp.factor(curvature[curvature_pivots, :].det())
    assert curvature_determinant == -sp.Rational(9943115, 3359232)
    return {
        "reference_vectors": [[str(v) for v in row] for row in references.tolist()],
        "references_and_origin_pairwise_timelike": True,
        "curvature_only_rank": 10,
        "curvature_only_minor_determinant": str(curvature_determinant),
        "jet_polynomial_shape": list(matrix.shape),
        "joint_rank": rank,
        "gradient_rank": 24,
        "weyl_rank_after_gradient": 10,
        "prime": prime,
        "joint_minor_rows_zero_based": pivots,
        "joint_minor_monomials": [list(monomials[i]) for i in pivots],
        "joint_minor_determinant_mod_prime": determinant,
        "scope": (
            "Generic local identifiability of the linearized fifth-degree vacuum jet "
            "from diagonal Schur residuals; no finite-count accuracy, conditioning, "
            "automatic landmark selection, or computational gain is established"
        ),
    }


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
    star_map, _ = reference_only_map(references, quartic)
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
        "skew_reference_only": skew_reference_certificate(quartic),
        "configuration": [[str(v) for v in row] for row in all_points.tolist()],
        "all_pairs_timelike": True,
        "scope": (
            "Existence/generic local identifiability, "
            "not precision or landmark-selection success"
        ),
    }


if __name__ == "__main__":
    result = certificate()
    result["code_version"] = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
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
