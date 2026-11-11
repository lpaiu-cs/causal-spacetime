"""Numerical check of the orders used in appendix_convergence.tex.

Deterministic, CPU only, a few seconds. Exact rotating vacuum plane waves
stand in for a general vacuum spacetime, and exact interval volumes for
counts. This checks orders in one family and one landmark configuration and
an exact rank for one rational configuration, and it evaluates the constants
of the proof quoted in the appendix. It proves nothing about the landmark
selector, the optimizer, finite counts or confidence coverage.

    .venv/Scripts/python.exe docs/paper/weyl_note/check_convergence_orders.py
"""

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "experiments/positive_control"))
from intrinsic_weyl import ETA, PAIRS, WEYL, gradient_basis  # noqa: E402
from intrinsic_weyl_joint import (  # noqa: E402
    BIVECTOR,
    TIDAL_GRAM,
    TIDAL_ROOT,
    axis_signs,
    fit_joint,
    model,
)
from intrinsic_weyl_square import correction  # noqa: E402
from rotating_plane_wave import RotatingWave  # noqa: E402

AMPLITUDE, ROTATION = 0.75, 0.5
SCALES = (1.0, 0.5, 0.25, 0.125)
PRIME = 2147483647


def configuration():
    """48 landmarks in three time bands, as in intrinsic_weyl_joint.selfcheck_jet.

    Returns Brinkmann coordinates (u, v, transverse), the clock mask (pairs
    whose flat squared interval is below -0.01, the same set at every scale)
    and the flat reference coordinates in units of T.
    """
    k = np.arange(48)
    z = 1 - 2 * (k + 0.5) / len(k)
    phi = k * np.pi * (3 - np.sqrt(5))
    directions = np.column_stack(
        (np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z)
    )
    flat = np.zeros((50, 4))
    flat[0, 0], flat[-1, 0] = -0.5, 0.5
    flat[1:-1, 0] = 0.27 * ((k % 3) - 1) + 0.015 * np.sin(k)
    flat[1:-1, 1:] = (0.12 + 0.08 * np.sin(k * 2.1) ** 2)[:, None] * directions
    delta = flat[:, None] - flat[None]
    mask = np.einsum("ija,ab,ijb->ij", delta, ETA, delta) < -0.01
    mask[0, -1] = mask[-1, 0] = False
    t, x, y, zz = flat.T
    brinkmann = np.column_stack((0.5 + t + zz, zz - t, np.sqrt(2) * x, np.sqrt(2) * y))
    return brinkmann, mask, flat


def adapted_frame(points):
    """Spatial rotation to the gauge of Sect. 3: landmark 1 on x, 2 in the xy plane."""
    first = points[1, 1:] / np.linalg.norm(points[1, 1:])
    second = points[2, 1:] - (points[2, 1:] @ first) * first
    second /= np.linalg.norm(second)
    frame = np.eye(4)
    frame[1:, 1:] = np.column_stack((first, second, np.cross(first, second)))
    return frame


def free_coordinates(count):
    """Landmark coordinates left free by the gauge; the two tips are fixed."""
    return [
        (i, a)
        for i in range(1, count - 1)
        for a in range(4)
        if not (i == 1 and a >= 2) and not (i == 2 and a == 3)
    ]


def flat_jacobian(points, pairs):
    """J(xi, 0, 0) of the implemented model for gauge-fixed points, in floats."""
    free = free_coordinates(len(points))
    _, dl, dr, dc = model(points, np.zeros(34), pairs)
    jacobian = np.zeros((len(pairs), len(free) + 34))
    for n, (i, a) in enumerate(free):
        jacobian[:, n] = (pairs[:, 0] == i) * dl[:, a] + (pairs[:, 1] == i) * dr[:, a]
    jacobian[:, len(free) :] = dc
    return jacobian


def proof_constants(jacobian):
    """mu and a lower bound for the constant C_2 given by the proof of Proposition 1.

    The proof takes lambda = w_- mu^2 / 4. With M_1 >= ||J||, M_2 >= 2 (the flat
    part of the model) and w_- <= w_+, Lemma 5 gives
    epsilon_1 <= mu^4 / (3072 ||J||^2), and C_2 = 4 |P| / epsilon_1^2.
    """
    singular = np.linalg.svd(jacobian, compute_uv=False)
    mu, norm = float(singular[-1]), float(singular[0])
    epsilon = mu**4 / (3072 * norm**2)
    return {
        "mu_smallest_singular_value": mu,
        "largest_singular_value": norm,
        "log10_lower_bound_of_C2": float(np.log10(4 * len(jacobian) / epsilon**2)),
    }


def truth(wave, brinkmann):
    """True gauge-frame coordinates and the 34 coefficients (T^2 C, T^3 grad C)."""
    normal = wave.normal_coordinates(brinkmann)
    k = np.array([1, 0, 0, 1]) / np.sqrt(2)
    transverse = np.eye(4)[1:3]

    def riemann(h):
        out = np.zeros((4, 4, 4, 4))
        for i, j in itertools.product(range(2), repeat=2):
            left = np.outer(k, transverse[i]) - np.outer(transverse[i], k)
            right = np.outer(k, transverse[j]) - np.outer(transverse[j], k)
            out -= h[i, j] * np.einsum("ab,cd->abcd", left, right)
        return out

    a, w = wave.amplitude, wave.rotation_rate
    curvature = 2 * riemann(np.diag([a, -a]))
    slope = riemann(2 * w * a * np.array([[0, 1.0], [1.0, 0]]))
    gradient = 2 * np.sqrt(2) * np.einsum("e,abcd->eabcd", k, slope)
    f = adapted_frame(normal)
    curvature = np.einsum("ia,jb,kc,ld,ijkl->abcd", f, f, f, f, curvature)
    gradient = np.einsum("ne,ia,jb,kc,ld,nijkl->eabcd", f, f, f, f, f, gradient)
    basis = WEYL.reshape(10, -1).T
    c = np.linalg.lstsq(basis, curvature.ravel(), rcond=None)[0]
    g = np.array([np.linalg.lstsq(basis, s.ravel(), rcond=None)[0] for s in gradient])
    b = gradient_basis().T @ g.ravel()
    assert np.max(abs(basis @ c - curvature.ravel())) < 1e-12
    assert np.max(abs(gradient_basis() @ b - g.ravel())) < 1e-12  # vacuum Bianchi
    return normal @ f, np.r_[c, b]


def orders():
    """Halving the window: model and volume residuals by 16, relative bias by 4."""
    brinkmann, mask, flat = configuration()
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    earlier = brinkmann[pairs[:, 0], 0] < brinkmann[pairs[:, 1], 0]
    pairs = np.where(earlier[:, None], pairs, pairs[:, ::-1])
    left, right = brinkmann[pairs[:, 0]], brinkmann[pairs[:, 1]]
    rows = []
    for scale in SCALES:
        wave = RotatingWave(AMPLITUDE * scale**2, ROTATION * scale)
        points, theta = truth(wave, brinkmann)
        du = right[:, 0] - left[:, 0]
        slack = left[:, 1] - right[:, 1] - wave.action(du, left[:, 2:], right[:, 2:])[0]
        interval = -du * slack  # signed squared interval over T^2 = 2
        volumes = wave.volumes(brinkmann)
        assert np.all(volumes[mask] > 0) and np.all(slack > 0)
        clock = -np.sqrt(volumes[pairs[:, 0], pairs[:, 1]] / volumes[0, -1])
        row = {
            "scale": scale,
            "world_function_residual": float(
                np.max(abs(interval - model(points, theta, pairs)[0]))
            ),
            "relative_volume_correction": float(
                np.max(abs((clock / interval) ** 2 - 1))
            ),
        }
        for name, square in (("fit_error", False), ("fit_error_sigma2", True)):
            fit = fit_joint(
                1e8 * volumes / volumes[0, -1],
                clock_mask=mask,
                starts=2,
                model_order=5,
                square_correction=square,
            )
            assert fit["optimizer_success"]
            assert fit["nuisance_rank"] == fit["nuisance_parameters"]
            inferred = np.array(fit["coordinates_inferred_from_counts"])
            slope = np.array(fit["covariant_gradient_components"])
            errors = []
            # The x and y signs follow the gauge. The z reflection is not
            # identifiable from the order (Sect. 2), so both images are compared.
            for mirror in (1, -1):
                signs = np.array(
                    [np.sign(inferred[1, 1]), np.sign(inferred[2, 2]), mirror]
                )
                weyl = axis_signs(signs) * np.array(fit["coefficients"]) - theta[:10]
                reflected = slope * np.r_[1, signs][:, None] * axis_signs(signs)
                gradient = gradient_basis().T @ reflected.ravel() - theta[10:]
                difference = np.r_[
                    (inferred * np.r_[1, signs] - points).ravel(), weyl, gradient
                ]
                errors.append(
                    (
                        np.linalg.norm(TIDAL_ROOT @ weyl),
                        np.linalg.norm(difference),
                        np.linalg.norm(gradient),
                    )
                )
            tidal, parameter, gradient = min(errors)
            row[name + "_relative"] = float(tidal / (2 * wave.amplitude))
            row[name.replace("fit", "parameter")] = float(parameter)
            row[name.replace("fit", "gradient")] = float(gradient)
        row["true_gradient_norm"] = float(np.linalg.norm(theta[10:]))
        rows.append(row)
    ratios = {
        key: [rows[i][key] / rows[i + 1][key] for i in range(len(rows) - 1)]
        for key in rows[0]
        if key != "scale"
    }
    for key, expected, checked in (
        ("world_function_residual", 16, 2),
        ("relative_volume_correction", 16, 2),
        ("fit_error_relative", 4, 2),
        ("fit_error_sigma2_relative", 4, 2),
        ("parameter_error", 16, 1),  # |theta_hat - theta*|, all coordinates, c and b
        ("parameter_error_sigma2", 16, 1),
    ):
        last = ratios[key][-checked:]
        assert all(abs(r / expected - 1) < 0.1 for r in last), (key, ratios)
    return {
        "clocks": len(pairs),
        "rows": rows,
        "ratios_per_halving": ratios,
        "constants": proof_constants(flat_jacobian(flat @ adapted_frame(flat), pairs)),
    }


def rank_mod_prime(matrix):
    """Rank over the field of PRIME elements; a lower bound for the rational rank."""
    a = np.array(matrix, dtype=np.int64) % PRIME
    rank = 0
    for column in range(a.shape[1]):
        nonzero = np.flatnonzero(a[rank:, column])
        if not len(nonzero):
            continue
        pivot = rank + nonzero[0]
        a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, column]), -1, PRIME) % PRIME
        others = np.flatnonzero(a[:, column])
        others = others[others != rank]
        a[others] = (a[others] - a[others, column][:, None] * a[rank]) % PRIME
        rank += 1
        if rank == a.shape[0]:
            break
    return rank


def rank_certificate(per_band=5):
    """Exact full column rank of J(xi, 0, 0) for one rational configuration.

    Three time bands; clocks between different bands and to both tips. The 24
    gradient components are the kernel of the integer Bianchi rows appended
    below, so no floating-point basis enters the exact part.
    """
    points = [(Fraction(-1, 2), 0, 0, 0)]
    for n in range(1, 3 * per_band + 1):
        t = 300 * ((n - 1) // per_band - 1) + (7 * n * n + 3 * n) % 41 - 20
        x = (37 * n + 11 * n * n) % 121 - 60
        y = (53 * n + 5 * n**3) % 121 - 60
        z = (71 * n + 13 * n * n) % 121 - 60
        if n == 1:  # gauge of Sect. 3: first landmark on the x axis
            x, y, z = abs(x) + 1, 0, 0
        if n == 2:  # second landmark in the xy plane
            y, z = abs(y) + 1, 0
        points.append(tuple(Fraction(v, 1000) for v in (t, x, y, z)))
    points.append((Fraction(1, 2), 0, 0, 0))
    count = len(points)
    band = [-1] + [(n - 1) // per_band for n in range(1, count - 1)] + [3]
    pairs = [
        (i, j)
        for i, j in itertools.combinations(range(count), 2)
        if band[i] != band[j] and (i, j) != (0, count - 1)
    ]
    column = {key: n for n, key in enumerate(free_coordinates(count))}
    size = len(column)
    eta = (-1, 1, 1, 1)
    bivector = np.rint(BIVECTOR).astype(int)
    assert np.array_equal(bivector, BIVECTOR)
    rows = []
    for i, j in pairs:
        x, y = points[i], points[j]
        delta = [x[a] - y[a] for a in range(4)]
        assert (
            delta[0] < 0 and sum(e * v * v for e, v in zip(eta, delta, strict=True)) < 0
        )
        row = [Fraction(0)] * (size + 50)
        for a in range(4):
            if (i, a) in column:
                row[column[i, a]] = 2 * eta[a] * delta[a]
            if (j, a) in column:
                row[column[j, a]] = -2 * eta[a] * delta[a]
        wedge = [x[a] * y[b] - x[b] * y[a] for a, b in PAIRS]
        for k in range(10):
            quartic = sum(
                wedge[p] * int(bivector[k, p, q]) * wedge[q]
                for p in range(6)
                for q in range(6)
            )
            row[size + k] = -quartic / 3
            for e in range(4):
                row[size + 10 + 10 * e + k] = -(x[e] + y[e]) * quartic / 12
        rows.append(row)
    clocks = len(rows)
    for a, b, c in itertools.combinations(range(4), 3):
        for d, e in PAIRS:
            row = [Fraction(0)] * (size + 50)
            for k in range(10):
                row[size + 10 + 10 * a + k] += int(WEYL[k, b, c, d, e])
                row[size + 10 + 10 * b + k] += int(WEYL[k, c, a, d, e])
                row[size + 10 + 10 * c + k] += int(WEYL[k, a, b, d, e])
            rows.append(row)
    residues = [
        [v.numerator * pow(v.denominator, -1, PRIME) % PRIME for v in row]
        for row in rows
    ]
    rank = rank_mod_prime(residues)
    assert rank == size + 50, (rank, size + 50)
    # Tie the exact rows to the implemented model: same Jacobian in floats.
    exact = np.array([[float(v) for v in row] for row in rows[:clocks]])
    exact = np.column_stack(
        (exact[:, : size + 10], exact[:, size + 10 :] @ gradient_basis())
    )
    floats = np.array([[float(v) for v in p] for p in points])
    implemented = flat_jacobian(floats, np.array(pairs))
    assert np.max(abs(exact - implemented)) < 1e-12
    return {
        "landmarks": count - 2,
        "clocks": clocks,
        "parameters": size + 34,
        "rank_modulo_prime_with_bianchi_rows": rank,
        "full_column_rank": True,
        "prime": PRIME,
        "constants": proof_constants(implemented),
    }


def volume_expansion():
    """The half-interval computation of Lemma 2, done by exact quadrature.

    (i) Without the vacuum condition its term of relative order tau^2 is
    Myrheim's formula, -R/180 + R_00/30, as stated after Lemma 2.
    (ii) In vacuum its term of relative order tau^4, built from the
    curvature-squared term of the model and the volume element
    1 - tr(E^2)/180, is (17 E^2 + 3 B^2)/10080, the volume-clock coefficient
    of Eq. (4). The proof of Proposition 1 uses neither coefficient.
    """
    radius, radial = np.polynomial.legendre.leggauss(6)
    radius, radial = (radius + 1) / 2, radial / 2
    cosine, polar = np.polynomial.legendre.leggauss(6)
    angle = 2 * np.pi * np.arange(12) / 12
    r, c, a = np.meshgrid(radius, cosine, angle, indexing="ij")
    s = np.sqrt(1 - c * c)
    nu = np.stack((np.ones_like(r), r * s * np.cos(a), r * s * np.sin(a), r * c), -1)
    nu = nu.reshape(-1, 4)  # directions (1, vec nu) with vec nu in the unit ball
    weight = (radial * radius**2)[:, None, None] * polar[None, :, None] / 12 + 0 * r
    weight = weight.ravel() * 3 / 2  # average over the unit ball
    assert abs(weight.sum() - 1) < 1e-13
    tip = np.array([1.0, 0, 0, 0])
    pairs = np.column_stack((np.arange(len(nu)), np.full(len(nu), len(nu))))
    signs = np.diag(ETA)
    rng = np.random.default_rng(2026)
    worst = 0.0
    for _ in range(4):
        # (i) generic algebraic curvature tensor: a Kulkarni-Nomizu product
        h, k = (x + x.T for x in rng.normal(size=(2, 4, 4)))
        riemann = (
            np.einsum("ac,bd->abcd", h, k)
            + np.einsum("bd,ac->abcd", h, k)
            - np.einsum("ad,bc->abcd", h, k)
            - np.einsum("bc,ad->abcd", h, k)
        )
        ricci = np.einsum("ac,abcd->bd", ETA, riemann)
        tidal = np.einsum("abcd,na,b,nc,d->n", riemann, nu, tip, nu, tip)
        second = weight @ (tidal / 3 - np.einsum("ab,na,nb->n", ricci, nu, nu) / 36)
        myrheim = -np.trace(ETA @ ricci) / 180 + ricci[0, 0] / 30
        worst = max(worst, abs(second - myrheim))
        # (ii) generic Weyl tensor, electric and magnetic parts
        coefficients = rng.normal(size=10)
        weyl = np.einsum("k,kabcd->abcd", coefficients, WEYL)
        tidal = np.einsum("abcd,na,b,nc,d->n", weyl, nu, tip, nu, tip)
        square = correction(
            np.vstack((nu / 2, tip)), coefficients, pairs, volume_clock=False
        )[0]
        jacobi = np.einsum("paqb,na,nb->npq", weyl, nu, nu)
        trace = np.einsum("p,q,npq,nqp->n", signs, signs, jacobi, jacobi)
        fourth = weight @ (7 * tidal**2 / 72 - 4 * square - trace / 5760)
        norm = TIDAL_GRAM[:5, :5]
        electric = coefficients[:5] @ norm @ coefficients[:5]
        magnetic = coefficients[5:] @ norm @ coefficients[5:]
        worst = max(worst, abs(fourth - (17 * electric + 3 * magnetic) / 10080))
    assert worst < 1e-12, worst
    return {
        "myrheim_coefficients_reproduced": True,
        "vacuum_coefficient_17_3_over_10080_reproduced": True,
        "largest_discrepancy": float(worst),
    }


if __name__ == "__main__":
    print(
        json.dumps(
            {
                "orders": orders(),
                "rank": rank_certificate(),
                "volume": volume_expansion(),
            },
            indent=2,
        )
    )
