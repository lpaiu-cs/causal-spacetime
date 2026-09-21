"""Run the committed G1 checks; no pilot or confirmation posets are generated."""

from __future__ import annotations

import json
import math
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import scipy
import window_response as wr
from p14_plane_wave import PlaneWaveGeometry, Slab, causal_relation, transverse_cost
from scipy.integrate import quad, solve_bvp
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / "docs/paper/paper_a/window_response_g1_results.json"


def clean_stamp():
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("commit first: the stage requires a clean worktree")
    return (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )


def ode_action(s, xp, xq, yp, yq, amplitude):
    """Solve the boundary-value geodesic ODE, then independently integrate action."""
    grid = np.linspace(0, 1, 33)
    guess = np.vstack(
        (
            xp + (xq - xp) * grid,
            yp + (yq - yp) * grid,
            np.full(33, xq - xp),
            np.full(33, yq - yp),
        )
    )

    def ode(_u, state):
        return np.vstack(
            (
                state[2],
                state[3],
                amplitude * s * s * state[0],
                -amplitude * s * s * state[1],
            )
        )

    def boundary(left, right):
        return np.array([left[0] - xp, right[0] - xq, left[1] - yp, right[1] - yq])

    result = solve_bvp(
        ode, boundary, grid, guess, tol=1e-10, bc_tol=1e-12, max_nodes=10000
    )
    if not result.success:
        raise ArithmeticError(result.message)

    def integrand(u):
        x, y, vx, vy = result.sol(u)
        return ((vx * vx + vy * vy) / s + amplitude * s * (x * x - y * y)) / 2

    cost, error = quad(integrand, 0, 1, epsabs=1e-10, epsrel=1e-10)
    return float(cost), float(error), float(np.max(result.rms_residuals))


def main():
    stamp = clean_stamp()
    started = time.perf_counter()
    rng = np.random.default_rng([20260921, 6101, 10])
    report = {
        "kind": "G1 instrument and analytic-coefficient validation",
        "code_version": stamp,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "coefficient_checks": {},
        "predicate_checks": [],
        "ode_checks": [],
        "diamond_checks": [],
        "near_cone_checks": [],
    }
    for region in wr.BOXES:
        exact = wr.box_coefficients(region)
        rows = []
        for order in (24, 40):
            numerical = wr.coefficient_quadrature(region, order)
            difference = float(np.max(abs(exact - numerical)))
            assert difference < 2e-13
            rows.append(
                {
                    "order": order,
                    "values": numerical.tolist(),
                    "maximum_absolute_difference": difference,
                }
            )
        report["coefficient_checks"][region] = {
            "closed_form": exact.tolist(),
            "quadratures": rows,
        }
    for amplitude in (0, 1 / 64, 1 / 4, 1, 2):
        geometry = PlaneWaveGeometry(Slab(1, 8, 16, 16), math.sqrt(amplitude))
        pairs = wr.sample_points("R+", amplitude, 512, rng).reshape(256, 2, 4)
        differences = []
        for pair in pairs:
            p, q = pair[np.argsort(pair[:, 0])]
            cost, error = wr.action(q[0] - p[0], p[2], q[2], p[3], q[3], amplitude)
            legacy = transverse_cost(
                q[0] - p[0], p[2], q[2], p[3], q[3], math.sqrt(amplitude)
            )
            differences.append(float(abs(cost - legacy)))
            assert abs(cost - legacy) <= 2 * error
            matrix, _ = wr.order_matrix(pair, amplitude)
            assert bool(matrix[0, 1]) == causal_relation(geometry, p, q).related
        report["predicate_checks"].append(
            {
                "A": amplitude,
                "pairs": len(pairs),
                "max_action_difference": max(differences),
                "classification_mismatches": 0,
            }
        )
        p = np.array([0.1, 0, 0.3, 0.8])
        s = 0.7
        c = transverse_cost(s, 0.3, -0.2, 0.8, -0.6, math.sqrt(amplitude))
        for shift in (-1e-9, -1e-14, 1e-14, 1e-9):
            q = np.array([0.8, -c + shift, -0.2, -0.6])
            reference = causal_relation(geometry, p, q)
            if reference.related is None:
                try:
                    wr.order_matrix(np.array([p, q]), amplitude)
                except ArithmeticError:
                    pass
                else:
                    raise AssertionError("unresolved pair was not stopped")
            else:
                assert (
                    bool(wr.order_matrix(np.array([p, q]), amplitude)[0][0, 1])
                    == reference.related
                )
            report["near_cone_checks"].append(
                {"A": amplitude, "shift": shift, "reference": reference.related}
            )
    for amplitude in (0, 0.25, 1, 2):
        rows = []
        for _ in range(24):
            s = rng.uniform(0.02, 1)
            xp, xq, yp, yq = rng.uniform(-3, 3, 4)
            cost, quad_error, residual = ode_action(s, xp, xq, yp, yq, amplitude)
            reference = transverse_cost(s, xp, xq, yp, yq, math.sqrt(amplitude))
            scaled_difference = abs(cost - reference) / (1 + abs(reference))
            assert scaled_difference < 1e-8
            rows.append([float(scaled_difference), quad_error, residual])
        report["ode_checks"].append(
            {
                "A": amplitude,
                "pairs": 24,
                "max_scaled_action_difference": max(r[0] for r in rows),
                "max_quad_error_estimate": max(r[1] for r in rows),
                "max_ode_rms_residual": max(r[2] for r in rows),
            }
        )
    for amplitude in (0, 1, 2):
        pts = wr.sample_points("D", amplitude, 4096, rng)
        geometry = PlaneWaveGeometry(Slab(1, 8, 16, 16), math.sqrt(amplitude))
        for point in pts:
            assert causal_relation(geometry, np.array([0, 0.5, 0, 0]), point).related
            assert causal_relation(geometry, point, np.array([1, -0.5, 0, 0])).related
        volume, volume_error = quad(
            lambda u, amplitude=amplitude: wr.diamond_slice(u, amplitude),
            0,
            1,
            epsabs=1e-12,
        )
        moments = {}
        for name, values, integrand in (
            ("u_squared", pts[:, 0] ** 2, lambda u: u * u),
            (
                "x_squared",
                pts[:, 2] ** 2,
                lambda u, amplitude=amplitude: 1 / (6 * wr.diamond_b(u, amplitude)[0]),
            ),
            (
                "y_squared",
                pts[:, 3] ** 2,
                lambda u, amplitude=amplitude: 1 / (6 * wr.diamond_b(u, amplitude)[1]),
            ),
        ):
            truth = quad(
                lambda u, integrand=integrand, amplitude=amplitude: (
                    integrand(u) * wr.diamond_slice(u, amplitude)
                ),
                0,
                1,
            )[0]
            truth /= volume
            se = float(values.std(ddof=1) / math.sqrt(len(values)))
            mean = float(values.mean())
            half = float(t.ppf(0.975, len(values) - 1) * se)
            moments[name] = {
                "analytic_integral": truth,
                "sample_mean": mean,
                "CI95": [mean - half, mean + half],
                "standard_error": se,
            }
        report["diamond_checks"].append(
            {
                "A": amplitude,
                "points": len(pts),
                "membership_failures": 0,
                "volume": volume,
                "volume_quad_error_estimate": volume_error,
                "moments_diagnostic_only": moments,
            }
        )
    report["elapsed_seconds"] = time.perf_counter() - started
    report["instrument_checks"] = "PASS"
    report["scope"] = (
        "Deterministic checks passed; sampler moments are diagnostics. "
        "No power or physical confirmation verdict."
    )
    RESULT.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "result": str(RESULT),
                "stamp": stamp,
                "seconds": report["elapsed_seconds"],
                "instrument_checks": "PASS",
            }
        )
    )


if __name__ == "__main__":
    main()
