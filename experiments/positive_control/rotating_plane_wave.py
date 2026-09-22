"""Generator-side rotating-polarization vacuum plane wave, no estimator imports.

Physical Brinkmann H(u)=R(w*(u-1/2)) diag(A,-A) R^T. Stored transverse
coordinates are y=R^T x. The exact propagator in (y, R^T x') is exp(K*du).
Only a bounded regular family is implemented: 0<=A<=.75, |w|<=.5, |du|<=1.
"""

import math

import numpy as np
from numpy.polynomial.polynomial import polyval


class RotatingWave:
    def __init__(self, amplitude, rotation_rate):
        if not (0 <= amplitude <= 0.75 and abs(rotation_rate) <= 0.5):
            raise ValueError("Outside the bounded rotating-wave generator family")
        self.amplitude, self.rotation_rate = float(amplitude), float(rotation_rate)
        rplus, rminus = amplitude - rotation_rate**2, -amplitude - rotation_rate**2
        self.coefficients = []
        for n in range(13):
            homogeneous = sum(rplus**j * rminus ** (n - j) for j in range(n + 1))
            self.coefficients.append(
                [
                    rplus**n / math.factorial(2 * n),
                    rminus**n / math.factorial(2 * n),
                    rplus**n / math.factorial(2 * n + 1),
                    rminus**n / math.factorial(2 * n + 1),
                    homogeneous / math.factorial(2 * n + 2),
                    homogeneous / math.factorial(2 * n + 3),
                ]
            )
        self.coefficients = np.array(self.coefficients)

    def blocks(self, length):
        """Return A, B/length, C, D blocks without subtracting close cosines.

        Taylor tails are <1e-24 in this bounded domain. Floating roundoff is
        handled separately by the action guard; this is not an interval library.
        """
        length = np.asarray(length, dtype=float)
        if not np.isfinite(length).all() or np.any(abs(length) > 1):
            raise ValueError("Finite propagator lengths in [-1,1] required")
        cp, cm, sp, sm, hc, hs = polyval(length**2, self.coefficients)
        w, a = self.rotation_rate, self.amplitude
        aa = np.zeros(length.shape + (2, 2))
        bb, cc, dd = (np.zeros_like(aa) for _ in range(3))
        aa[..., 0, 0] = dd[..., 0, 0] = cp
        aa[..., 1, 1] = dd[..., 1, 1] = cm
        aa[..., 0, 1], aa[..., 1, 0] = w * length * sm, -w * length * sp
        dd[..., 0, 1], dd[..., 1, 0] = w * length * sp, -w * length * sm
        bb[..., 0, 0] = sp - 2 * w * w * length**2 * hs
        bb[..., 1, 1] = sm - 2 * w * w * length**2 * hs
        bb[..., 0, 1], bb[..., 1, 0] = 2 * w * length * hc, -2 * w * length * hc
        cc[..., 0, 0], cc[..., 1, 1] = a * length * sp, -a * length * sm
        return aa, bb, cc, dd

    @staticmethod
    def inverse(matrix):
        det = (
            matrix[..., 0, 0] * matrix[..., 1, 1]
            - matrix[..., 0, 1] * matrix[..., 1, 0]
        )
        if np.any(det <= 0) or not np.isfinite(det).all():
            raise ArithmeticError("Nonpositive scaled Jacobi determinant")
        result = np.empty_like(matrix)
        result[..., 0, 0], result[..., 1, 1] = matrix[..., 1, 1], matrix[..., 0, 0]
        result[..., 0, 1], result[..., 1, 0] = -matrix[..., 0, 1], -matrix[..., 1, 0]
        return result / det[..., None, None]

    def action(self, length, left, right):
        """Exact-action formula and a conservative floating arithmetic guard."""
        length = np.asarray(length)
        if np.any(length <= 0):
            raise ValueError("Positive action lengths required")
        aa, bb, _, dd = self.blocks(length)
        inverse = self.inverse(bb)
        first, last = inverse @ aa, dd @ inverse
        terms = np.stack(
            [
                np.einsum("...i,...ij,...j->...", left, first, left),
                -2 * np.einsum("...i,...ij,...j->...", left, inverse, right),
                np.einsum("...i,...ij,...j->...", right, last, right),
            ]
        )
        value = terms.sum(axis=0) / (2 * length)
        norm = lambda m: np.max(np.sum(abs(m), axis=-1), axis=-1)  # noqa: E731
        g, b, an, dn = norm(inverse), norm(bb), norm(aa), norm(dd)
        eps, block_error = np.finfo(float).eps, 1e-12
        residual = 2 * block_error * g + 64 * eps * (1 + b * g)
        if np.any(residual >= 0.5):
            raise ArithmeticError("Jacobi inversion error cannot be controlled")
        inverse_error = g * residual / (1 - residual)
        efirst = (
            inverse_error * an
            + 2 * block_error * (g + inverse_error)
            + 64 * eps * g * an
        )
        elast = (
            inverse_error * dn
            + 2 * block_error * (g + inverse_error)
            + 64 * eps * g * dn
        )
        lp, rp = np.max(abs(left), axis=-1), np.max(abs(right), axis=-1)
        error = (
            lp * lp * efirst + 2 * lp * rp * inverse_error + rp * rp * elast
        ) / length
        error += 64 * eps * (1 + abs(terms).sum(axis=0)) / length
        return value, np.nextafter(4 * error, np.inf)

    def relations(self, left, right):
        """Broadcastable pair queries; unresolved margins stop the generator."""
        du, dv = right[..., 0] - left[..., 0], right[..., 1] - left[..., 1]
        future = du > 0
        first = np.where(future[..., None], left[..., 2:], right[..., 2:])
        last = np.where(future[..., None], right[..., 2:], left[..., 2:])
        cost, error = self.action(np.where(du == 0, 1.0, abs(du)), first, last)
        margin = -np.sign(du) * dv - cost
        guard = error + 256 * np.finfo(float).eps * (1 + abs(dv))
        near = (abs(margin) <= guard) & (du != 0)
        if near.any():
            raise ArithmeticError(f"{int(near.sum())} unresolved rotating-wave margins")
        related = (margin > guard) & (du != 0)
        return related & future, related & ~future

    def order(self, points):
        indices = np.triu_indices(len(points), 1)
        result = np.zeros((len(points), len(points)), dtype=bool)
        for start in range(0, len(indices[0]), 262144):
            i, j = (k[start : start + 262144] for k in indices)
            future, past = self.relations(points[i], points[j])
            result[i, j], result[j, i] = future, past
        return result

    def slice_matrices(self, u):
        _, left_b, _, left_d = self.blocks(u)
        right_a, right_b, _, _ = self.blocks(1 - u)
        left = (left_d @ self.inverse(left_b)) / np.asarray(u)[..., None, None]
        right = (self.inverse(right_b) @ right_a) / (1 - np.asarray(u))[..., None, None]
        q = left + right
        return left, right, (q + q.swapaxes(-1, -2)) / 2

    def sample(self, n, rng):
        """Physical-volume rejection sample in the actual curved diamond."""
        if n < 2:
            raise ValueError("At least two events required")
        # One-endpoint Poincare: slice <= pi*u*(1-u)/(1-A/2).
        ceiling = np.nextafter(np.pi / (4 * (1 - self.amplitude / 2)), np.inf)
        accepted, count = [], 0
        while count < n:
            u = rng.random(max(16, math.ceil(2.5 * (n - count))))
            u = u[u > 0]
            q = self.slice_matrices(u)[2]
            density = np.pi / np.sqrt(np.linalg.det(q))
            if np.max(density) > ceiling:
                raise ArithmeticError("Analytic rejection ceiling violated")
            kept = u[rng.random(len(u)) * ceiling < density]
            accepted.append(kept)
            count += len(kept)
        u = np.concatenate(accepted)[:n]
        left, right, q = self.slice_matrices(u)
        uniform = rng.random((n, 3))
        radius2 = uniform[:, 0] / (1 + np.sqrt(1 - uniform[:, 0]))
        angle = 2 * np.pi * uniform[:, 1]
        z = np.sqrt(2 * radius2)[:, None] * np.column_stack(
            (np.cos(angle), np.sin(angle))
        )
        y = np.linalg.solve(np.linalg.cholesky(q).swapaxes(-1, -2), z[..., None])[
            ..., 0
        ]
        lo = -0.5 + 0.5 * np.einsum("ni,nij,nj->n", y, right, y)
        hi = 0.5 - 0.5 * np.einsum("ni,nij,nj->n", y, left, y)
        if np.any(hi <= lo):
            raise ArithmeticError("Nonpositive conditional sampling width")
        return np.column_stack((u, lo + uniform[:, 2] * (hi - lo), y))

    def volumes(self, points, quadrature_order=64):
        nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
        nodes, weights = (nodes + 1) / 2, weights / 2
        output = np.zeros((len(points), len(points)))
        for i in range(len(points) - 1):
            for j in range(i + 1, len(points)):
                p, q = points[i], points[j]
                if p[0] > q[0]:
                    p, q = q, p
                du = q[0] - p[0]
                if du == 0:
                    continue
                slack = p[1] - q[1] - self.action(du, p[2:], q[2:])[0]
                if slack <= 0:
                    continue
                b0 = np.linalg.det(self.blocks(du)[1])
                bl = np.linalg.det(self.blocks(du * nodes)[1])
                br = np.linalg.det(self.blocks(du * (1 - nodes))[1])
                integral = np.dot(weights, nodes * (1 - nodes) * np.sqrt(bl * br / b0))
                output[i, j] = output[j, i] = np.pi * du**2 * slack**2 * integral
        return output

    def normal_coordinates(self, points):
        """Truth-only inverse exponential map in the parallel midpoint frame."""
        u = points[:, 0] - 0.5
        x = np.einsum("nij,nj->ni", self.inverse(self.blocks(u)[1]), points[:, 2:])
        nodes, weights = np.polynomial.legendre.leggauss(32)
        nodes, weights = (nodes + 1) / 2, weights / 2
        path = nodes[None, :, None] * np.einsum(
            "nsij,nj->nsi", self.blocks(u[:, None] * nodes)[1], x
        )
        integrand = (
            path[..., 0] ** 2
            - path[..., 1] ** 2
            - 2
            * self.rotation_rate
            * u[:, None]
            * (1 - nodes)
            * path[..., 0]
            * path[..., 1]
        )
        v = points[:, 1] + u * self.amplitude * (integrand @ weights)
        return np.column_stack(((u - v) / 2, x / np.sqrt(2), (u + v) / 2))


def selfcheck():
    from intrinsic_weyl_joint_audit import normal_coordinates, volume_table
    from scipy.integrate import solve_ivp
    from scipy.linalg import expm
    from window_response import action, diamond_moments

    j = np.array([[0.0, -1.0], [1.0, 0.0]])
    symplectic = np.block(
        [[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]]
    )
    points = np.array(
        [
            [0.0, 0.5, 0, 0],
            [0.2, 0.17, 0.04, -0.03],
            [0.5, 0.01, -0.06, 0.05],
            [0.73, -0.13, 0.03, 0.02],
            [1.0, -0.5, 0, 0],
        ]
    )
    maximum = normal_error = action_error = 0.0
    nodes, weights = np.polynomial.legendre.leggauss(64)
    nodes, weights = (nodes + 1) / 2, weights / 2
    for a in (0.0, 0.25, 0.75):
        for w in (-0.5, 0.0, 0.25, 0.5):
            wave = RotatingWave(a, w)
            k = np.block([[-w * j, np.eye(2)], [np.diag([a, -a]), -w * j]])
            for length in (-1.0, -0.25, -1e-8, 0.0, 1e-8, 0.25, 1.0):
                aa, bs, cc, dd = wave.blocks(length)
                f = np.block([[aa, length * bs], [cc, dd]])
                error = float(np.max(abs(f - expm(k * length))))
                maximum = max(maximum, error)
                assert error < 2e-14
                assert np.max(abs(f.T @ symplectic @ f - symplectic)) < 2e-14
            v = wave.volumes(points)
            assert np.max(abs(v - wave.volumes(points, 32))) < 1e-12
            density = np.pi / np.sqrt(np.linalg.det(wave.slice_matrices(nodes)[2]))
            assert abs(weights @ density - v[0, -1]) < 1e-12
            bound = np.pi * nodes * (1 - nodes) / (1 - a / 2)
            assert np.all(density <= bound * (1 + 64 * np.finfo(float).eps))

            def rotation(u, w=w):
                c, s = np.cos(w * (u - 0.5)), np.sin(w * (u - 0.5))
                return np.array([[c, -s], [s, c]])

            def profile(u, a=a, w=w):
                r = rotation(u)
                h = r @ np.diag([a, -a]) @ r.T
                return h, w * (j @ h - h @ j)

            # Independent ODE in original Brinkmann x coordinates; the production
            # propagator instead uses co-rotating y and a constant 4x4 matrix.
            def geodesic_rhs(_, state):
                u, _, x, y, du, dv, dx, dy = state
                transverse, velocity = np.array([x, y]), np.array([dx, dy])
                h, dh = profile(u)
                acceleration = du**2 * h @ transverse
                vdd = -0.5 * du**2 * (transverse @ dh @ transverse)
                vdd -= 2 * du * (velocity @ h @ transverse)
                return [du, dv, dx, dy, 0, vdd, *acceleration]

            for p, normal in zip(points, wave.normal_coordinates(points), strict=True):
                du, dv = normal[0] + normal[3], normal[3] - normal[0]
                initial = [0.5, 0, 0, 0, du, dv, *(np.sqrt(2) * normal[1:3])]
                result = solve_ivp(
                    geodesic_rhs, [0, 1], initial, rtol=1e-12, atol=1e-14
                )
                assert result.success
                expected = np.r_[p[:2], rotation(p[0]) @ p[2:]]
                error = float(np.max(abs(result.y[:4, -1] - expected)))
                normal_error = max(normal_error, error)
                assert error < 1e-11
            for p, q in zip(points[:-1], points[1:], strict=True):
                length = q[0] - p[0]
                aa, bs, _, _ = wave.blocks(length)
                initial_velocity = (
                    rotation(p[0]) @ (wave.inverse(bs) @ (q[2:] - aa @ p[2:])) / length
                )

                def transverse_rhs(u, state):
                    h = profile(u)[0]
                    xx, xp = state[:2], state[2:4]
                    return np.r_[xp, h @ xx, 0.5 * (xp @ xp + xx @ h @ xx)]

                result = solve_ivp(
                    transverse_rhs,
                    p[:1].tolist() + q[:1].tolist(),
                    np.r_[rotation(p[0]) @ p[2:], initial_velocity, 0],
                    rtol=1e-12,
                    atol=1e-14,
                )
                assert result.success
                assert np.max(abs(result.y[:2, -1] - rotation(q[0]) @ q[2:])) < 1e-11
                cost, guard = wave.action(length, p[2:], q[2:])
                error = abs(float(result.y[-1, -1]) - cost)
                action_error = max(action_error, error)
                assert error < guard
                forward, backward = wave.relations(p, q)
                assert forward == (p[1] - q[1] > result.y[-1, -1]) and not backward
                reverse_forward, reverse_backward = wave.relations(q, p)
                assert not reverse_forward and reverse_backward == forward
            if w == 0:
                assert np.max(abs(v - volume_table(points, a))) < 1e-12
                assert abs(v[0, -1] - diamond_moments(a)["volume"]) < 1e-12
                assert (
                    np.max(
                        abs(
                            wave.normal_coordinates(points)
                            - normal_coordinates(points, a)
                        )
                    )
                    < 1e-12
                )
                for p, q in zip(points[:-1], points[1:], strict=True):
                    got, guard = wave.action(q[0] - p[0], p[2:], q[2:])
                    expected, _ = action(q[0] - p[0], p[2], q[2], p[3], q[3], a)
                    assert abs(got - expected) < guard
    return {
        "propagator_max_error": maximum,
        "symplectic_identity": "passed",
        "constant_action_volume_and_normal_map": "passed",
        "volume_refinement": "passed",
        "independent_brinkmann_geodesic_roundtrip_max_error": normal_error,
        "independent_action_integral_max_error": float(action_error),
        "physical_slice_density_and_rejection_bound": "passed",
    }


if __name__ == "__main__":
    import json

    print(json.dumps(selfcheck()))
