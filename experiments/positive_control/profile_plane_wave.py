"""Generator-side linearly polarized vacuum plane wave with a u-dependent profile.

Physical Brinkmann H(u) = a(u) diag(1, -1) on the chart, tips and sign
conventions of rotating_plane_wave:

    ds^2 = 2 du dv + a(u) (x^2 - y^2) du^2 + dx^2 + dy^2,
    tips (0, 1/2, 0, 0) and (1, -1/2, 0, 0).

Any profile a(u) is an exact vacuum solution. The transverse Jacobi equations
decouple, x'' = a x and y'' = -a y. The profile is given by its Taylor
coefficients at u = 0; the fundamental solutions from u = 0 follow from the
power-series recurrence in 60-digit decimal arithmetic (standard library only)
and are stored as Chebyshev series, so every two-point block is algebraic in
eight one-variable functions. No estimator imports.

Only weak profiles without conjugate points on [0, 1] are supported; positivity
of every Jacobi block is checked at use, not assumed. The floating guard is an
explicit error budget, not interval arithmetic. Margins inside the guard are
decided by evaluating the power series in decimal arithmetic and are counted.
"""

import math
from decimal import Decimal, localcontext

import numpy as np

DEGREE = 48
PRECISION = 60
ORDER = 144
# Assumed absolute bound on every stored fundamental value in float64. The
# selfcheck measures the actual error against the decimal power series.
BLOCK_ERROR = 1e-13


def _horner(series, u):
    total = Decimal(0)
    for coefficient in reversed(series):
        total = total * u + coefficient
    return total


def decimal_pi():
    """Machin's formula; no stored constant."""
    with localcontext() as context:
        context.prec = PRECISION + 10

        def arctan_inverse(n):
            power = total = 1 / Decimal(n)
            square, j = Decimal(n) * n, 0
            while abs(power) > Decimal(10) ** -(PRECISION + 8):
                power, j = -power / square, j + 1
                total += power / (2 * j + 1)
            return total

        return 16 * arctan_inverse(5) - 4 * arctan_inverse(239)


class ProfileWave:
    def __init__(self, taylor, bound, label):
        """taylor: Decimal Taylor coefficients of a(u) at 0; bound >= max |a|."""
        if not 0 <= bound <= 0.75:
            raise ValueError("Outside the bounded weak-profile generator family")
        self.bound, self.label = float(bound), str(label)
        self.resolved_margins = 0
        count = DEGREE + 1
        with localcontext() as context:
            context.prec = PRECISION
            taylor = [Decimal(value) for value in taylor]
            tiny = Decimal(10) ** -45
            # Axes: sign (+a, -a), function (c, c', s, s'), Taylor order.
            self._taylor = []
            for sign in (1, -1):
                functions = []
                for start in ((1, 0), (0, 1)):
                    x = [Decimal(start[0]), Decimal(start[1])]
                    for n in range(ORDER - 1):
                        total = Decimal(0)
                        for k in range(min(n, len(taylor) - 1) + 1):
                            total += taylor[k] * x[n - k]
                        x.append(sign * total / ((n + 2) * (n + 1)))
                    if max(abs(value) for value in x[-8:]) > tiny:
                        raise ArithmeticError("Taylor order too low for this profile")
                    functions.append(x)
                    functions.append([(n + 1) * x[n + 1] for n in range(ORDER)])
                self._taylor.append(functions)
            pi = decimal_pi()
            cosines = []
            for k in range(count):
                angle = pi * (2 * k + 1) / (2 * count)
                term = total = Decimal(1)
                j = 0
                while abs(term) > Decimal(10) ** -(PRECISION - 2):
                    j += 2
                    term = -term * angle * angle / (j * (j - 1))
                    total += term
                cosines.append(total)
            # T_j(cos theta_k) by the three-term recurrence.
            chebyshev = [[Decimal(1)] * count, cosines]
            for _ in range(2, count):
                chebyshev.append(
                    [
                        2 * c * a - b
                        for c, a, b in zip(
                            cosines, chebyshev[-1], chebyshev[-2], strict=True
                        )
                    ]
                )
            self.coefficients = np.empty((2, 4, count))
            for sign in range(2):
                for function in range(4):
                    series = self._taylor[sign][function]
                    values = [_horner(series, (c + 1) / 2) for c in cosines]
                    for j in range(count):
                        total = Decimal(0)
                        for value, basis in zip(values, chebyshev[j], strict=True):
                            total += value * basis
                        self.coefficients[sign, function, j] = float(
                            total * 2 / count / (2 if j == 0 else 1)
                        )

    def fundamental(self, u, xp=np):
        """Return c, c', s, s' for both signs; shape (2, 4) + u.shape."""
        u = xp.asarray(u, dtype=xp.float64)
        if not bool(xp.isfinite(u).all()) or bool(((u < 0) | (u > 1)).any()):
            raise ValueError("Finite Brinkmann u in [0, 1] required")
        coefficients = xp.asarray(self.coefficients)
        shape = (2, 4) + (1,) * u.ndim
        t = (2 * u - 1)[None, None]
        first = xp.zeros((2, 4) + u.shape)
        second = xp.zeros_like(first)
        for j in range(DEGREE, 0, -1):
            first, second = (
                coefficients[:, :, j].reshape(shape) + 2 * t * first - second,
                first,
            )
        return coefficients[:, :, 0].reshape(shape) + t * first - second

    def _pair(self, left, right, xp=np):
        """Action between two events, its guard, and the unusable-block mask.

        With M = M(right <- left) = [[A, B], [C, D]] for each sign, the action
        of the time-ordered pair is (A l^2 - 2 l r + D r^2) / (2 |B|) in either
        time order, because M(left <- right) is the symplectic inverse.
        """
        fl = self.fundamental(left[..., 0], xp)
        fr = self.fundamental(right[..., 0], xp)
        a = fr[:, 0] * fl[:, 3] - fr[:, 2] * fl[:, 1]
        b = fr[:, 2] * fl[:, 0] - fr[:, 0] * fl[:, 2]
        d = fr[:, 3] * fl[:, 0] - fr[:, 1] * fl[:, 2]
        lt = xp.moveaxis(left[..., 2:], -1, 0)
        rt = xp.moveaxis(right[..., 2:], -1, 0)
        numerator = a * lt * lt - 2 * lt * rt + d * rt * rt
        terms = abs(a) * lt * lt + 2 * abs(lt * rt) + abs(d) * rt * rt
        eps, block = np.finfo(float).eps, 8 * BLOCK_ERROR
        size = abs(b)
        unusable = (size <= 4 * block).any(axis=0)
        size = xp.where(size > 4 * block, size, 1.0)
        cost = (numerator / (2 * size)).sum(axis=0)
        error = (
            (block * (lt * lt + rt * rt) + block * abs(numerator) / size)
            / (1 - block / size)
            + 64 * eps * terms
        ) / (2 * size)
        return cost, 4 * error.sum(axis=0) * (1 + 2.0**-40), unusable

    def _resolve(self, left, right):
        """Decide one guarded margin from the power series in decimal arithmetic."""
        with localcontext() as context:
            context.prec = PRECISION
            u1, u2 = Decimal(float(left[0])), Decimal(float(right[0]))
            dv = Decimal(float(right[1])) - Decimal(float(left[1]))
            if u1 == u2:
                return False, False
            cost = Decimal(0)
            for sign in range(2):
                c1, cp1, s1, sp1 = (_horner(f, u1) for f in self._taylor[sign])
                c2, cp2, s2, sp2 = (_horner(f, u2) for f in self._taylor[sign])
                a, b, d = c2 * sp1 - s2 * cp1, s2 * c1 - c2 * s1, sp2 * c1 - cp2 * s1
                lt, rt = Decimal(float(left[2 + sign])), Decimal(float(right[2 + sign]))
                cost += (a * lt * lt - 2 * lt * rt + d * rt * rt) / (2 * abs(b))
            margin = (-dv if u2 > u1 else dv) - cost
            if abs(margin) < Decimal(10) ** -30:
                raise ArithmeticError("Margin unresolved even in decimal arithmetic")
            related = margin > 0
            return bool(related and u2 > u1), bool(related and u2 < u1)

    def relations(self, left, right, xp=np):
        """Broadcastable pair queries with the conventions of rotating_plane_wave."""
        left = xp.asarray(left, dtype=xp.float64)
        right = xp.asarray(right, dtype=xp.float64)
        du, dv = right[..., 0] - left[..., 0], right[..., 1] - left[..., 1]
        # The one-variable series are evaluated before broadcasting to pairs.
        cost, error, unusable = self._pair(left, right, xp)
        margin = -xp.sign(du) * dv - cost
        guard = error + 256 * np.finfo(float).eps * (1 + abs(dv))
        near = ((abs(margin) <= guard) | unusable) & (du != 0)
        related = (margin > guard) & (du != 0) & ~unusable
        future, past = related & (du > 0), related & (du < 0)
        if bool(near.any()):

            def host(value):
                return value.get() if hasattr(value, "get") else np.asarray(value)

            lefts, rights = xp.broadcast_arrays(left, right)
            if near.ndim == 0:
                self.resolved_margins += 1
                decided = self._resolve(host(lefts), host(rights))
                return xp.asarray(decided[0]), xp.asarray(decided[1])
            index = xp.nonzero(near)
            keys = np.column_stack([host(i) for i in index])
            lefts, rights = host(lefts[index]), host(rights[index])
            future, past = xp.array(future), xp.array(past)
            for key, one, other in zip(keys, lefts, rights, strict=True):
                key = tuple(int(i) for i in key)
                future[key], past[key] = self._resolve(one, other)
            self.resolved_margins += len(keys)
        return future, past

    def packed_relations(self, landmarks, events, xp=np):
        if len(events) % 8:
            raise ValueError("Measurement batches must be divisible by eight")
        future, past = self.relations(
            xp.asarray(landmarks, dtype=xp.float64)[:, None],
            xp.asarray(events, dtype=xp.float64)[None],
            xp,
        )
        packed = [
            xp.packbits(value).reshape(len(landmarks), -1) for value in (future, past)
        ]
        return tuple(p.get() if hasattr(p, "get") else p for p in packed)

    def order(self, points):
        indices = np.triu_indices(len(points), 1)
        result = np.zeros((len(points), len(points)), dtype=bool)
        for start in range(0, len(indices[0]), 262144):
            i, j = (k[start : start + 262144] for k in indices)
            future, past = self.relations(points[i], points[j])
            result[i, j], result[j, i] = future, past
        return result

    def slice_matrices(self, u):
        """Diagonal Hessians of the two tip actions at Brinkmann time u."""
        f, end = self.fundamental(u), self.fundamental(np.array(1.0))
        left_b = f[:, 2]
        right_b = end[:, 2, None] * f[:, 0] - end[:, 0, None] * f[:, 2]
        right_a = end[:, 0, None] * f[:, 3] - end[:, 2, None] * f[:, 1]
        if np.any(left_b <= 0) or np.any(right_b <= 0):
            raise ArithmeticError("Nonpositive Jacobi block; conjugate point or tip")
        left, right = f[:, 3] / left_b, right_a / right_b
        if np.any(left <= 0) or np.any(right <= 0):
            raise ArithmeticError("Nonpositive tip-action Hessian")
        return left, right, left + right

    def sample(self, n, rng):
        """Physical-volume rejection sample in the actual curved diamond.

        Mirrors RotatingWave.sample draw for draw; det g = -1 in this chart.
        """
        if n < 2:
            raise ValueError("At least two events required")
        ceiling = np.nextafter(np.pi / (4 * (1 - self.bound / 2)), np.inf)
        accepted, count = [], 0
        while count < n:
            u = rng.random(max(16, math.ceil(2.5 * (n - count))))
            u = u[u > 0]
            q = self.slice_matrices(u)[2]
            density = np.pi / np.sqrt(q[0] * q[1])
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
        z = np.sqrt(2 * radius2) * np.array([np.cos(angle), np.sin(angle)])
        y = z / np.sqrt(q)
        lo = -0.5 + 0.5 * (right * y * y).sum(axis=0)
        hi = 0.5 - 0.5 * (left * y * y).sum(axis=0)
        if np.any(hi <= lo):
            raise ArithmeticError("Nonpositive conditional sampling width")
        return np.column_stack((u, lo + uniform[:, 2] * (hi - lo), y.T))

    def volumes(self, points, quadrature_order=64):
        """Exact-geometry interval volumes by one-dimensional quadrature."""
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
                cost, _, unusable = self._pair(p, q)
                slack = p[1] - q[1] - cost
                if unusable or slack <= 0:
                    continue
                fp, fq = self.fundamental(p[0]), self.fundamental(q[0])
                fm = self.fundamental(p[0] + du * nodes)
                total = fq[:, 2] * fp[:, 0] - fq[:, 0] * fp[:, 2]
                first = fm[:, 2] * fp[:, 0, None] - fm[:, 0] * fp[:, 2, None]
                last = fq[:, 2, None] * fm[:, 0] - fq[:, 0, None] * fm[:, 2]
                if np.any(total <= 0) or np.any(first <= 0) or np.any(last <= 0):
                    raise ArithmeticError("Nonpositive Jacobi block in a volume")
                integrand = np.sqrt(np.prod(first * last / total[:, None], axis=0))
                output[i, j] = output[j, i] = (
                    np.pi * slack**2 * du * (weights @ integrand)
                )
        return output

    def normal_coordinates(self, points):
        """Truth-only inverse exponential map in the parallel midpoint frame.

        With the Jacobi solution based at u = 1/2, x(u) = s_m(u) x'(1/2), and
        the conserved Lagrangian gives the normal v-component in closed form:
        v_n = v - (U/2) |x'_0|^2 + (x x' + y y') / 2 at the endpoint.
        """
        u = points[:, 0] - 0.5
        if np.any(abs(u) < 1e-7):
            raise ArithmeticError("Normal map is ill-conditioned on the mid slice")
        f, mid = self.fundamental(points[:, 0]), self.fundamental(np.array(0.5))
        based = f[:, 2] * mid[:, 0, None] - f[:, 0] * mid[:, 2, None]
        slope = f[:, 3] * mid[:, 0, None] - f[:, 1] * mid[:, 2, None]
        transverse = points[:, 2:].T
        velocity = transverse / based
        v = (
            points[:, 1]
            - 0.5 * u * (velocity**2).sum(axis=0)
            + 0.5 * (transverse * slope * velocity).sum(axis=0)
        )
        x = velocity * u
        return np.column_stack(
            ((u - v) / 2, x[0] / np.sqrt(2), x[1] / np.sqrt(2), (u + v) / 2)
        )


def half_cycle(amplitude):
    """a(u) = amplitude * sin(pi u): one crest-centered half cycle in the diamond."""
    with localcontext() as context:
        context.prec = PRECISION
        amplitude, pi = Decimal(amplitude), decimal_pi()
        taylor, term = [], pi
        for k in range(ORDER):
            if k % 2:
                # Odd order k holds amplitude * (-1)^m pi^k / k!, k = 2m + 1.
                taylor.append(amplitude * term)
                term = -term * pi * pi / ((k + 1) * (k + 2))
            else:
                taylor.append(Decimal(0))
    return ProfileWave(taylor, abs(float(amplitude)), "half_cycle")


def constant(amplitude):
    return ProfileWave([Decimal(amplitude)], abs(float(amplitude)), "constant")


def selfcheck():
    """Deterministic engineering checks; no benchmark sample is drawn here."""
    from rotating_plane_wave import RotatingWave
    from scipy.integrate import solve_ivp

    report = {"decimal_pi_matches_float": float(decimal_pi()) == math.pi}
    assert report["decimal_pi_matches_float"]
    # 1. The constant profile must reproduce the audited closed-form generator.
    reference, wave = RotatingWave(0.25, 0.0), constant("0.25")
    ours = wave.sample(20000, np.random.default_rng(7))
    theirs = reference.sample(20000, np.random.default_rng(7))
    report["constant_sample_max_difference"] = float(np.max(abs(ours - theirs)))
    assert report["constant_sample_max_difference"] < 1e-12
    left, right = theirs[:10000], theirs[10000:]
    cost, guard, unusable = wave._pair(left, right)
    du = right[:, 0] - left[:, 0]
    first = np.where((du > 0)[:, None], left[:, 2:], right[:, 2:])
    last = np.where((du > 0)[:, None], right[:, 2:], left[:, 2:])
    expected, _ = reference.action(abs(du), first, last)
    report["constant_action_max_difference"] = float(np.max(abs(cost - expected)))
    assert not unusable.any() and np.all(abs(cost - expected) < guard)
    tips = np.array([[0, 0.5, 0, 0], [1, -0.5, 0, 0.0]])
    points = np.vstack((tips[:1], theirs[:400], tips[1:]))
    assert np.array_equal(wave.order(points), reference.order(points))
    landmarks = np.vstack((tips[:1], theirs[400:424], tips[1:]))
    report["constant_volume_max_difference"] = float(
        np.max(abs(wave.volumes(landmarks) - reference.volumes(landmarks)))
    )
    assert report["constant_volume_max_difference"] < 1e-12
    report["constant_normal_max_difference"] = float(
        np.max(
            abs(
                wave.normal_coordinates(landmarks)
                - reference.normal_coordinates(landmarks)
            )
        )
    )
    assert report["constant_normal_max_difference"] < 1e-11
    flat, flat_reference = constant("0"), RotatingWave(0.0, 0.0)
    assert np.array_equal(flat.order(points), flat_reference.order(points))
    report["flat_sample_max_difference"] = float(
        np.max(
            abs(
                flat.sample(4096, np.random.default_rng(3))
                - flat_reference.sample(4096, np.random.default_rng(3))
            )
        )
    )
    assert report["flat_sample_max_difference"] < 1e-12
    # 2. Stored series against the decimal power series, and the Wronskian.
    with localcontext() as context:
        context.prec = PRECISION
        amplitude = decimal_pi() ** 2 / 200
    wave, a0 = half_cycle(amplitude), float(amplitude)
    grid = np.random.default_rng(20261111).random(400)
    values = wave.fundamental(grid)
    worst = 0.0
    with localcontext() as context:
        context.prec = PRECISION
        for sign in range(2):
            for function in range(4):
                for k, u in enumerate(grid):
                    exact = _horner(wave._taylor[sign][function], Decimal(float(u)))
                    worst = max(worst, abs(float(exact) - values[sign, function, k]))
    report["half_cycle_series_max_error"] = worst
    assert worst < BLOCK_ERROR / 10
    wronskian = values[:, 0] * values[:, 3] - values[:, 1] * values[:, 2]
    report["half_cycle_wronskian_max_error"] = float(np.max(abs(wronskian - 1)))
    assert report["half_cycle_wronskian_max_error"] < 1e-14

    # 3. Independent ODE checks of the series, action, relation and normal map.
    def profile(u):
        return a0 * np.sin(np.pi * u), a0 * np.pi * np.cos(np.pi * u)

    ode_error = 0.0
    for sign_index, sign in enumerate((1, -1)):
        result = solve_ivp(
            lambda u, y, sign=sign: [
                y[1],
                sign * profile(u)[0] * y[0],
                y[3],
                sign * profile(u)[0] * y[2],
            ],
            [0, 1],
            [1, 0, 0, 1],
            rtol=1e-13,
            atol=1e-15,
            dense_output=True,
            method="DOP853",
        )
        assert result.success
        dense = result.sol(np.sort(grid))
        ode_error = max(
            ode_error,
            float(np.max(abs(dense - wave.fundamental(np.sort(grid))[sign_index]))),
        )
    report["half_cycle_independent_ode_max_difference"] = ode_error
    assert ode_error < 1e-9
    sample = wave.sample(64, np.random.default_rng(11))
    action_error = normal_error = 0.0
    for p, q in zip(sample[:32], sample[32:], strict=True):
        if p[0] > q[0]:
            p, q = q, p
        fp, fq = wave.fundamental(p[0]), wave.fundamental(q[0])
        b = fq[:, 2] * fp[:, 0] - fq[:, 0] * fp[:, 2]
        a = fq[:, 0] * fp[:, 3] - fq[:, 2] * fp[:, 1]
        initial_velocity = (q[2:] - a * p[2:]) / b

        def transverse_rhs(u, state):
            h = profile(u)[0]
            x, y, xp, yp = state[:4]
            lagrangian = 0.5 * (xp * xp + yp * yp + h * (x * x - y * y))
            return [xp, yp, h * x, -h * y, lagrangian]

        result = solve_ivp(
            transverse_rhs,
            [p[0], q[0]],
            [*p[2:], *initial_velocity, 0],
            rtol=1e-12,
            atol=1e-14,
        )
        assert result.success and np.max(abs(result.y[:2, -1] - q[2:])) < 1e-10
        cost, guard, unusable = wave._pair(p, q)
        action_error = max(action_error, abs(float(result.y[-1, -1]) - float(cost)))
        assert not unusable
        future, past = wave.relations(p, q)
        assert bool(future) == bool(p[1] - q[1] > result.y[-1, -1]) and not past
        back_future, back_past = wave.relations(q, p)
        assert not back_future and bool(back_past) == bool(future)
    report["half_cycle_independent_action_max_error"] = float(action_error)
    assert action_error < 1e-10

    def geodesic_rhs(_, state):
        u, _, x, y, du, dv, dx, dy = state
        h, dh = profile(u)
        vdd = -0.5 * du**2 * dh * (x * x - y * y) - 2 * du * h * (dx * x - dy * y)
        return [du, dv, dx, dy, 0, vdd, du**2 * h * x, -(du**2) * h * y]

    for p, normal in zip(sample, wave.normal_coordinates(sample), strict=True):
        du, dv = normal[0] + normal[3], normal[3] - normal[0]
        initial = [0.5, 0, 0, 0, du, dv, *(np.sqrt(2) * normal[1:3])]
        result = solve_ivp(geodesic_rhs, [0, 1], initial, rtol=1e-12, atol=1e-14)
        assert result.success
        normal_error = max(normal_error, float(np.max(abs(result.y[:4, -1] - p))))
    report["half_cycle_geodesic_roundtrip_max_error"] = normal_error
    assert normal_error < 1e-10
    # 4. Volume quadrature, the diamond volume and the rejection ceiling.
    landmarks = np.vstack((tips[:1], sample[:24], tips[1:]))
    fine, coarse = wave.volumes(landmarks), wave.volumes(landmarks, 32)
    report["half_cycle_volume_refinement"] = float(np.max(abs(fine - coarse)))
    assert report["half_cycle_volume_refinement"] < 1e-12
    nodes, weights = np.polynomial.legendre.leggauss(64)
    nodes, weights = (nodes + 1) / 2, weights / 2
    q = wave.slice_matrices(nodes)[2]
    density = np.pi / np.sqrt(q[0] * q[1])
    assert abs(weights @ density - fine[0, -1]) < 1e-12
    assert np.all(density <= np.pi * nodes * (1 - nodes) / (1 - a0 / 2) * (1 + 1e-13))
    report["half_cycle_diamond_volume_over_flat"] = float(fine[0, -1] / (np.pi / 6))
    # 5. A guarded margin is decided in decimal arithmetic, not guessed.
    p = sample[0]
    q = p.copy()
    q[0] = min(p[0] + 0.25, 0.999)
    cost, _, _ = wave._pair(p, q)
    q[1] = p[1] - float(cost)
    before = wave.resolved_margins
    wave.relations(p, q)
    report["guarded_margin_resolved_in_decimal_arithmetic"] = (
        wave.resolved_margins == before + 1
    )
    assert report["guarded_margin_resolved_in_decimal_arithmetic"]
    return report


if __name__ == "__main__":
    import json

    print(json.dumps(selfcheck(), indent=2))
