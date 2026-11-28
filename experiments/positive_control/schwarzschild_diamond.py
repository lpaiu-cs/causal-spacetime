"""Generator side of a Schwarzschild test geometry for the order-only Weyl inverse.

Isotropic coordinates (t, X): ds^2 = -A^2 dt^2 + B^4 dX.dX with A = (1-m/rho)/(1+m/rho),
B = 1 + m/rho, rho = |X| and m = M/2. The geometry is static, so q is causally after
p exactly when t_q - t_p is at least the optical distance of their spatial points:
the Fermat minimum of the integral of n = B^2/A along paths. That minimum is found
here by Newton iteration on a low-order polynomial deviation from the chord; the
error of the time is second order in the error of the path.

The interval is the Alexandrov interval of two events on the radial geodesic that
is momentarily at rest at the apex, a proper time T/2 before and after it. Points
are (t, x, y, z); the apex is (0, rho_a, 0, 0). No estimator module imports this.
"""

import numpy as np
from numpy.polynomial import chebyshev, legendre
from scipy.integrate import solve_ivp

MODES, NODES = 10, 32
DIRECT_GUARD = 1e-13
SURROGATE_DEGREE = 24
BOUNDARY_LAYER = 1e-10


def _basis(modes=MODES, nodes=NODES):
    """Deviation modes s(1-s)P_k(2s-1) and their derivatives at Gauss nodes."""
    x, w = legendre.leggauss(nodes)
    s, w = (x + 1) / 2, w / 2
    p = np.stack([legendre.legval(x, np.eye(modes)[k]) for k in range(modes)], axis=1)
    dp = np.stack(
        [legendre.legval(x, legendre.legder(np.eye(modes)[k])) for k in range(modes)],
        axis=1,
    )
    phi = (s * (1 - s))[:, None] * p
    dphi = (1 - 2 * s)[:, None] * p + (2 * s * (1 - s))[:, None] * dp
    return s, w, phi, dphi


class SchwarzschildDiamond:
    def __init__(self, mass, apex_areal_radius, proper_time):
        if not (mass >= 0 and proper_time > 0 and apex_areal_radius > 4 * max(mass, 0)):
            raise ValueError("Outside the supported exterior domain")
        self.mass, self.m = float(mass), float(mass) / 2
        self.proper_time = float(proper_time)
        r = float(apex_areal_radius)
        self.apex_areal_radius = r
        self.apex = (r - mass + np.sqrt(r * r - 2 * mass * r)) / 2
        self.resolved_margins = 0
        self._modes = _basis()
        self._cache = {}
        half = self._tip()
        self.half_time, self.base = float(half[0]), np.array([half[1], 0.0, 0.0])
        self.tips = np.array(
            [[-self.half_time, *self.base], [self.half_time, *self.base]]
        )
        # n >= 1, so the optical ball of the tips lies in this coordinate ball.
        self.radius = self.half_time / self.index(self.base[0] + self.half_time)[0]

    # -- metric functions -------------------------------------------------
    def index(self, rho, xp=np):
        """n, dn/drho and d2n/drho2 of the Fermat index n = B^2/A."""
        u = self.m / rho
        n = (1 + u) ** 3 / (1 - u)
        first = -(u / rho) * (3 / (1 + u) + 1 / (1 - u))
        second = (2 * u / rho**2) * (3 / (1 + u) + 1 / (1 - u)) - (u / rho) ** 2 * (
            3 / (1 + u) ** 2 - 1 / (1 - u) ** 2
        )
        return n, n * first, n * (first**2 + second)

    def lapse_volume(self, rho):
        """A and sqrt(-g) = A B^6 in these coordinates."""
        u = self.m / rho
        lapse = (1 - u) / (1 + u)
        return lapse, lapse * (1 + u) ** 6

    # -- optical distance -------------------------------------------------
    def travel_time(self, left, right, modes=None, invariant=False):
        """Minimal coordinate time of light between spatial points, broadcast.

        With invariant=True also return the spread along the path of Bouguer's
        n rho sin(angle to the radius), which is constant on a true ray.
        """
        s, w, phi, dphi = self._modes if modes is None else modes
        left, right = np.broadcast_arrays(
            np.asarray(left, dtype=float), np.asarray(right, dtype=float)
        )
        shape = left.shape[:-1]
        left, right = left.reshape(-1, 3), right.reshape(-1, 3)
        chord = right - left
        d = np.linalg.norm(chord, axis=1)
        same = d == 0
        d = np.where(same, 1.0, d)
        e = chord / d[:, None]
        along = -np.einsum("pi,pi->p", left, e)  # distance from left to the foot
        h = np.linalg.norm(left + along[:, None] * e, axis=1)
        xi = d[:, None] * s[None] - along[:, None]
        axis = xi * xi
        a = np.zeros((len(d), phi.shape[1]))
        pp = (phi[:, :, None] * phi[:, None, :]).reshape(len(s), -1)
        pd = (
            phi[:, :, None] * dphi[:, None, :] + dphi[:, :, None] * phi[:, None, :]
        ).reshape(len(s), -1)
        dd = (dphi[:, :, None] * dphi[:, None, :]).reshape(len(s), -1)
        k = phi.shape[1]
        for _ in range(8):
            y, dy = h[:, None] + a @ phi.T, a @ dphi.T
            rho = np.sqrt(axis + y * y)
            speed = np.sqrt(d[:, None] ** 2 + dy * dy)
            n, n1, n2 = self.index(rho)
            f_y = n1 * y / rho * speed
            f_d = n * dy / speed
            gradient = (w * f_y) @ phi + (w * f_d) @ dphi
            f_yy = (n2 * (y / rho) ** 2 + n1 * axis / rho**3) * speed
            f_yd = n1 * y / rho * dy / speed
            f_dd = n * d[:, None] ** 2 / speed**3
            hessian = ((w * f_yy) @ pp + (w * f_yd) @ pd + (w * f_dd) @ dd).reshape(
                -1, k, k
            )
            a = a - np.linalg.solve(hessian, gradient[:, :, None])[:, :, 0]
        y, dy = h[:, None] + a @ phi.T, a @ dphi.T
        speed = np.sqrt(d[:, None] ** 2 + dy * dy)
        n = self.index(np.sqrt(axis + y * y))[0]
        time = np.where(same, 0.0, (n * speed) @ w).reshape(shape)
        if not invariant:
            return time
        bouguer = n * (xi * dy - y * d[:, None]) / speed
        return time, np.ptp(bouguer, axis=1).reshape(shape)

    # -- tips and frames --------------------------------------------------
    def _geodesic(self, _, state):
        """Affinely parametrized geodesic in (t, X) with velocity (t', X')."""
        x, v = state[1:4], state[5:8]
        rho = np.linalg.norm(x)
        u = self.m / rho
        lapse = (1 - u) / (1 + u)
        d_lapse = 2 * u / (rho * (1 + u) ** 2)  # dA/drho
        d_log_b = -u / (rho * (1 + u))  # d ln B / drho
        radial = x / rho
        grad_phi = 2 * d_log_b * radial
        accel_t = -2 * (d_lapse / lapse) * (radial @ v) * state[4]
        accel_x = (
            -(lapse * d_lapse / (1 + u) ** 4) * radial * state[4] ** 2
            - 2 * v * (v @ grad_phi)
            + (v @ v) * grad_phi
        )
        return np.r_[state[4:], accel_t, accel_x]

    def _shoot(self, velocity):
        start = np.r_[0.0, self.apex, 0.0, 0.0, velocity]
        return solve_ivp(
            self._geodesic, (0, 1), start, method="DOP853", rtol=1e-12, atol=1e-14
        ).y[:4, -1]

    def _tip(self):
        """Coordinate time and radius reached a proper time T/2 after the apex."""
        lapse = self.lapse_volume(self.apex)[0]
        end = self._shoot(np.r_[self.proper_time / 2 / lapse, 0.0, 0.0, 0.0])
        return end[0], end[1]

    def normal_coordinates(self, points):
        """Normal coordinates at the apex in its static frame, in units of T."""
        u = self.m / self.apex
        scale = np.r_[(1 - u) / (1 + u), np.full(3, (1 + u) ** 2)]
        result = []
        for target in np.asarray(points, dtype=float):
            velocity = target - np.r_[0.0, self.apex, 0.0, 0.0]
            for _ in range(8):
                miss = self._shoot(velocity) - target
                if np.max(abs(miss)) < 1e-11:
                    break
                jacobian = np.column_stack(
                    [
                        (self._shoot(velocity + step) - self._shoot(velocity - step))
                        / 2e-6
                        for step in 1e-6 * np.eye(4)
                    ]
                )
                velocity = velocity - np.linalg.solve(jacobian, miss)
            else:
                raise ArithmeticError("Exponential map did not converge")
            result.append(velocity * scale / self.proper_time)
        return np.array(result)

    def tidal(self):
        """T^2 E_ij at the apex in (radial, y, z) axes; the magnetic part is zero."""
        strength = self.proper_time**2 * self.mass / self.apex_areal_radius**3
        return strength * np.diag([-2.0, 1.0, 1.0])

    # -- surrogates of the squared optical distance -------------------------
    def _surrogates(self, points):
        """Per-point Chebyshev fit in (a, b) = (axial, squared transverse) position.

        The squared distance from a point is smooth there; the distance is not.
        """
        points = np.ascontiguousarray(points, dtype=float)
        key = points.tobytes()
        if key not in self._cache:
            degree = SURROGATE_DEGREE
            nodes = np.cos(np.pi * (np.arange(degree + 1) + 0.5) / (degree + 1))
            inverse = np.linalg.inv(chebyshev.chebvander(nodes, degree))
            rng = np.random.default_rng(20261112)
            direction = rng.normal(size=(512, 3))
            direction /= np.linalg.norm(direction, axis=1)[:, None]
            test = self.base + (
                self.radius * rng.random((512, 1)) ** (1 / 3) * direction
            )
            axes, bounds, coefficients, guards = [], [], [], []
            for point in points:
                axis = point / np.linalg.norm(point)
                offset = self.base - (self.base @ axis) * axis
                reach = np.linalg.norm(offset) + self.radius
                side = (
                    offset if reach > self.radius else np.cross(axis, [0.0, 0.0, 1.0])
                )
                if np.linalg.norm(side) < 1e-9:
                    side = np.cross(axis, [0.0, 1.0, 0.0])
                side = side / np.linalg.norm(side)
                centre = self.base @ axis
                bound = np.array(
                    [centre - self.radius, centre + self.radius, 0.0, reach**2]
                )
                a = (bound[0] + bound[1]) / 2 + (bound[1] - bound[0]) / 2 * nodes
                b = (bound[2] + bound[3]) / 2 + (bound[3] - bound[2]) / 2 * nodes
                grid = (
                    a[:, None, None] * axis + np.sqrt(b)[None, :, None] * side
                ).reshape(-1, 3)
                values = self.travel_time(point, grid).reshape(len(a), len(b)) ** 2
                fit = inverse @ values @ inverse.T
                error = abs(
                    self._evaluate(axis[None], bound[None], fit[None], test)[0]
                    - self.travel_time(point, test) ** 2
                )
                axes.append(axis)
                bounds.append(bound)
                coefficients.append(fit)
                guards.append(8 * error.max() + 1e-12)
            self._cache[key] = tuple(
                map(np.array, (axes, bounds, coefficients, guards))
            )
        return self._cache[key]

    @staticmethod
    def _evaluate(axes, bounds, coefficients, positions, xp=np):
        """Squared optical distance from every fitted point to every position."""
        a = xp.asarray(axes) @ positions.T
        b = (positions * positions).sum(axis=1)[None] - a * a
        bounds, coefficients = xp.asarray(bounds), xp.asarray(coefficients)
        terms = []
        for value, lo, hi in (
            (a, bounds[:, 0:1], bounds[:, 1:2]),
            (b, bounds[:, 2:3], bounds[:, 3:4]),
        ):
            x = (2 * value - (hi + lo)) / (hi - lo)
            series = [xp.ones_like(x), x]
            for _ in range(coefficients.shape[1] - 2):
                series.append(2 * x * series[-1] - series[-2])
            terms.append(xp.stack(series, axis=2))
        return (xp.matmul(terms[0], coefficients) * terms[1]).sum(axis=2)

    # -- relations ----------------------------------------------------------
    def relations(self, left, right):
        """Direct pair queries; an unresolved margin stops the generator."""
        left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
        dt = right[..., 0] - left[..., 0]
        margin = abs(dt) - self.travel_time(left[..., 1:], right[..., 1:])
        distinct = (dt != 0) | np.any(left[..., 1:] != right[..., 1:], axis=-1)
        if np.any((abs(margin) <= DIRECT_GUARD) & distinct):
            raise ArithmeticError("Unresolved Schwarzschild causal margin")
        related = (margin > DIRECT_GUARD) & distinct
        return related & (dt > 0), related & (dt < 0)

    def order(self, points):
        """Strict causal order matrix: entry (i, j) says i precedes j."""
        points = np.asarray(points, dtype=float)
        upper = np.column_stack(np.triu_indices(len(points), 1))
        result = np.zeros((len(points), len(points)), dtype=bool)
        for start in range(0, len(upper), 200_000):
            i, j = upper[start : start + 200_000].T
            future, past = self.relations(points[i], points[j])
            result[i[future], j[future]] = True
            result[j[past], i[past]] = True
        return result

    def packed_relations(self, landmarks, events, xp=np):
        if len(events) % 8:
            raise ValueError("Measurement batches must be divisible by eight")
        landmarks, events = np.asarray(landmarks, float), np.asarray(events, float)
        axes, bounds, coefficients, guards = self._surrogates(landmarks[:, 1:])
        dt = xp.asarray(events[:, 0])[None] - xp.asarray(landmarks[:, 0])[:, None]
        squared = self._evaluate(
            axes, bounds, coefficients, xp.asarray(events[:, 1:]), xp
        )
        margin = dt * dt - squared
        guard = xp.asarray(guards)[:, None] + 64 * np.finfo(float).eps * (
            dt * dt + abs(squared)
        )
        related, near = margin > guard, abs(margin) <= guard
        future, past = related & (dt > 0), related & (dt < 0)
        if bool(near.any()):
            host = (lambda v: v.get()) if hasattr(near, "get") else np.asarray
            rows, columns = (host(v) for v in xp.nonzero(near))
            decided = self.relations(landmarks[rows], events[columns])
            future, past = xp.array(future), xp.array(past)
            future[rows, columns] = xp.asarray(decided[0])
            past[rows, columns] = xp.asarray(decided[1])
            self.resolved_margins += len(rows)
        packed = [xp.packbits(v).reshape(len(landmarks), -1) for v in (future, past)]
        return tuple(p.get() if hasattr(p, "get") else p for p in packed)

    # -- sampling -----------------------------------------------------------
    def sample(self, n, rng):
        """Rejection sample, uniform in proper volume, of the tips' interval.

        Events closer than BOUNDARY_LAYER (in coordinate time) to the null
        boundary of the tips are not drawn, so every event is decidably inside.
        """
        if n < 2:
            raise ValueError("At least two events required")
        fit = self._surrogates(self.base[None])[:3]
        ceiling = self.lapse_volume(self.base[0] - self.radius)[1]
        accepted, count = [], 0
        while count < n:
            k = max(64, int(4 * (n - count)))
            radius = self.radius * rng.beta(3, 2, k)
            direction = rng.normal(size=(k, 3))
            direction /= np.linalg.norm(direction, axis=1)[:, None]
            position = self.base + radius[:, None] * direction
            distance = np.sqrt(np.maximum(self._evaluate(*fit, position)[0], 0))
            height = self.half_time - distance - BOUNDARY_LAYER
            density = self.lapse_volume(np.linalg.norm(position, axis=1))[1]
            ratio = (density / ceiling) * (
                height / (self.half_time * (1 - radius / self.radius))
            )
            if np.max(ratio) > 1 + 1e-9:
                raise ArithmeticError("Rejection envelope violated")
            keep = rng.random(k) < ratio
            time = (2 * rng.random(k) - 1) * height
            accepted.append(np.column_stack((time, position))[keep])
            count += int(keep.sum())
        return np.concatenate(accepted)[:n]

    def interval_nodes(self, nodes=48):
        """Gauss nodes on the spatial section of the interval, direct solver only.

        Returns offsets from the tips' position, coordinate-volume weights and the
        half-height in time. The proper volume of the interval is the sum of
        weight * sqrt(-g) * 2 * half-height.
        """
        mu, w_mu = legendre.leggauss(nodes)
        omega = np.column_stack((mu, np.sqrt(1 - mu * mu), np.zeros_like(mu)))
        edge, outside = np.zeros_like(mu), np.full_like(mu, self.radius)
        for _ in range(70):
            middle = (edge + outside) / 2
            near = (
                self.travel_time(self.base, self.base + middle[:, None] * omega)
                < self.half_time
            )
            edge, outside = (
                np.where(near, middle, edge),
                np.where(near, outside, middle),
            )
        x, w_r = legendre.leggauss(nodes)
        r = edge[:, None] * (x + 1) / 2
        offset = r[:, :, None] * omega[:, None, :]
        height = self.half_time - self.travel_time(self.base, self.base + offset)
        return offset, np.pi * (w_mu * edge)[:, None] * w_r * r * r, height


def selfcheck():
    """Deterministic checks of the solver, the frames and the generator path."""
    from scipy.optimize import brentq

    mass, apex, proper = 1.0, 6.0, 3.0
    g = SchwarzschildDiamond(mass, apex, proper)
    rng = np.random.default_rng(20261112)

    def ball(k):
        v = rng.normal(size=(k, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        return g.base + g.radius * rng.random((k, 1)) ** (1 / 3) * v

    rho = np.array([3.6, 4.9, 6.3])
    n, n1, n2 = g.index(rho)
    step = 1e-30j
    index_error = max(
        np.max(abs(g.index(rho + step)[0].imag / 1e-30 - n1)),
        np.max(abs(g.index(rho + step)[1].imag / 1e-30 - n2)),
    )
    assert index_error < 1e-13

    # Radial pairs: exact tortoise-coordinate difference in the areal radius.
    inner, outer = np.array([3.7, 4.2, 5.5]), np.array([6.2, 4.25, 3.6])
    areal = lambda x: x * (1 + g.m / x) ** 2  # noqa: E731
    tortoise = lambda r: r + 2 * mass * np.log(r / (2 * mass) - 1)  # noqa: E731
    zeros = np.zeros(3)
    radial = g.travel_time(
        np.column_stack((inner, zeros, zeros)), np.column_stack((outer, zeros, zeros))
    )
    radial_error = np.max(
        abs(radial - abs(tortoise(areal(outer)) - tortoise(areal(inner))))
    )
    assert radial_error < 1e-13

    a, b = ball(2000), ball(2000)
    time, bouguer = g.travel_time(a, b, invariant=True)
    bouguer = float(np.max(bouguer))
    refinement = np.max(abs(time - g.travel_time(a, b, modes=_basis(16, 64))))
    symmetry = np.max(abs(time - g.travel_time(b, a)))
    c = ball(2000)
    triangle = np.min(g.travel_time(a, c) + g.travel_time(c, b) - time)
    assert refinement < 1e-14 and symmetry < 1e-14 and triangle > -1e-12
    assert bouguer < 1e-7

    def ray_time(x1, x2):
        """Independent: integrate the ray equation and shoot for the endpoint."""
        e1 = x1 / np.linalg.norm(x1)
        e2 = x2 - (x2 @ e1) * e1
        e2 /= np.linalg.norm(e2)
        p1, p2 = np.array([x1 @ e1, 0.0]), np.array([x2 @ e1, x2 @ e2])
        d = np.linalg.norm(p2 - p1)
        along = (p2 - p1) / d
        across = np.array([-along[1], along[0]])

        def rhs(_, s):
            r = np.linalg.norm(s[:2])
            m0, m1, _ = g.index(np.array([r]))
            return np.r_[s[2:4] / m0[0], m1[0] * s[:2] / r, m0[0]]

        def arrive(_, s):
            return (s[:2] - p1) @ along - d

        arrive.terminal, arrive.direction = True, 1

        def end(angle):
            k0 = g.index(np.array([p1[0]]))[0][0] * np.array(
                [np.cos(angle), np.sin(angle)]
            )
            return solve_ivp(
                rhs,
                (0, 3 * d),
                np.r_[p1, k0, 0.0],
                method="DOP853",
                rtol=1e-13,
                atol=1e-15,
                events=arrive,
            ).y_events[0][0]

        heading = np.arctan2(*along[::-1])
        angle = brentq(
            lambda x: (end(x)[:2] - p1) @ across,
            heading - 0.3,
            heading + 0.3,
            xtol=1e-15,
        )
        return end(angle)[4]

    ray_error = max(abs(ray_time(a[i], b[i]) - time[i]) for i in range(12))
    assert ray_error < 1e-10

    # Geodesic integrator: norm, energy and angular momentum are conserved on the
    # radial tip geodesic and on generic timelike, spacelike and null geodesics.
    def invariants(state):
        u = g.m / np.linalg.norm(state[1:4])
        lapse, conformal = (1 - u) / (1 + u), (1 + u) ** 4
        return np.r_[
            -((lapse * state[4]) ** 2) + conformal * state[5:8] @ state[5:8],
            lapse**2 * state[4],
            conformal * np.cross(state[1:4], state[5:8]),
        ]

    lapse = g.lapse_volume(g.apex)[0]
    geodesic_error = 0.0
    for velocity in (
        [proper / 2 / lapse, 0.0, 0.0, 0.0],
        [1.9, 0.4, -0.7, 0.5],
        [0.6, -0.5, 0.8, -0.6],
        [g.index(g.apex)[0] * np.sqrt(0.77), 0.5, 0.6, -0.4],
    ):
        start = np.r_[0.0, g.apex, 0.0, 0.0, velocity]
        final = solve_ivp(
            g._geodesic, (0, 1), start, method="DOP853", rtol=1e-12, atol=1e-14
        ).y[:, -1]
        geodesic_error = max(
            geodesic_error, np.max(abs(invariants(final) - invariants(start)))
        )
    assert abs(invariants(start)[0]) < 1e-14 and geodesic_error < 1e-9

    # Tidal tensor from the relative motion of particles released at rest. The
    # separation is measured in the rest frame of the first particle: a bare
    # coordinate difference would add the squared acceleration of the static
    # observers to the radial component.
    tau, shift, measured = 0.05, 1e-4, []
    scale = (1 + g.m / g.apex) ** 2
    for axis in range(3):
        ends = []
        for sign in (0.0, 1.0):
            x0 = np.r_[g.apex, 0.0, 0.0] + sign * shift / scale * np.eye(3)[axis]
            rate = tau / g.lapse_volume(np.linalg.norm(x0))[0]
            ends.append(
                solve_ivp(
                    g._geodesic,
                    (0, 1),
                    np.r_[0.0, x0, rate, 0.0, 0.0, 0.0],
                    method="DOP853",
                    rtol=1e-13,
                    atol=1e-16,
                ).y[:, -1]
            )
        u = g.m / np.linalg.norm(ends[0][1:4])
        metric = np.r_[-(((1 - u) / (1 + u)) ** 2), np.full(3, (1 + u) ** 4)]
        gap, velocity = ends[1][:4] - ends[0][:4], ends[0][4:] / tau
        separation = np.sqrt((metric * gap) @ gap + ((metric * gap) @ velocity) ** 2)
        measured.append(-2 * (separation / shift - 1) / tau**2 * proper**2)
    tidal_error = np.max(abs(np.array(measured) - np.diag(g.tidal()))) / abs(
        g.tidal()[0, 0]
    )
    assert tidal_error < 2e-3

    tips = g.normal_coordinates(g.tips)
    tip_error = np.max(abs(tips - [[-0.5, 0, 0, 0], [0.5, 0, 0, 0]]))
    assert tip_error < 1e-9

    stream = np.random.default_rng(20261113)
    events = np.vstack([g.sample(4096, stream) for _ in range(2)])
    inside = (
        g.relations(g.tips[0][None], events)[0]
        & g.relations(events, g.tips[1][None])[0]
    )
    landmarks = np.vstack([g.tips[0], events[:32], g.tips[1]])
    packed = g.packed_relations(landmarks, events)
    direct = g.relations(landmarks[:, None], events[None])
    bits = all(
        np.array_equal(np.unpackbits(p, axis=1).astype(bool), d)
        for p, d in zip(packed, direct, strict=True)
    )
    order = g.order(np.vstack([g.tips[0], events[:200], g.tips[1]]))
    closure = (order.astype(np.int32) @ order.astype(np.int32)) > 0
    assert inside.all() and bits and order[0, 1:].all() and order[:-1, -1].all()
    assert not (order & order.T).any() and not (closure & ~order).any()

    # Sampler: moments under the proper-volume measure of the interval, by a
    # quadrature that uses the direct solver only, against a large sample.
    node, measure, height = g.interval_nodes()
    isotropic = np.linalg.norm(g.base + node, axis=2)
    weight = measure * g.lapse_volume(isotropic)[1] * height
    expected = [
        float((weight * value).sum() / weight.sum())
        for value in (
            node[..., 0],
            (node * node).sum(axis=2),
            height * height / 3,
            isotropic,
        )
    ]
    drawn = g.sample(2_000_000, np.random.default_rng(20261115))
    offset = drawn[:, 1:] - g.base
    statistics = (
        offset[:, 0],
        (offset * offset).sum(axis=1),
        drawn[:, 0] ** 2,
        np.linalg.norm(drawn[:, 1:], axis=1),
    )
    sampler_z = [
        float((value.mean() - target) / (value.std() / np.sqrt(len(value))))
        for value, target in zip(statistics, expected, strict=True)
    ]
    assert max(map(abs, sampler_z)) < 4.5

    # Volume: a small interval exceeds the flat volume pi T^4 / 24 by the
    # curvature-square law 17 T^4 E.E / 10080 of a purely electric field. This
    # ties metric, solver, volume element and tip proper time together.
    def volume(geometry):
        node, measure, height = geometry.interval_nodes()
        radius = np.linalg.norm(geometry.base + node, axis=2)
        return float(2 * (measure * geometry.lapse_volume(radius)[1] * height).sum())

    small = SchwarzschildDiamond(mass, apex, 0.5)
    volume_law = (volume(small) / (np.pi * 0.5**4 / 24) - 1) / (
        17 * np.sum(small.tidal() ** 2) / 10080
    )
    assert abs(volume_law - 1) < 1e-3

    flat = SchwarzschildDiamond(0.0, apex, proper)
    flat_volume = abs(volume(flat) / (np.pi * proper**4 / 24) - 1)
    assert flat_volume < 1e-12
    flat_events = flat.sample(4096, np.random.default_rng(20261114))
    flat_error = max(
        np.max(abs(flat.travel_time(a, b) - np.linalg.norm(a - b, axis=1))),
        abs(flat.half_time - proper / 2),
        np.max(
            abs(
                flat.normal_coordinates(flat_events[:4])
                - (flat_events[:4] - [0, apex, 0, 0]) / proper
            )
        ),
    )
    reach = abs(flat_events[:, 0]) + np.linalg.norm(
        flat_events[:, 1:] - flat.base, axis=1
    )
    assert flat_error < 1e-13 and reach.max() < proper / 2
    return {
        "index_derivative_error": float(index_error),
        "radial_closed_form_error": float(radial_error),
        "mode_refinement_error": float(refinement),
        "bouguer_invariant_spread": bouguer,
        "symmetry_error": float(symmetry),
        "triangle_inequality_margin": float(triangle),
        "independent_ray_integration_error": float(ray_error),
        "geodesic_first_integral_error": float(geodesic_error),
        "tidal_tensor_relative_error_from_released_particles": float(tidal_error),
        "tip_normal_coordinate_error": float(tip_error),
        "sampled_events_inside_tips_interval": bool(inside.all()),
        "sampler_moment_z_scores": sampler_z,
        "interval_volume_over_flat": volume(g) / (np.pi * proper**4 / 24),
        "small_interval_volume_excess_over_curvature_square_law": float(volume_law),
        "flat_interval_volume_error": float(flat_volume),
        "surrogate_bits_equal_direct_solver": bool(bits),
        "surrogate_margins_resolved_directly": int(g.resolved_margins),
        "order_antisymmetric_and_transitive": True,
        "flat_limit_error": float(flat_error),
        "tidal_norm": float(np.sqrt(np.sum(g.tidal() ** 2))),
    }


if __name__ == "__main__":
    import json

    print(json.dumps(selfcheck(), indent=2))
