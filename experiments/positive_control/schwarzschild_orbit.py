"""Generator side of a Schwarzschild interval whose tips move through the apex.

This extends schwarzschild_diamond.SchwarzschildDiamond, whose isotropic metric,
optical-distance solver, geodesic integrator, Chebyshev surrogates and guarded
relation decisions are reused unchanged. Everything that class assumed about
radial tips at rest is replaced:

- Tips. The geodesic through the apex (t, X) = (0, rho_a, 0, 0) with
  four-velocity u = gamma (e_t + v e_y) in the static orthonormal frame, at
  proper times -/+ T/2. With v = sqrt(M / (r_a - 2M)) it is the circular
  geodesic of areal radius r_a; with v = 0 it is the radial geodesic at rest at
  the apex and the old geometry is reproduced. The past tip is the image of the
  future one under (t, x, y, z) -> (-t, x, -y, z), an isometry that fixes the
  apex and reverses u.
- Region. With foci P and Q, the tips' spatial points, the spatial section is
  the optical ellipsoid d(P, X) + d(X, Q) <= t_Q - t_P. It lies in the
  coordinate ball about c = (P + Q) / 2 of radius t_h / n(|c| + t_h), with
  t_h = (t_Q - t_P) / 2. That ball is the parent's `base` and `radius`, so it
  bounds the surrogates and the sampler.
- Sampling. Spatial density A B^6 h(X), h = t_Q - d(X, Q) - t_P - d(P, X);
  then t uniform on [t_P + d(P, X), t_Q - d(X, Q)].
- Quadrature. Moments of the measure split into one term per focus, each
  smooth along rays from its own focus.
- Normal coordinates and truth. The comoving orthonormal frame
  (u, e_x, gamma (v e_t + e_y), e_z), i.e. the static frame boosted along y.

Points are (t, x, y, z). No estimator module imports this.
"""

import numpy as np
from numpy.polynomial import legendre
from schwarzschild_diamond import BOUNDARY_LAYER, SchwarzschildDiamond, _basis
from scipy.integrate import solve_ivp

EPSILON = np.zeros((3, 3, 3))
for _i, _j, _k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    EPSILON[_i, _j, _k], EPSILON[_i, _k, _j] = 1.0, -1.0


def circular_velocity(mass, areal_radius):
    """Speed of the circular geodesic relative to the static observer."""
    return float(np.sqrt(mass / (areal_radius - 2 * mass)))


def static_riemann(mass, areal_radius):
    """R_abcd of Schwarzschild in the static orthonormal frame (t, r, y, z).

    MTW sign convention, so that R_0i0j is the tidal tensor E_ij:
    R_0r0r = -2 M/r^3, R_0y0y = R_0z0z = M/r^3, R_yzyz = 2 M/r^3,
    R_ryry = R_rzrz = -M/r^3.
    """
    psi = mass / areal_radius**3
    tensor = np.zeros((4, 4, 4, 4))
    for (a, b, c, d), value in (
        ((0, 1, 0, 1), -2 * psi),
        ((0, 2, 0, 2), psi),
        ((0, 3, 0, 3), psi),
        ((2, 3, 2, 3), 2 * psi),
        ((1, 2, 1, 2), -psi),
        ((1, 3, 1, 3), -psi),
    ):
        for p, q, s in (((a, b), (c, d), 1), ((c, d), (a, b), 1)):
            tensor[p[0], p[1], q[0], q[1]] = s * value
            tensor[p[1], p[0], q[0], q[1]] = -s * value
            tensor[p[0], p[1], q[1], q[0]] = -s * value
            tensor[p[1], p[0], q[1], q[0]] = s * value
    return tensor


def boost(velocity):
    """Rows: the comoving frame vectors in the static orthonormal frame."""
    gamma = 1 / np.sqrt(1 - velocity * velocity)
    frame = np.eye(4)
    frame[0, [0, 2]] = gamma, gamma * velocity
    frame[2, [0, 2]] = gamma * velocity, gamma
    return frame


def electric_magnetic(tensor):
    """E_ij = C_0i0j and B_ij = eps_ikl C_kl0j / 2 of a four-index tensor."""
    electric = tensor[0, 1:, 0, 1:]
    magnetic = 0.5 * np.einsum("ikl,klj->ij", EPSILON, tensor[1:, 1:, 0, 1:])
    return electric, magnetic


class MovingSchwarzschildDiamond(SchwarzschildDiamond):
    def __init__(self, mass, apex_areal_radius, proper_time, velocity):
        v = float(velocity)
        if not (
            mass >= 0
            and proper_time > 0
            and apex_areal_radius > 4 * max(mass, 0)
            and 0 <= v < 1
        ):
            raise ValueError("Outside the supported exterior domain")
        self.mass, self.m = float(mass), float(mass) / 2
        self.proper_time = float(proper_time)
        r = float(apex_areal_radius)
        self.apex_areal_radius = r
        self.apex = (r - mass + np.sqrt(r * r - 2 * mass * r)) / 2
        self.velocity, self.gamma = v, 1 / np.sqrt(1 - v * v)
        self.resolved_margins = 0
        self._modes = _basis()
        self._cache = {}
        lapse = self.lapse_volume(self.apex)[0]
        conformal = (1 + self.m / self.apex) ** 2
        half = self.proper_time / 2 * self.gamma
        # Coordinate components of (T/2) u at the apex; at v = 0 these are the
        # old class's operands bit for bit.
        self.tip_velocity = np.r_[half / lapse, 0.0, half * v / conformal, 0.0]
        future = self._shoot(self.tip_velocity)
        past = future * np.r_[-1.0, 1.0, -1.0, 1.0]
        self.tips = np.array([past, future])
        self.half_time = float((future[0] - past[0]) / 2)
        self.foci = self.tips[:, 1:].copy()
        # Centre and radius of the coordinate ball that holds the interval.
        self.base = self.foci.mean(axis=0)
        self.radius = (
            self.half_time / self.index(np.linalg.norm(self.base) + self.half_time)[0]
        )

    # -- frames and truth -----------------------------------------------------
    def tetrad(self):
        """Rows: coordinate components of the comoving frame at the apex."""
        u = self.m / self.apex
        lapse, conformal = (1 - u) / (1 + u), (1 + u) ** 2
        static = np.diag([1 / lapse, 1 / conformal, 1 / conformal, 1 / conformal])
        return boost(self.velocity) @ static

    def normal_coordinates(self, points):
        """Normal coordinates at the apex in the comoving frame, in units of T."""
        static = super().normal_coordinates(points)
        g, v = self.gamma, self.velocity
        result = static.copy()
        result[:, 0] = g * (static[:, 0] - v * static[:, 2])
        result[:, 2] = g * (static[:, 2] - v * static[:, 0])
        return result

    def weyl(self):
        """T^2 C_abcd at the apex in the comoving frame (u, radial, motion, z)."""
        frame = boost(self.velocity)
        tensor = np.einsum(
            "ae,bf,cg,dh,efgh->abcd",
            frame,
            frame,
            frame,
            frame,
            static_riemann(self.mass, self.apex_areal_radius),
        )
        return self.proper_time**2 * tensor

    def tidal(self):
        """T^2 E_ij at the apex, axes (radial, direction of motion, z)."""
        return electric_magnetic(self.weyl())[0]

    def magnetic(self):
        """T^2 B_ij = T^2 eps_ikl C_kl0j / 2 at the apex, same axes."""
        return electric_magnetic(self.weyl())[1]

    # -- sampling -------------------------------------------------------------
    def height(self, positions):
        """Coordinate-time extent of the interval above spatial points (direct)."""
        positions = np.asarray(positions, dtype=float)
        return (self.tips[1, 0] - self.travel_time(positions, self.foci[1])) - (
            self.tips[0, 0] + self.travel_time(self.foci[0], positions)
        )

    def sample(self, n, rng):
        """Rejection sample, uniform in proper volume, of the tips' interval.

        Proposal: radius R Beta(3, 2) about the ball's centre, a density
        proportional to 1 - |X - c| / R. Because the index exceeds n_min =
        n(|c| + t_h) on the region and |X - P| + |X - Q| >= 2 |X - c|, the height
        is at most 2 t_h (1 - |X - c| / R); a violated envelope raises, and the
        runner records that case as failed.
        Events closer than BOUNDARY_LAYER (in coordinate time) to the null
        boundary of the tips are not drawn.
        """
        if n < 2:
            raise ValueError("At least two events required")
        fit = self._surrogates(self.foci)[:3]
        ceiling = self.lapse_volume(np.linalg.norm(self.base) - self.radius)[1]
        accepted, count = [], 0
        while count < n:
            k = max(64, int(4 * (n - count)))
            radius = self.radius * rng.beta(3, 2, k)
            direction = rng.normal(size=(k, 3))
            direction /= np.linalg.norm(direction, axis=1)[:, None]
            position = self.base + radius[:, None] * direction
            past, future = np.sqrt(np.maximum(self._evaluate(*fit, position), 0))
            low, high = self.tips[0, 0] + past, self.tips[1, 0] - future
            height = (high - low) / 2 - BOUNDARY_LAYER
            density = self.lapse_volume(np.linalg.norm(position, axis=1))[1]
            ratio = (density / ceiling) * (
                height / (self.half_time * (1 - radius / self.radius))
            )
            if np.max(ratio) > 1 + 1e-9:
                raise ArithmeticError("Rejection envelope violated")
            keep = rng.random(k) < ratio
            time = (high + low) / 2 + (2 * rng.random(k) - 1) * height
            accepted.append(np.column_stack((time, position))[keep])
            count += int(keep.sum())
        return np.concatenate(accepted)[:n]

    # -- quadrature -----------------------------------------------------------
    def directions(self, polar, azimuth):
        """Gauss(cos) x trapezoid(azimuth) directions about the focal axis."""
        axis = self.foci[1] - self.foci[0]
        size = np.linalg.norm(axis)
        axis = axis / size if size > 1e-12 else np.array([1.0, 0.0, 0.0])
        third = np.array([0.0, 0.0, 1.0])
        side = np.cross(third, axis)
        side /= np.linalg.norm(side)
        third = np.cross(axis, side)
        mu, w_mu = legendre.leggauss(polar)
        phi = 2 * np.pi * (np.arange(azimuth) + 0.5) / azimuth
        sine = np.sqrt(1 - mu * mu)
        omega = (
            mu[:, None, None] * axis
            + (sine[:, None] * np.cos(phi))[..., None] * side
            + (sine[:, None] * np.sin(phi))[..., None] * third
        ).reshape(-1, 3)
        weight = np.repeat(w_mu, azimuth) * 2 * np.pi / azimuth
        return omega, weight

    def ray_edges(self, focus, omega, steps=64):
        """Distance along each ray from a focus to the null boundary (direct)."""
        lo = np.zeros(len(omega))
        hi = np.full(len(omega), np.linalg.norm(focus - self.base) + self.radius)
        if np.any(self.height(focus + hi[:, None] * omega) > 0):
            raise ArithmeticError("The interval is not inside its bounding ball")
        for _ in range(steps):
            middle = (lo + hi) / 2
            inside = self.height(focus + middle[:, None] * omega) > 0
            lo, hi = np.where(inside, middle, lo), np.where(inside, hi, middle)
        return (lo + hi) / 2

    def quadrature(self, polar=24, azimuth=32, radial=24):
        """Node sets for the proper-volume measure, one per focus.

        For f(X) g(t) with an antiderivative G of g,
            int f g dvol = sum_Q w f sqrt(-g) G(t_Q - d(X, Q))
                         - sum_P w f sqrt(-g) G(t_P + d(P, X)),
        since the time range above X is [t_P + d(P, X), t_Q - d(X, Q)]. The
        P sum is smooth along rays from P and the Q sum along rays from Q.
        Returns ((points, weights, lower) on P's rays, (points, weights, upper)
        on Q's rays); weights are coordinate volumes.
        """
        omega, w_omega = self.directions(polar, azimuth)
        x, w_r = legendre.leggauss(radial)
        result = []
        for index, focus in enumerate(self.foci):
            edge = self.ray_edges(focus, omega)
            r = edge[:, None] * (x + 1) / 2
            points = (focus + r[..., None] * omega[:, None, :]).reshape(-1, 3)
            weight = ((w_omega * edge / 2)[:, None] * w_r * r * r).ravel()
            if index == 0:
                bound = self.tips[0, 0] + self.travel_time(focus, points)
            else:
                bound = self.tips[1, 0] - self.travel_time(points, focus)
            result.append((points, weight, bound))
        return tuple(result)

    def integrate(self, spatial, antiderivative, nodes=None):
        """Integral of spatial(X) g(t) over the interval's proper volume."""
        (p, wp, lower), (q, wq, upper) = nodes or self.quadrature()
        total = 0.0
        for points, weight, bound, sign in ((q, wq, upper, 1), (p, wp, lower, -1)):
            density = self.lapse_volume(np.linalg.norm(points, axis=1))[1]
            total += sign * float(
                np.sum(weight * density * spatial(points) * antiderivative(bound))
            )
        return total

    def volume(self, nodes=None):
        return self.integrate(lambda x: np.ones(len(x)), lambda t: t, nodes)

    # -- geodesic tools for the checks ----------------------------------------
    def _flow(self, start, velocity, span=1.0):
        state = np.r_[np.asarray(start, float), np.asarray(velocity, float)]
        return solve_ivp(
            self._geodesic, (0, span), state, method="DOP853", rtol=1e-13, atol=1e-15
        ).y[:, -1]

    def metric(self, point):
        u = self.m / np.linalg.norm(point[1:4])
        return np.diag(
            [-(((1 - u) / (1 + u)) ** 2), (1 + u) ** 4, (1 + u) ** 4, (1 + u) ** 4]
        )

    def world_function(self, first, second):
        """Squared geodesic interval sigma(first, second) by shooting (Newton)."""
        first, second = np.asarray(first, float), np.asarray(second, float)
        velocity = second - first
        for _ in range(12):
            miss = self._flow(first, velocity)[:4] - second
            if np.max(abs(miss)) < 1e-13:
                break
            jacobian = np.column_stack(
                [
                    (
                        self._flow(first, velocity + step)[:4]
                        - self._flow(first, velocity - step)[:4]
                    )
                    / 2e-7
                    for step in 1e-7 * np.eye(4)
                ]
            )
            velocity = velocity - np.linalg.solve(jacobian, miss)
        else:
            raise ArithmeticError("Boundary-value shooting did not converge")
        return float(velocity @ self.metric(first) @ velocity)


def _quartic_basis():
    """Ten algebraic Weyl tensors built here from (E, B) pairs, for the checks.

    Independent of the estimator's basis: C_0i0j = E_ij, C_0jkl = eps_klm B_jm,
    C_ijkl = -eps_ijm eps_kln E_mn, with all pair symmetries.
    """
    symmetric = [np.diag([1.0, -1.0, 0.0]), np.diag([0.0, 1.0, -1.0])]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        value = np.zeros((3, 3))
        value[i, j] = value[j, i] = 1.0
        symmetric.append(value)
    result = []
    for channel in range(2):
        for value in symmetric:
            electric = value if channel == 0 else np.zeros((3, 3))
            magnetic = value if channel == 1 else np.zeros((3, 3))
            tensor = np.zeros((4, 4, 4, 4))
            tensor[0, 1:, 0, 1:] = electric
            tensor[1:, 0, 1:, 0] = electric
            tensor[1:, 0, 0, 1:] = -electric
            tensor[0, 1:, 1:, 0] = -electric
            spatial = -np.einsum("ijm,kln,mn->ijkl", EPSILON, EPSILON, electric)
            tensor[1:, 1:, 1:, 1:] = spatial
            mixed = np.einsum("klm,jm->jkl", EPSILON, magnetic)  # C_0jkl
            tensor[0, 1:, 1:, 1:] = mixed
            tensor[1:, 0, 1:, 1:] = -mixed
            tensor[1:, 1:, 0, 1:] = np.transpose(mixed, (1, 2, 0))
            tensor[1:, 1:, 1:, 0] = -np.transpose(mixed, (1, 2, 0))
            result.append(tensor)
    return np.array(result)


def coordinate_riemann(geometry, point, step=2e-3):
    """R_abcd in isotropic coordinates by nested central differences of the metric.

    Shares nothing with the boost: Christoffel symbols from fourth-order
    differences of g, their derivatives from fourth-order differences of those.
    """

    def metric(x):
        u = geometry.m / np.linalg.norm(x[1:])
        return np.diag(
            [-(((1 - u) / (1 + u)) ** 2), (1 + u) ** 4, (1 + u) ** 4, (1 + u) ** 4]
        )

    def derivative(function, x, h):
        out = []
        for a in range(4):
            e = np.eye(4)[a] * h
            out.append(
                (
                    -function(x + 2 * e)
                    + 8 * function(x + e)
                    - 8 * function(x - e)
                    + function(x - 2 * e)
                )
                / (12 * h)
            )
        return np.array(out)  # derivative index first

    def christoffel(x):
        g = metric(x)
        dg = derivative(metric, x, step)  # dg[c, a, b] = d_c g_ab
        lower = 0.5 * (
            np.einsum("cab->abc", dg) + np.einsum("bac->abc", dg) - dg
        )  # Gamma_abc = (d_c g_ab + d_b g_ac - d_a g_bc) / 2
        return np.einsum("ad,dbc->abc", np.linalg.inv(g), lower)

    gamma = christoffel(point)
    d_gamma = derivative(christoffel, point, step)  # d_gamma[c, a, b, d]
    # R^a_bcd = d_c G^a_bd - d_d G^a_bc + G^a_ce G^e_bd - G^a_de G^e_bc
    up = (
        np.einsum("cabd->abcd", d_gamma)
        - np.einsum("dabc->abcd", d_gamma)
        + np.einsum("ace,ebd->abcd", gamma, gamma)
        - np.einsum("ade,ebc->abcd", gamma, gamma)
    )
    return np.einsum("ae,ebcd->abcd", metric(point), up)


def selfcheck():
    """Deterministic checks of tips, frames, truth, relations, sampler and volume."""
    from schwarzschild_diamond import selfcheck as parent_selfcheck

    mass, apex, proper = 1.0, 6.0, 3.0
    v = circular_velocity(mass, apex)
    orbit = MovingSchwarzschildDiamond(mass, apex, proper, v)
    rest = MovingSchwarzschildDiamond(mass, apex, proper, 0.0)
    flat = MovingSchwarzschildDiamond(0.0, apex, proper, v)
    old = SchwarzschildDiamond(mass, apex, proper)
    report = {"circular_velocity": v}

    # 1. At v = 0 the class is the old generator: tips, region, samples, frames.
    same_tips = np.array_equal(rest.tips, old.tips) and rest.radius == old.radius
    draws = [g.sample(4096, np.random.default_rng(20261121)) for g in (rest, old)]
    points = draws[0][:6]
    frame_gap = float(
        np.max(abs(rest.normal_coordinates(points) - old.normal_coordinates(points)))
    )
    truth_gap = float(np.max(abs(rest.tidal() - old.tidal())))
    assert same_tips and np.array_equal(*draws) and frame_gap == 0
    assert truth_gap < 1e-15 and np.max(abs(rest.magnetic())) == 0
    report.update(
        rest_tips_and_ball_equal_old_bitwise=bool(same_tips),
        rest_samples_equal_old_bitwise=True,
        rest_normal_coordinates_equal_old=frame_gap,
        rest_tidal_difference_from_old=truth_gap,
    )

    # 2. The orbit tips lie on the circular geodesic, in the orbital plane.
    rho = np.linalg.norm(orbit.foci, axis=1)
    omega = np.sqrt(mass / apex**3)
    expected_time = proper / 2 / np.sqrt(1 - 3 * mass / apex)
    angle = np.arctan2(orbit.foci[1, 1], orbit.foci[1, 0])
    tip_error = max(
        float(np.max(abs(rho - orbit.apex))),
        abs(orbit.half_time - expected_time),
        abs(angle - omega * expected_time),
        float(np.max(abs(orbit.foci[:, 2]))),
    )
    # The reflected future tip is where the geodesic arrives backwards.
    reversal = float(np.max(abs(orbit._shoot(-orbit.tip_velocity) - orbit.tips[0])))
    assert tip_error < 1e-11 and reversal < 1e-9
    normal = orbit.normal_coordinates(orbit.tips)
    tip_normal = float(np.max(abs(normal - [[-0.5, 0, 0, 0], [0.5, 0, 0, 0]])))
    assert tip_normal < 1e-9
    report.update(
        orbit_tip_error=tip_error,
        orbit_tip_reflection_vs_integration=reversal,
        orbit_tips_normal_coordinate_error=tip_normal,
    )

    # 3. Truth, closed form: Psi = M/r^3, K = gamma^2 v^2.
    psi = proper**2 * mass / apex**3
    k = v * v / (1 - v * v)
    electric, magnetic = orbit.tidal(), orbit.magnetic()
    closed_e = psi * np.diag([-(2 + 3 * k), 1.0, 1 + 3 * k])
    closed_b = np.zeros((3, 3))
    closed_b[0, 2] = closed_b[2, 0] = 3 * np.sqrt(k * (1 + k)) * psi
    closed = max(
        float(np.max(abs(electric - closed_e))), float(np.max(abs(magnetic - closed_b)))
    )
    invariants = (
        float(np.sum(electric**2) - np.sum(magnetic**2) - 6 * psi**2),
        float(np.sum(electric * magnetic)),
        float(np.trace(electric)),
        float(np.max(abs(magnetic - magnetic.T))),
    )
    assert closed < 1e-15 and max(map(abs, invariants)) < 1e-15
    report.update(
        truth_closed_form_difference=closed,
        truth_invariant_residuals=invariants,
        truth_electric=electric.tolist(),
        truth_magnetic=magnetic.tolist(),
    )

    # 4. Truth, independently: coordinate Riemann tensor by finite differences,
    # contracted with the comoving tetrad.
    tetrad = orbit.tetrad()
    coordinate = coordinate_riemann(orbit, np.r_[0.0, orbit.apex, 0.0, 0.0])
    framed = proper**2 * np.einsum(
        "ae,bf,cg,dh,efgh->abcd", tetrad, tetrad, tetrad, tetrad, coordinate
    )
    finite_difference = float(np.max(abs(framed - orbit.weyl())))
    assert finite_difference < 1e-7
    report["truth_vs_finite_difference_riemann"] = finite_difference

    # 5. Truth, operationally: relative motion of particles comoving with the
    # tips, released at the apex. The squared separation in the rest frame of
    # the first particle, averaged over proper times +tau and -tau, is
    # |xi|^2 (1 - tau^2 n.E.n) + O(tau^4). Separations are 1e-4 apart.
    tau, shift = 0.04, 1e-4
    start = np.r_[0.0, orbit.apex, 0.0, 0.0]
    u0 = tetrad[0]

    def acceleration(velocity):
        return orbit._geodesic(0, np.r_[start, velocity])[4:]

    def quadratic(n):
        offset = shift * (n @ tetrad[1:])
        # Parallel transport of u to the displaced event, to first order.
        connection = (acceleration(u0 + offset) - acceleration(u0 - offset)) / 4
        moved = u0 + connection
        total = 0.0
        for sign in (1.0, -1.0):
            one = orbit._flow(start, sign * tau * u0)
            two = orbit._flow(start + offset, sign * tau * moved)
            gap = two[:4] - one[:4]
            g = orbit.metric(one[:4])
            velocity = one[4:] / (sign * tau)
            total += gap @ g @ gap + (gap @ g @ velocity) ** 2
        return (2 - total / shift**2) / (2 * tau * tau) * proper**2

    axes = np.eye(3)
    diagonal = np.array([quadratic(e) for e in axes])
    measured = np.diag(diagonal)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        both = quadratic((axes[i] + axes[j]) / np.sqrt(2))
        measured[i, j] = measured[j, i] = both - (diagonal[i] + diagonal[j]) / 2
    deviation = float(np.max(abs(measured - electric)) / np.linalg.norm(electric))
    assert deviation < 2e-3
    report["electric_from_comoving_particles_relative_error"] = deviation
    report["electric_from_comoving_particles"] = measured.tolist()

    # 6. Truth in the estimator's own convention: the even part of the world
    # function between exp(T lambda x) and exp(T lambda y) is
    # eta(D, D) - C(x, y, x, y) lambda^4 / 3 + O(lambda^6); two scales remove the
    # lambda^6 term. Fitting the quartic over 24 pairs gives the ten components
    # in the normal-coordinate frame, magnetic sign included.
    basis = _quartic_basis()
    stream = np.random.default_rng(20261122)
    pairs = stream.uniform(-0.5, 0.5, size=(24, 2, 4))
    eta = np.diag([-1.0, 1, 1, 1])

    def event(x):
        return orbit._shoot(proper * (x @ tetrad))

    def quartic_values(scale):
        values = []
        for x, y in pairs:
            even = 0.0
            for sign in (1.0, -1.0):
                a, b = sign * scale * x, sign * scale * y
                even += orbit.world_function(event(a), event(b)) / proper**2 / 2
            d = scale * (x - y)
            values.append(3 * (d @ eta @ d - even) / scale**4)
        return np.array(values)

    richardson = (4 * quartic_values(0.05) - quartic_values(0.1)) / 3
    design = np.einsum(
        "pa,pb,pc,pd,kabcd->pk", *pairs.transpose(1, 0, 2)[[0, 1, 0, 1]], basis
    )
    fitted, *_ = np.linalg.lstsq(design, richardson, rcond=None)
    rebuilt = np.einsum("k,kabcd->abcd", fitted, basis)
    world = electric_magnetic(rebuilt)
    world_error = float(
        max(np.max(abs(world[0] - electric)), np.max(abs(world[1] - magnetic)))
        / np.sqrt(np.sum(electric**2) + np.sum(magnetic**2))
    )
    assert np.linalg.matrix_rank(design) == 10 and world_error < 2e-3
    report.update(
        world_function_relative_error=world_error,
        world_function_electric=world[0].tolist(),
        world_function_magnetic=world[1].tolist(),
    )

    # 6b. The frame of normal_coordinates is the tetrad's: events shot along
    # x @ tetrad come back as x. The truth in the runner is rotated with
    # normal_coordinates, the checks above use the tetrad; this ties them.
    probes = np.random.default_rng(20261129).uniform(-0.3, 0.3, size=(6, 4))
    returned = orbit.normal_coordinates(np.array([event(x) for x in probes]))
    frame_link = float(np.max(abs(returned - probes)))
    assert frame_link < 1e-9
    report["normal_coordinates_return_tetrad_vectors"] = frame_link

    # 7. Volume: a small interval exceeds pi T^4 / 24 by T^4 (17 E.E + 3 B.B) /
    # 10080 for the observer along the tips. For the orbit both terms enter, and
    # the boost-invariant E.E - B.B then fixes |E| and |B| separately.
    laws = {}
    for name, speed in (("orbit", v), ("rest", 0.0)):
        small = MovingSchwarzschildDiamond(mass, apex, 0.5, speed)
        e, b = small.tidal(), small.magnetic()
        excess = small.volume(small.quadrature(32, 48, 32)) / (np.pi * 0.5**4 / 24) - 1
        laws[name] = excess / ((17 * np.sum(e * e) + 3 * np.sum(b * b)) / 10080)
    assert all(abs(x - 1) < 1e-3 for x in laws.values())
    report["small_interval_volume_over_curvature_square_law"] = laws

    # 7b. The quadrature: the orbit volume on two grids, and at v = 0 against the
    # old axially symmetric quadrature of the parent class.
    coarse, fine = (
        orbit.volume(orbit.quadrature(*n)) for n in ((24, 32, 24), (32, 48, 32))
    )
    node, measure, height = old.interval_nodes()
    radius = np.linalg.norm(old.base + node, axis=2)
    parent = float(2 * (measure * old.lapse_volume(radius)[1] * height).sum())
    against_parent = abs(rest.volume() / parent - 1)
    convergence = abs(coarse / fine - 1)
    assert convergence < 1e-12 and against_parent < 1e-12
    report.update(
        orbit_volume_grid_change=convergence,
        rest_volume_against_parent_quadrature=against_parent,
    )

    # 8. Relations: surrogate bits equal the direct solver, the order is a
    # partial order, and sampled events lie inside the tips' interval.
    for name, geometry in (("orbit", orbit), ("flat_moving", flat)):
        stream = np.random.default_rng(20261123)
        events = np.vstack([geometry.sample(4096, stream) for _ in range(2)])
        inside = (
            geometry.relations(geometry.tips[0][None], events)[0]
            & geometry.relations(events, geometry.tips[1][None])[0]
        )
        landmarks = np.vstack([geometry.tips[0], events[:32], geometry.tips[1]])
        packed = geometry.packed_relations(landmarks, events)
        direct = geometry.relations(landmarks[:, None], events[None])
        bits = all(
            np.array_equal(np.unpackbits(p, axis=1).astype(bool), d)
            for p, d in zip(packed, direct, strict=True)
        )
        order = geometry.order(
            np.vstack([geometry.tips[0], events[:200], geometry.tips[1]])
        )
        closure = (order.astype(np.int32) @ order.astype(np.int32)) > 0
        assert inside.all() and bits and order[0, 1:].all() and order[:-1, -1].all()
        assert not (order & order.T).any() and not (closure & ~order).any()
        report[f"{name}_surrogate_bits_equal_direct_solver"] = bool(bits)
        report[f"{name}_sampled_events_inside"] = bool(inside.all())
        report[f"{name}_margins_resolved_directly"] = int(geometry.resolved_margins)

    # 9. Region: on a 16 x 24 direction grid sampled at 64 radii, every ray meets
    # the boundary once, inside the bounding ball; the rejection envelope holds
    # on 20,000 probe points.
    omega, _ = orbit.directions(16, 24)
    worst_ratio, crossings = 0.0, 0
    for focus in orbit.foci:
        edge = orbit.ray_edges(focus, omega)
        assert (
            np.max(np.linalg.norm(focus + edge[:, None] * omega - orbit.base, axis=1))
            < orbit.radius
        )
        reach = np.linalg.norm(focus - orbit.base) + orbit.radius
        grid = np.linspace(0, 1, 65)[1:] * reach
        h = orbit.height(focus + grid[None, :, None] * omega[:, None, :])
        sign = h > 0
        crossings = max(
            crossings, int(np.max(np.sum(sign[:, 1:] != sign[:, :-1], axis=1)))
        )
    probe = orbit.base + orbit.radius * np.random.default_rng(20261124).random(
        (20000, 1)
    ) ** (1 / 3) * (lambda d: d / np.linalg.norm(d, axis=1)[:, None])(
        np.random.default_rng(20261125).normal(size=(20000, 3))
    )
    h = np.maximum(orbit.height(probe), 0)
    density = orbit.lapse_volume(np.linalg.norm(probe, axis=1))[1]
    ceiling = orbit.lapse_volume(np.linalg.norm(orbit.base) - orbit.radius)[1]
    r = np.linalg.norm(probe - orbit.base, axis=1)
    worst_ratio = float(
        np.max(
            (density / ceiling) * (h / 2 / (orbit.half_time * (1 - r / orbit.radius)))
        )
    )
    assert crossings == 1 and worst_ratio < 1
    report.update(
        largest_boundary_crossings_per_ray=crossings, largest_envelope_ratio=worst_ratio
    )

    # 10. Sampler: moments of 2e6 events against the quadrature of the measure.
    nodes = orbit.quadrature()
    size = orbit.volume(nodes)
    c = orbit.base
    spatial = (
        ("x", lambda x: x[:, 0] - c[0], lambda t: t),
        ("y", lambda x: x[:, 1] - c[1], lambda t: t),
        ("squared_offset", lambda x: np.sum((x - c) ** 2, axis=1), lambda t: t),
        ("isotropic_radius", lambda x: np.linalg.norm(x, axis=1), lambda t: t),
        ("t_squared", lambda x: np.ones(len(x)), lambda t: t**3 / 3),
        ("t_times_y", lambda x: x[:, 1] - c[1], lambda t: t**2 / 2),
    )
    drawn = orbit.sample(2_000_000, np.random.default_rng(20261126))
    offset = drawn[:, 1:] - c
    samples = {
        "x": offset[:, 0],
        "y": offset[:, 1],
        "squared_offset": np.sum(offset**2, axis=1),
        "isotropic_radius": np.linalg.norm(drawn[:, 1:], axis=1),
        "t_squared": drawn[:, 0] ** 2,
        "t_times_y": drawn[:, 0] * offset[:, 1],
    }
    z = {}
    for name, f, g in spatial:
        target = orbit.integrate(f, g, nodes) / size
        value = samples[name]
        z[name] = float((value.mean() - target) / (value.std() / np.sqrt(len(value))))
    assert max(map(abs, z.values())) < 4.5
    report["sampler_moment_z_scores"] = z
    report["orbit_volume_over_flat"] = size / (np.pi * proper**4 / 24)

    # 11. Flat moving control: Minkowski relations, boosted tips and frames,
    # exact volume, zero curvature.
    events = flat.sample(4096, np.random.default_rng(20261127))
    g, w = flat.gamma, flat.velocity
    expected_tips = np.array(
        [
            [-g * proper / 2, apex, -g * w * proper / 2, 0],
            [g * proper / 2, apex, g * w * proper / 2, 0],
        ]
    )
    shifted = events[:4] - [0, apex, 0, 0]
    lorentz = (
        np.column_stack(
            (
                g * (shifted[:, 0] - w * shifted[:, 2]),
                shifted[:, 1],
                g * (shifted[:, 2] - w * shifted[:, 0]),
                shifted[:, 3],
            )
        )
        / proper
    )
    flat_error = max(
        float(np.max(abs(flat.tips - expected_tips))),
        float(np.max(abs(flat.normal_coordinates(events[:4]) - lorentz))),
        abs(flat.volume(flat.quadrature()) / (np.pi * proper**4 / 24) - 1),
        float(np.max(abs(flat.weyl()))),
    )
    assert flat_error < 1e-11
    report["flat_moving_error"] = flat_error

    # 12. The parent class still passes its own checks.
    parent = parent_selfcheck()
    report["parent_selfcheck_tidal_error"] = parent[
        "tidal_tensor_relative_error_from_released_particles"
    ]
    report["tidal_norm"] = float(np.sqrt(np.sum(electric**2) + np.sum(magnetic**2)))
    report["electric_norm"] = float(np.linalg.norm(electric))
    report["magnetic_norm"] = float(np.linalg.norm(magnetic))
    return report


if __name__ == "__main__":
    import json

    print(json.dumps(selfcheck(), indent=2))
