"""Deterministic landmark-design study for the sixth-degree jet.

Development diagnostic: no pre-specified criterion and no measurement sample.
For one landmark layout and one fixed selection sample it computes

1. the landmarks and the clock mask, by an order-only band rule that
   generalizes ``intrinsic_weyl_joint.select_landmarks`` (designs "connected"
   and "nested"); for the present layout it reproduces that rule exactly;
2. the flat expected clock covariance at fixed N,
   Cov(n_a/N, n_b/N) = (pi_ab - pi_a pi_b)/N, propagated to the clocks
   -sqrt(n/N) by the delta method as in ``read_order_clocks``; pi_ab is
   estimated from a deterministic scrambled-Sobol flat sample of the diamond,
   with the landmarks at their true normal coordinates;
3. the fifth- and sixth-degree continuum fits (curvature-square terms on) on
   exact interval volumes of the half-cycle benchmark wave, which give the
   continuum bias and the nominal 95% statistical radius at that N.

The continuum fit starts at the true gauge-fixed normal coordinates; its
objective, whitening and covariance are those of ``fit_joint``, and the
validation mode checks it against ``fit_joint`` on an archived benchmark run.
Nothing here is a finite-N result: radii and powers are nominal Gaussian
quantities under the flat expected covariance.
"""

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_SEED = 2026111100
SELECTION_N = 4096
REFERENCE_N = 2**24
FALSE_POSITIVE = 1e-3
POWER = 0.9
TARGET_RADIUS = 0.10

# Bands are (t_low, t_high, r_low, r_high) in the flat proxies of
# select_landmarks; the first and last band are the "end" bands.
PRESENT = {
    "name": "present",
    "bands": [
        (-0.36, -0.24, 0.075, 0.2),
        (-0.075, 0.075, 0.075, 0.2),
        (0.24, 0.36, 0.075, 0.2),
    ],
    "score_radius": 0.2,
    "maximum": 64,
    "nested": True,
    "selection_n": SELECTION_N,
}


def layout(name, bands, **options):
    value = dict(PRESENT, name=name, bands=[tuple(b) for b in bands])
    value.update(options)
    return value


# ---------------------------------------------------------------------------
# Selection sample and order-only landmark rule
# ---------------------------------------------------------------------------


def selection_sample(seed, size, wave):
    """The benchmark's selection draw: same RNG stream, tips and relabelling."""
    rng = np.random.default_rng([seed, 0])
    points = np.vstack(([0, 0.5, 0, 0], wave.sample(size, rng), [1, -0.5, 0, 0]))
    permutation = np.r_[
        0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
    ]
    points = points[permutation]
    return points, wave.order(points)


def cached_selection(cache, case, size, wave):
    """Selection sample of case ``case`` (seed BENCHMARK_SEED + case), cached.

    Written under a temporary name and renamed, so that parallel preparation
    never exposes a partial file.
    """
    import os

    path = Path(cache) / f"selection_case{case}_n{size}.npz"
    if path.exists():
        with np.load(path) as saved:
            return saved["points"], saved["order"]
    points, order = selection_sample(BENCHMARK_SEED + case, size, wave)
    partial = path.with_name(f"{path.stem}.{os.getpid()}.partial.npz")
    np.savez(partial, points=points, order=order)
    os.replace(partial, path)
    return points, order


def selection_counts(order):
    order = np.asarray(order, dtype=bool)
    counts = order.astype(np.float32) @ order.astype(np.float32)
    counts = counts.astype(float)
    counts += counts.T.copy()
    return counts


def proxies(counts):
    total = counts[0, -1]
    tau_p, tau_q = np.sqrt(counts[0] / total), np.sqrt(counts[-1] / total)
    times = (tau_p - tau_q) / 2
    radius2 = 0.25 + times**2 - (tau_p + tau_q) / 2
    return times, radius2


def band_index(times, radius2, bands):
    """Index of the band holding each candidate, or -1."""
    index = np.full(len(times), -1)
    for k, (t_low, t_high, r_low, r_high) in enumerate(bands):
        inside = (
            (times >= t_low)
            & (times <= t_high)
            & (radius2 >= r_low**2)
            & (radius2 <= r_high**2)
            & (index < 0)
        )
        index[inside] = k
    return index


def select_bands(counts, spec, maximum, minimum_count=2):
    """Generalized 'connected' rule; equals select_landmarks for PRESENT.

    End bands (first, last) keep the top maximum//3 by related-count times
    capped squared radius. Middle candidates need four related clocks to each
    retained end band; the remaining quota is split equally over the middle
    bands (one middle band: the original rule). Peeling to degree eight and
    the cap, then the identifiability checks of the original rule, extended
    to every middle band.

    Options (absent for the present layout): ``end_fraction`` replaces 1/3
    for the end-band quota; ``connectivity="cross"`` lets a middle band whose
    centre is off zero connect to the FAR end band and to the middle bands on
    the other side of t = 0 instead of to both end bands, so that bands too
    close to one end band can still be retained while every middle landmark
    stays linked across the centre.
    """
    bands = spec["bands"]
    cross = spec.get("connectivity", "ends") == "cross"
    times, radius2 = proxies(counts)
    index = band_index(times, radius2, bands)
    # Optional fifth band entry: that band's score-radius cap (default
    # spec["score_radius"]); a cap of 0 ranks the band by degree alone.
    caps = np.array([b[4] if len(b) > 4 else spec["score_radius"] for b in bands])
    ids = np.flatnonzero(
        (index >= 0) & (counts[0] >= minimum_count) & (counts[-1] >= minimum_count)
    )
    cap2 = caps[index[ids]] ** 2
    radial = np.where(cap2 > 0, np.minimum(radius2[ids], cap2), 1.0)
    related = counts[ids[:, None], ids] >= minimum_count
    score = related.sum(axis=1) * radial
    fraction = spec.get("end_fraction")
    third = maximum // 3 if fraction is None else int(maximum * fraction)
    # Optional ``end_strata``: split each end band into that many equal t
    # strata with equal quotas (remainder to the inner strata), so that the
    # radius score cannot pile every end landmark at the band's inner edge.
    strata = spec.get("end_strata", 1)
    ends = []
    for k in (0, len(bands) - 1):
        choices = np.flatnonzero(index[ids] == k)
        if strata == 1:
            ranked = choices[np.argsort(-score[choices], kind="stable")]
            ends.append(ids[ranked[:third]])
            continue
        edges = np.linspace(bands[k][0], bands[k][1], strata + 1)
        inner_first = range(strata) if k == len(bands) - 1 else range(strata)[::-1]
        part = [third // strata] * strata
        for j, s in enumerate(inner_first):
            if j < third - sum([third // strata] * strata):
                part[s] += 1
        picked = []
        for s in range(strata):
            low, high = edges[s], edges[s + 1]
            inside = choices[
                (times[ids[choices]] >= low)
                & ((times[ids[choices]] < high) | (s == strata - 1))
            ]
            ranked = inside[np.argsort(-score[inside], kind="stable")]
            picked.append(ids[ranked[: part[s]]])
        ends.append(np.concatenate(picked))
    middle_bands = list(range(1, len(bands) - 1))
    centres = {k: (bands[k][0] + bands[k][1]) / 2 for k in range(len(bands))}

    def targets(k, end_sets, member):
        """Landmark sets a middle band must connect to (four clocks each)."""
        if not cross or abs(centres[k]) < 1e-9:
            return end_sets
        far = end_sets[1] if centres[k] < 0 else end_sets[0]
        other = [j for j in middle_bands if centres[j] * centres[k] < 0]
        return [far, member[np.isin(index[member], other)]]

    quota = maximum - 2 * third
    shares = [quota // len(middle_bands)] * len(middle_bands)
    # The remainder goes to the most central band(s).
    centre = np.argsort([abs(sum(bands[k][:2])) for k in middle_bands], kind="stable")
    for j in range(quota - sum(shares)):
        shares[centre[j % len(centre)]] += 1
    middle_parts = []
    for k, share in zip(middle_bands, shares, strict=True):
        middle = ids[index[ids] == k]
        degrees = np.column_stack(
            [
                (counts[middle[:, None], end] >= minimum_count).sum(axis=1)
                for end in targets(k, ends, ids)
            ]
        )
        connected = np.all(degrees >= 4, axis=1)
        middle, degrees = middle[connected], degrees[connected]
        band_cap2 = caps[k] ** 2
        middle_score = degrees.min(axis=1) * (
            np.minimum(radius2[middle], band_cap2) if band_cap2 > 0 else 1.0
        )
        ranked = np.argsort(-middle_score, kind="stable")
        middle_parts.append(middle[ranked[:share]])
    ids = np.r_[ends[0], np.concatenate(middle_parts), ends[1]]
    while len(ids):
        graph = counts[ids[:, None], ids] >= minimum_count
        degrees = graph.sum(axis=1)
        if degrees.min() >= 8 and len(ids) <= maximum:
            break
        score = degrees + np.log1p(counts[ids[:, None], ids]).sum(axis=1) * 1e-8
        ids = np.delete(ids, int(np.argmin(score)))
    if len(ids) < 16:
        raise ValueError("Insufficient mutually informative landmarks")
    groups = [ids[index[ids] == k] for k in range(len(bands))]
    if len(groups[0]) < 16 or len(groups[-1]) < 16:
        raise ValueError("End bands need at least 16 retained landmarks each")
    for k in middle_bands:
        group = groups[k]
        if len(group) < 8:
            raise ValueError("Every middle band needs at least 8 retained landmarks")
        if any(
            np.any((counts[group[:, None], end] >= minimum_count).sum(axis=1) < 4)
            for end in targets(k, (groups[0], groups[-1]), ids)
        ):
            raise ValueError("Middle landmarks need four retained clocks to each side")
    return np.r_[0, ids, len(counts) - 1]


def select_layout(counts, spec):
    maximum = spec["maximum"]
    if not spec.get("nested", True):
        return select_bands(counts, spec, maximum)
    base = select_bands(counts, spec, maximum)
    expanded = select_bands(counts, spec, 2 * maximum)
    extra = expanded[~np.isin(expanded, base)]
    return np.r_[base[:-1], extra, base[-1]]


def clock_set(order, counts, ids, minimum_count=2):
    restricted = counts[ids[:, None], ids]
    mask = restricted >= minimum_count
    mask[0, -1] = mask[-1, 0] = False
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    forward = np.asarray(order)[ids[:, None], ids][pairs[:, 0], pairs[:, 1]]
    a = np.where(forward, pairs[:, 0], pairs[:, 1])
    b = np.where(forward, pairs[:, 1], pairs[:, 0])
    return mask, pairs, a, b


# ---------------------------------------------------------------------------
# Flat expected clock covariance
# ---------------------------------------------------------------------------


def diamond_sobol(log2_size, seed):
    """Measure-preserving map of a scrambled Sobol sequence onto the diamond.

    |x| has density proportional to rho^2 (1 - 2 rho) on [0, 1/2], CDF
    32 rho^3 - 48 rho^4; t is uniform on |t| <= 1/2 - rho; the direction is
    uniform on the sphere.
    """
    from scipy.stats import qmc

    u = qmc.Sobol(d=4, scramble=True, rng=np.random.default_rng(seed)).random_base2(
        log2_size
    )
    grid = np.linspace(0, 0.5, 1 << 17)
    cdf = 32 * grid**3 - 48 * grid**4
    rho = np.interp(u[:, 0], cdf, grid)
    t = (2 * u[:, 1] - 1) * (0.5 - rho)
    z = 1 - 2 * u[:, 2]
    phi = 2 * np.pi * u[:, 3]
    s = np.sqrt(np.maximum(1 - z * z, 0))
    return np.column_stack((t, rho * s * np.cos(phi), rho * s * np.sin(phi), rho * z))


def flat_membership_moments(points, a, b, log2_size=20, seed=2026112101, chunk=16384):
    """First and second membership moments of the clocks in flat spacetime."""
    events = diamond_sobol(log2_size, seed)
    k = len(a)
    first = np.zeros(k)
    second = np.zeros((k, k))
    for start in range(0, len(events), chunk):
        e = events[start : start + chunk]
        dt = e[None, :, 0] - points[:, None, 0]
        dx2 = ((e[None, :, 1:] - points[:, None, 1:]) ** 2).sum(axis=-1)
        light = dt * dt > dx2
        after, before = light & (dt > 0), light & (dt < 0)
        member = (after[a] & before[b]).astype(np.float32)
        first += member.sum(axis=1)
        # Products of 0/1 summed over at most 2**14 events are exact in float32.
        second += member @ member.T
    return first, second, len(events)


def flat_clock_covariance(first, second, size, n):
    p = first / size
    if np.any(p <= 0):
        raise ValueError("A clock has no quasi-random event")
    indicator = (second - size * np.outer(p, p)) / (size - 1)
    derivative = -0.5 / np.sqrt(p)
    return indicator * np.outer(derivative, derivative) / n, p


# ---------------------------------------------------------------------------
# Continuum fit from the truth start (objective of fit_joint)
# ---------------------------------------------------------------------------


def continuum_fit(counts, clock_mask, clock_covariance, model_order, start):
    """fit_joint's objective, whitening and covariance, started at the truth.

    ``start`` are gauge-fixed normal coordinates of all points (tips first and
    last). Only the curved Levenberg-Marquardt stage is run.
    """
    from intrinsic_weyl import WEYL
    from intrinsic_weyl_joint import model
    from scipy.optimize import least_squares

    curvature_parameters = {5: 34, 6: 76}[model_order]
    counts = np.asarray(counts, dtype=float)
    n = len(counts)
    available = np.triu(clock_mask, 1)
    available[0, -1] = False
    pairs = np.column_stack(np.nonzero(available))
    target = -np.sqrt(counts[pairs[:, 0], pairs[:, 1]] / counts[0, -1])
    covariance = np.asarray(clock_covariance)
    scales = np.sqrt(np.diag(covariance))
    eigenvalues, vectors = np.linalg.eigh(covariance / scales[:, None] / scales)
    keep = eigenvalues > eigenvalues[-1] * 1e-8
    whiten = (vectors[:, keep] / np.sqrt(eigenvalues[keep])).T / scales
    mask = np.ones((n, 4), dtype=bool)
    mask[[0, -1]] = False
    mask[1, 2:] = False
    mask[2, 3] = False
    parameter_index = np.full((n, 4), -1)
    parameter_index[mask] = np.arange(mask.sum())
    dimension = int(mask.sum())

    def unpack(parameters):
        points = np.zeros((n, 4))
        points[0, 0], points[-1, 0] = -0.5, 0.5
        points[mask] = parameters[:dimension]
        return points, parameters[dimension:]

    def residual(parameters):
        value = model(*unpack(parameters), pairs, square_correction=True)[0] - target
        return whiten @ value

    def jacobian(parameters):
        _, dl, dr, dc = model(*unpack(parameters), pairs, square_correction=True)
        result = np.zeros((len(pairs), len(parameters)))
        for side, derivative in enumerate((dl, dr)):
            for component in range(4):
                index = parameter_index[pairs[:, side], component]
                valid = index >= 0
                result[np.flatnonzero(valid), index[valid]] = derivative[
                    valid, component
                ]
        result[:, dimension:] = dc
        return whiten @ result

    initial = np.r_[start[mask], np.zeros(curvature_parameters)]
    curved = least_squares(
        residual,
        initial,
        jac=jacobian,
        method="lm",
        max_nfev=600,
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-10,
    )
    points, c = unpack(curved.x)
    jac = jacobian(curved.x)
    geometry = jac[:, :dimension]
    curvature = jac[:, dimension : dimension + 10]
    nuisance = np.column_stack((geometry, jac[:, dimension + 10 :]))
    projected = (
        curvature - nuisance @ np.linalg.lstsq(nuisance, curvature, rcond=None)[0]
    )
    singular = np.linalg.svd(projected, compute_uv=False)
    influence = np.linalg.pinv(jac)[dimension : dimension + 10]
    tensor = np.einsum("k,kabcd->abcd", c[:10], WEYL)
    return {
        "coefficients": c[:10].tolist(),
        "electric_eigenvalues": np.linalg.eigvalsh(tensor[0, 1:, 0, 1:]).tolist(),
        "coordinates_inferred_from_counts": points.tolist(),
        "landmarks": n - 2,
        "measured_intervals": len(pairs),
        "coordinate_rank": int(np.linalg.matrix_rank(geometry)),
        "coordinate_parameters": dimension,
        "projected_curvature_rank": int(np.linalg.matrix_rank(projected)),
        "projected_curvature_singular_values": singular.tolist(),
        "curved_cost": float(curved.cost),
        "optimizer_success": bool(curved.success),
        "function_evaluations": int(curved.nfev),
        "boundary_excess": float(max(0, np.max(abs(points[:, 0])) - 0.5)),
        "model_order": model_order,
        "nuisance_rank": int(np.linalg.matrix_rank(nuisance)),
        "nuisance_parameters": nuisance.shape[1],
        "nominal_component_covariance": (influence @ influence.T).tolist(),
        "retained_clock_covariance_rank": len(whiten),
    }


def required_noncentrality():
    from scipy.optimize import brentq
    from scipy.stats import chi2, ncx2

    threshold = float(chi2.ppf(1 - FALSE_POSITIVE, 10))
    value = brentq(lambda lam: ncx2.sf(threshold, 10, lam) - POWER, 1e-6, 1e4)
    return threshold, value


def evaluate(fit, expected, amplitude, n):
    from intrinsic_weyl_joint import TIDAL_ROOT, confidence_ball, tidal_distance

    truth_norm = 2 * amplitude
    ball = confidence_ball(fit)
    coefficients = np.asarray(fit["coefficients"])
    covariance = np.asarray(fit["nominal_component_covariance"])
    noncentrality = float(coefficients @ np.linalg.solve(covariance, coefficients))
    threshold, needed = required_noncentrality()
    radius = ball["statistical_radius"] / truth_norm
    # TIDAL_ROOT is block diagonal: five electric, then five magnetic.
    tidal = TIDAL_ROOT @ covariance @ TIDAL_ROOT.T
    values, vectors = np.linalg.eigh(tidal)
    return {
        "n_times_trace_electric": n * float(np.trace(tidal[:5, :5])),
        "n_times_trace_magnetic": n * float(np.trace(tidal[5:, 5:])),
        "worst_direction_electric_share": float(np.sum(vectors[:5, -1] ** 2)),
        "relative_continuum_bias": tidal_distance(
            ball["coefficients_xy_positive"], expected
        )
        / truth_norm,
        "relative_statistical_radius_95": radius,
        "continuum_noncentrality": noncentrality,
        "events_for_radius_10pct": n * (radius / TARGET_RADIUS) ** 2,
        "events_for_90pct_detection": n * needed / noncentrality,
        "detection_threshold": threshold,
        "noncentrality_for_90pct": needed,
        "projected_curvature_rank": fit["projected_curvature_rank"],
        "nuisance_rank": fit["nuisance_rank"],
        "nuisance_parameters": fit["nuisance_parameters"],
        "optimizer_success": fit["optimizer_success"],
        "function_evaluations": fit["function_evaluations"],
    }


def expected_coefficients(normal, amplitude):
    from intrinsic_weyl_joint import gauge
    from intrinsic_weyl_joint_audit import aligned_tensor_error

    return np.asarray(
        aligned_tensor_error(
            {
                "coefficients": np.zeros(10),
                "coordinates_inferred_from_counts": gauge(normal),
            },
            normal,
            amplitude,
        )["aligned_truth_coefficients"]
    )


# ---------------------------------------------------------------------------
# One layout on one selection sample
# ---------------------------------------------------------------------------


def study(spec, points, order, counts, wave, amplitude, log2_size, n=REFERENCE_N):
    from intrinsic_weyl_joint import gauge

    begin = time.time()
    ids = select_layout(counts, spec)
    mask, pairs, a, b = clock_set(order, counts, ids)
    landmarks = points[ids]
    normal = wave.normal_coordinates(landmarks)
    volume = wave.volumes(landmarks)
    usable = volume[pairs[:, 0], pairs[:, 1]] > 0
    if not np.all(usable):
        raise ArithmeticError("A selected clock has no exact wave volume")
    times = normal[1:-1, 0]
    radii = np.linalg.norm(normal[1:-1, 1:], axis=1)
    selection_seconds = time.time() - begin
    begin = time.time()
    first, second, size = flat_membership_moments(normal, a, b, log2_size)
    covariance, p = flat_clock_covariance(first, second, size, n)
    covariance_seconds = time.time() - begin
    exact_p = (
        (normal[b, 0] - normal[a, 0]) ** 2
        - ((normal[b, 1:] - normal[a, 1:]) ** 2).sum(axis=1)
    ) ** 2
    expected = expected_coefficients(normal, amplitude)
    start = gauge(normal)
    continuum = n * volume / volume[0, -1]
    result = {
        "layout": spec,
        "landmarks": len(ids) - 2,
        "clocks": len(pairs),
        "chord_tt_spread": chord_spread(normal, a, b),
        "landmark_t_range": [float(times.min()), float(times.max())],
        "landmark_r_range": [float(radii.min()), float(radii.max())],
        "landmark_t_quantiles": np.quantile(times, [0.1, 0.5, 0.9]).tolist(),
        "landmark_r_quantiles": np.quantile(radii, [0.1, 0.5, 0.9]).tolist(),
        "band_counts": [
            int(np.sum((times >= lo - 0.03) & (times <= hi + 0.03)))
            for lo, hi, *_ in spec["bands"]
        ],
        "clock_probability_range": [float(p.min()), float(p.max())],
        "quasi_random_vs_exact_flat_probability_max_relative": float(
            np.max(abs(p / exact_p - 1))
        ),
        "selection_seconds": round(selection_seconds, 1),
        "covariance_seconds": round(covariance_seconds, 1),
        "reference_n": n,
    }
    for model_order in (5, 6):
        begin = time.time()
        fit = continuum_fit(continuum, mask, covariance, model_order, start)
        row = evaluate(fit, expected, amplitude, n)
        row["seconds"] = round(time.time() - begin, 1)
        result[f"order_{model_order}"] = row
    result["radius_ratio_6_to_5"] = (
        result["order_6"]["relative_statistical_radius_95"]
        / result["order_5"]["relative_statistical_radius_95"]
    )
    result["principal_angles"] = principal_angles(normal, mask, covariance)
    result["memory"] = {
        "relation_bits_per_event": 2 * len(ids),
        "clock_covariance_gib": len(pairs) ** 2 * 8 / 2**30,
    }
    return result


def ideal_geometry(config, wave, amplitude, log2_size, n=REFERENCE_N):
    """Mechanism test on landmarks placed directly (NOT an order-only design).

    Bands are (t_centre, half_width, r_low, r_high, count) in flat units;
    directions are isotropic. Points are put in the Brinkmann chart through
    the flat map u' = t + z, v = z - t, transverse = sqrt(2) * spatial. A clock
    is a pair whose exact flat expected selection count nsel * tau^4 is at
    least two, as the selection-count rule gives on average; every landmark is
    kept inside the tip-clock floor.
    """
    from intrinsic_weyl_joint import gauge

    rng = np.random.default_rng(config.get("seed", 5))
    nsel = config.get("nsel", SELECTION_N)
    tau_min = (2 / nsel) ** 0.25
    times, radii = [], []
    for centre, width, r_low, r_high, count in config["bands"]:
        t = centre + rng.uniform(-width, width, count)
        r = np.sqrt(rng.uniform(r_low**2, r_high**2, count))
        ceiling = np.sqrt(np.maximum((0.5 - abs(t)) ** 2 - (1.02 * tau_min) ** 2, 0))
        times.append(t)
        radii.append(np.minimum(r, 0.97 * ceiling))
    t, r = np.concatenate(times), np.concatenate(radii)
    direction = rng.normal(size=(len(t), 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    flat = np.column_stack((t, r[:, None] * direction))
    brinkmann = np.vstack(
        (
            [0, 0.5, 0, 0],
            np.column_stack(
                (
                    flat[:, 0] + flat[:, 3] + 0.5,
                    flat[:, 3] - flat[:, 0],
                    np.sqrt(2) * flat[:, 1],
                    np.sqrt(2) * flat[:, 2],
                )
            ),
            [1, -0.5, 0, 0],
        )
    )
    normal = wave.normal_coordinates(brinkmann)
    dt = normal[None, :, 0] - normal[:, None, 0]
    dx2 = ((normal[None, :, 1:] - normal[:, None, 1:]) ** 2).sum(axis=-1)
    tau2 = np.where(dt > 0, dt * dt - dx2, -1.0)
    tau = np.sqrt(np.maximum(np.maximum(tau2, tau2.T), 0))
    mask = tau >= tau_min
    np.fill_diagonal(mask, False)
    mask[0, -1] = mask[-1, 0] = False
    if not (np.all(mask[0, 1:-1]) and np.all(mask[1:-1, -1])):
        raise ArithmeticError("A placed landmark lacks a tip clock")
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    forward = dt[pairs[:, 0], pairs[:, 1]] > 0
    a = np.where(forward, pairs[:, 0], pairs[:, 1])
    b = np.where(forward, pairs[:, 1], pairs[:, 0])
    volume = wave.volumes(brinkmann)
    first, second, size = flat_membership_moments(normal, a, b, log2_size)
    covariance, _ = flat_clock_covariance(first, second, size, n)
    expected = expected_coefficients(normal, amplitude)
    start = gauge(normal)
    result = {
        "ideal": config,
        "landmarks": len(t),
        "clocks": len(pairs),
        "reference_n": n,
    }
    for model_order in (5, 6):
        fit = continuum_fit(
            n * volume / volume[0, -1], mask, covariance, model_order, start
        )
        result[f"order_{model_order}"] = evaluate(fit, expected, amplitude, n)
    result["radius_ratio_6_to_5"] = (
        result["order_6"]["relative_statistical_radius_95"]
        / result["order_5"]["relative_statistical_radius_95"]
    )
    return result


def chord_weight(x, y):
    """Weight of the second-derivative term in intrinsic_weyl_joint: mm/20 + DD/240."""
    m, step = (x + y) / 2, y - x
    return (
        np.einsum("...e,...f->...ef", m, m) / 20
        + np.einsum("...e,...f->...ef", step, step) / 240
    )


def principal_angles(normal, clock_mask, clock_covariance):
    """Collinearity of the Weyl and second-derivative columns, nuisance removed.

    Whitened sixth-degree Jacobian at the true gauge-fixed coordinates and zero
    curvature. After projecting out coordinates and the 24 gradients, the
    cosines of the principal angles between the ten Weyl columns and the 42
    second-derivative columns give the variance inflation 1/(1 - cos^2) of the
    sixth- relative to a jet with the second derivatives known. For the most
    collinear directions, the share of each derivative pair (e, f).
    """
    from intrinsic_weyl import SYMMETRIC_PAIRS, second_derivative_basis
    from intrinsic_weyl_joint import gauge, model

    points = gauge(normal)
    n = len(points)
    pairs = np.column_stack(np.nonzero(np.triu(clock_mask, 1)))
    _, dl, dr, dc = model(points, np.zeros(76), pairs, square_correction=True)
    free = np.ones((n, 4), dtype=bool)
    free[[0, -1]] = False
    free[1, 2:] = False
    free[2, 3] = False
    index = np.full((n, 4), -1)
    index[free] = np.arange(free.sum())
    geometry = np.zeros((len(pairs), int(free.sum())))
    for side, derivative in enumerate((dl, dr)):
        for component in range(4):
            column = index[pairs[:, side], component]
            valid = column >= 0
            geometry[np.flatnonzero(valid), column[valid]] = derivative[
                valid, component
            ]
    covariance = np.asarray(clock_covariance)
    scales = np.sqrt(np.diag(covariance))
    values, vectors = np.linalg.eigh(covariance / scales[:, None] / scales)
    keep = values > values[-1] * 1e-8
    whiten = (vectors[:, keep] / np.sqrt(values[keep])).T / scales
    nuisance, _ = np.linalg.qr(whiten @ np.column_stack((geometry, dc[:, 10:34])))
    weyl = whiten @ dc[:, :10]
    second = whiten @ dc[:, 34:]
    weyl -= nuisance @ (nuisance.T @ weyl)
    second -= nuisance @ (nuisance.T @ second)
    q_weyl, _ = np.linalg.qr(weyl)
    q_second, r_second = np.linalg.qr(second)
    _, cosines, right = np.linalg.svd(q_weyl.T @ q_second)
    directions = np.linalg.solve(r_second, right.T)
    packed = second_derivative_basis().reshape(10, 10, 42)
    names = ["".join("txyz"[i] for i in pair) for pair in SYMMETRIC_PAIRS]
    shares = []
    for j in range(3):
        part = (np.einsum("skj,j->sk", packed, directions[:, j]) ** 2).sum(axis=1)
        shares.append(
            dict(zip(names, (part / part.sum()).round(4).tolist(), strict=True))
        )
    return {
        "cosines": cosines.tolist(),
        "inflation": (1 / (1 - cosines**2)).tolist(),
        "derivative_pair_shares_top3": shares,
    }


def chord_spread(normal, a, b):
    """Coefficient of variation of the time-time chord weight across clocks.

    Weighted by |x^y|^4, the square of the clock's bivector size, which is how
    strongly the clock responds to curvature. Tip clocks (one endpoint a tip)
    respond to the electric part only and are reported separately.
    """
    x, y = normal[a], normal[b]
    weight = chord_weight(x, y)[:, 0, 0]
    bivector = np.stack(
        [
            x[:, i] * y[:, j] - x[:, j] * y[:, i]
            for i in range(4)
            for j in range(i + 1, 4)
        ],
        axis=1,
    )
    signal = (bivector**2).sum(axis=1) ** 2
    tip = (a == 0) | (b == len(normal) - 1)
    result = {}
    for name, chosen in (("landmark_clocks", ~tip), ("tip_clocks", tip)):
        w, s = weight[chosen], signal[chosen]
        mean = np.average(w, weights=s)
        spread = np.sqrt(np.average((w - mean) ** 2, weights=s))
        result[name] = {"clocks": int(chosen.sum()), "cv": float(spread / mean)}
    return result


def selfcheck():
    """Deterministic checks that need no archive."""
    from intrinsic_weyl_joint import select_landmarks

    rng = np.random.default_rng(2026112102)
    x, y = rng.normal(size=(2, 4))
    nodes, weights = np.polynomial.legendre.leggauss(8)
    nodes, weights = (nodes + 1) / 2, weights / 2
    chord = sum(
        w * np.outer(x + s * (y - x), x + s * (y - x))
        for s, w in zip(nodes, weights, strict=True)
    )
    identity = float(np.max(abs(chord_weight(x, y) - chord / 20)))
    assert identity < 1e-15

    # Three slices: the inter-slice time-time weights coincide iff t1+t2+t3 = 0.
    def tt(a, b):
        return chord_weight(np.array([a, 0, 0, 0]), np.array([b, 0, 0, 0]))[0, 0]

    symmetric = [tt(-0.3, 0.0), tt(-0.3, 0.3), tt(0.0, 0.3)]
    shifted = [tt(-0.3, 0.1), tt(-0.3, 0.3), tt(0.1, 0.3)]
    assert np.ptp(symmetric) < 1e-17 and np.ptp(shifted) > 1e-4
    events = diamond_sobol(16, 7)
    inside = np.all(np.linalg.norm(events[:, 1:], axis=1) <= 0.5 - abs(events[:, 0]))
    moments = [
        float(np.mean(events[:, 0] ** 2)),
        float(np.mean((events[:, 1:] ** 2).sum(1))),
    ]
    assert inside and abs(moments[0] - 1 / 60) < 1e-4 and abs(moments[1] - 0.1) < 1e-3
    # Flat selection sample with the exact flat order: the generalized rule
    # equals the repository rule for the present layout.
    sample = diamond_sobol(12, 11)
    points = np.vstack(([-0.5, 0, 0, 0], sample, [0.5, 0, 0, 0]))
    dt = points[None, :, 0] - points[:, None, 0]
    dx2 = ((points[None, :, 1:] - points[:, None, 1:]) ** 2).sum(axis=-1)
    order = (dt > 0) & (dt * dt > dx2)
    counts = selection_counts(order)
    same = []
    for maximum in (64, 128):
        mine = select_bands(counts, PRESENT, maximum)
        theirs = select_landmarks(counts, maximum, 2, "connected")
        same.append(bool(np.array_equal(mine, theirs)))
    nested = select_layout(counts, PRESENT)
    same.append(bool(np.array_equal(nested, select_landmarks(counts, 64, 2, "nested"))))
    assert all(same)
    # Quasi-random membership of a clock against exact flat probability tau^4.
    landmarks = points[nested]
    _, _, a, b = clock_set(order, counts, nested)
    exact = (
        (landmarks[b, 0] - landmarks[a, 0]) ** 2
        - ((landmarks[b, 1:] - landmarks[a, 1:]) ** 2).sum(axis=1)
    ) ** 2
    largest = np.argsort(-exact, kind="stable")[:64]
    first, _, size = flat_membership_moments(landmarks, a[largest], b[largest], 16)
    worst = float(np.max(abs(first / size / exact[largest] - 1)))
    assert worst < 0.02
    return {
        "chord_weight_identity_error": identity,
        "three_slice_tt_weights_symmetric": symmetric,
        "three_slice_tt_weights_shifted": shifted,
        "diamond_moments_t2_r2": moments,
        "diamond_moments_exact": [1 / 60, 0.1],
        "generalized_rule_equals_repository_rule": same,
        "flat_selection_landmarks": len(nested) - 2,
        "quasi_random_clock_probability_worst_relative_error_2e16": worst,
    }


def validate(archive, cache, wave, amplitude, log2_sizes=(20, 22)):
    """Replay against an archived benchmark run (cases 0 and 1, N = 2**24)."""
    from intrinsic_weyl_joint import gauge, read_order_clocks, selection

    rows = []
    for case in (0, 1):
        points, order = cached_selection(cache, case, SELECTION_N, wave)
        name = archive / f"case_{case}_n{REFERENCE_N}.npz"
        with np.load(name) as saved:
            archived_order, archived_ids = (
                saved["selection_order"],
                saved["selected_ids"],
            )
        with np.load(archive / f"truth_{case}.npz") as saved:
            truth = saved["points"]
        counts = selection_counts(order)
        ids = select_layout(counts, PRESENT)
        row = {
            "case": case,
            "selection_order_equal": bool(np.array_equal(order, archived_order)),
            "repository_rule_equal": bool(
                np.array_equal(selection(order, design="nested")[0], archived_ids)
            ),
            "generalized_rule_equal": bool(np.array_equal(ids, archived_ids)),
            "truth_points_equal": bool(np.array_equal(points[ids], truth)),
        }
        clocks = cache / f"run_clocks_case{case}_n{REFERENCE_N}.npz"
        if clocks.exists():
            with np.load(clocks) as saved:
                run_counts, run_mask, run_covariance = (
                    saved["counts"],
                    saved["mask"],
                    saved["cov"],
                )
        else:
            run_counts, run_mask, run_covariance, _, _ = read_order_clocks(name)
            np.savez(clocks, counts=run_counts, mask=run_mask, cov=run_covariance)
        mask, pairs, a, b = clock_set(order, counts, ids)
        row["clock_mask_equal"] = bool(np.array_equal(mask, run_mask))
        landmarks = points[ids]
        normal = wave.normal_coordinates(landmarks)
        volume = wave.volumes(landmarks)
        expected = expected_coefficients(normal, amplitude)
        start = gauge(normal)
        continuum = REFERENCE_N * volume / volume[0, -1]
        for model_order in (5, 6):
            fit = continuum_fit(continuum, mask, run_covariance, model_order, start)
            result = evaluate(fit, expected, amplitude, REFERENCE_N)
            row[f"run_covariance_order_{model_order}"] = {
                k: result[k]
                for k in (
                    "relative_continuum_bias",
                    "relative_statistical_radius_95",
                    "events_for_90pct_detection",
                )
            }
        for log2_size in log2_sizes:
            first, second, size = flat_membership_moments(normal, a, b, log2_size)
            covariance, _ = flat_clock_covariance(first, second, size, REFERENCE_N)
            for model_order in (5, 6):
                fit = continuum_fit(continuum, mask, covariance, model_order, start)
                result = evaluate(fit, expected, amplitude, REFERENCE_N)
                row[f"flat_qmc{log2_size}_order_{model_order}"] = {
                    k: result[k]
                    for k in (
                        "relative_continuum_bias",
                        "relative_statistical_radius_95",
                        "events_for_90pct_detection",
                    )
                }
        rows.append(row)
        print(json.dumps(row), flush=True)
    return rows


def code_state():
    commit = subprocess.run(
        ["git", "rev-parse", "--short=7", "HEAD"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    return commit, not dirty


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--layouts", type=Path)
    parser.add_argument("--cases", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--qmc-log2", type=int, default=20)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--validate",
        type=Path,
        help="benchmark archive folder; replay the checks against that run",
    )
    parser.add_argument(
        "--ideal",
        type=Path,
        help="JSON list of idealized geometries (mechanism test, not order-only)",
    )
    parser.add_argument(
        "--prepare",
        type=int,
        help="only cache the selection samples of this size for --cases",
    )
    args = parser.parse_args()
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from profile_plane_wave import half_cycle

    amplitude = benchmark_amplitude()
    wave = half_cycle(amplitude)
    commit, clean = code_state()
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.ideal:
        rows = []
        for config in json.loads(args.ideal.read_text(encoding="utf-8")):
            row = ideal_geometry(config, wave, float(amplitude), args.qmc_log2)
            rows.append(row)
            print(
                config["name"],
                row["landmarks"],
                row["clocks"],
                f"r5 {row['order_5']['relative_statistical_radius_95']:.3f}",
                f"r6 {row['order_6']['relative_statistical_radius_95']:.3f}",
                f"ratio {row['radius_ratio_6_to_5']:.2f}",
                flush=True,
            )
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(
                    {
                        "diagnostic": "idealized landmark geometries; not order-only",
                        "code_version": commit,
                        "worktree_clean": clean,
                        "python": platform.python_version(),
                        "numpy": np.__version__,
                        "qmc_log2": args.qmc_log2,
                        "rows": rows,
                    },
                    indent=1,
                )
                + "\n",
                encoding="utf-8",
            )
        return
    if args.validate:
        result = {
            "diagnostic": "landmark design harness: checks against the benchmark run",
            "code_version": commit,
            "worktree_clean": clean,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "selfcheck": selfcheck(),
            "archive": str(args.validate),
            "rows": validate(args.validate, args.cache, wave, float(amplitude)),
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
        return
    if args.prepare:
        for case in args.cases:
            cached_selection(args.cache, case, args.prepare, wave)
        return
    specs = json.loads(args.layouts.read_text(encoding="utf-8"))
    rows = []
    if args.output.exists():
        rows = json.loads(args.output.read_text(encoding="utf-8"))["rows"]
    done = {(r["case"], r["layout"]["name"]) for r in rows}
    samples = {}
    for spec in specs:
        spec = dict(PRESENT, **spec)
        for case in args.cases:
            if (case, spec["name"]) in done:
                continue
            key = (case, spec["selection_n"])
            if key not in samples:
                points, order = cached_selection(
                    args.cache, case, spec["selection_n"], wave
                )
                samples = {key: (points, order, selection_counts(order))}
            points, order, counts = samples[key]
            try:
                row = study(
                    spec, points, order, counts, wave, float(amplitude), args.qmc_log2
                )
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                row = {"layout": spec, "failure": f"{type(error).__name__}: {error}"}
            row["case"] = case
            rows.append(row)
            summary = {
                "diagnostic": "landmark design study; deterministic; no criterion",
                "code_version": commit,
                "worktree_clean": clean,
                "python": platform.python_version(),
                "numpy": np.__version__,
                "qmc_log2": args.qmc_log2,
                "benchmark_amplitude": str(amplitude),
                "rows": rows,
            }
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(summary, indent=1) + "\n", encoding="utf-8"
            )
            if "failure" in row:
                print(case, spec["name"], row["failure"], flush=True)
            else:
                print(
                    case,
                    spec["name"],
                    row["landmarks"],
                    row["clocks"],
                    f"r5 {row['order_5']['relative_statistical_radius_95']:.3f}",
                    f"r6 {row['order_6']['relative_statistical_radius_95']:.3f}",
                    f"ratio {row['radius_ratio_6_to_5']:.2f}",
                    f"bias5 {row['order_5']['relative_continuum_bias']:.4f}",
                    f"bias6 {row['order_6']['relative_continuum_bias']:.4f}",
                    flush=True,
                )


if __name__ == "__main__":
    sys.exit(main())
