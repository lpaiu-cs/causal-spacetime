"""Deterministic multi-window study of the fifth-degree truncation bias.

Development diagnostic: no pre-specified criterion, no sampling noise and no
measurement sample. It asks two things about the 12% bias of the fifth-degree
jet on the half-cycle benchmark:

1. Can windows displaced along the observer's worldline remove it at a smaller
   cost than the sixth-degree jet (radius x5-6, events x23-37)?
2. How large is an explicit first-order bias bound under a stated regularity
   assumption, for plane waves and in general?

Setting. The window centred at proper-time offset ``shift`` on the observer's
geodesic (x = y = 0 of the Brinkmann chart) is isometric to the standard
window of the wave with profile a(u + shift). The Christoffel symbols vanish on
that axis, so the coordinate frame is parallel along it and the normal frames
of all windows are parallel transports of one another (the "generator frame").
Each window reuses the landmark Brinkmann coordinates, clock mask and clock
covariance of the benchmark run; its clocks are exact interval volumes, and it
is fitted by the landmark harness's continuum fit started at the truth.

Parts:
1. Window fits for the half-cycle wave (wavelength 2 diamond heights) and for
   longer sinusoids, and one fit of the quadratic profile for validation.
2. The linear alias map L of the 100 packed second-derivative components onto
   the ten Weyl coefficients of the fifth-degree fit: Weyl rows of J5^+ Jd for
   the whitened Jacobian at the window-0 inferred coordinates.
3. Corrections c0 - L K D2[c] (+ (L K)^2 D4[c]), K = k k (x) (k the wave
   covector: oracle, or the principal null direction of the window-0 fit),
   D2 and D4 finite-difference stencils along the worldline; and a
   minimum-norm variant that uses only the measured time-time slice.
4. Bias bounds: worst case over propagation direction and polarization of the
   first-order plane-wave bias, and the operator norm of L under a stated
   regularity norm on the second derivatives.
5. Which displacement directions determine the part of the second
   derivatives that L sees, when no wave type is assumed.
6. Robustness: side windows from the other benchmark layout (each window with
   its own landmarks, clocks and covariance), and side windows whose frames
   are misaligned with the centre window by a fixed angle.

Not established here: frame alignment between windows is oracle (each
window's gauge frame is mapped to the generator frame with the truth; part 6
measures how much alignment the correction needs, not how to obtain it), the
windows are independent sprinkles, the radii are nominal Gaussian quantities
under the run's clock covariance, and nothing here is a finite-N result.
"""

import argparse
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_N = 2**24
TARGET_RADIUS = 0.10
MAIN_SHIFTS = (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0)
LONG_SHIFTS = (-1.0, -0.5, 0.0, 0.5, 1.0)
LONG_WAVELENGTHS = (3.0, 4.0, 6.0)
# Stencils on the worldline: name -> (offsets in units of h, D2 weights * h^2,
# D4 weights * h^4 or None).
STENCILS = {
    "three_point": ((-1, 0, 1), (1, -2, 1), None),
    "five_point": (
        (-2, -1, 0, 1, 2),
        tuple(np.array([-1, 16, -30, 16, -1]) / 12),
        (1, -4, 6, -4, 1),
    ),
}


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------


def decimal_sin_cos(x):
    """sin and cos of a Decimal argument by the Taylor series."""
    from profile_plane_wave import PRECISION

    with localcontext() as context:
        context.prec = PRECISION + 20
        x = Decimal(x)
        sine = cosine = Decimal(0)
        term, k = Decimal(1), 0
        while k < 8 or abs(term) > Decimal(10) ** -(PRECISION + 15):
            if k % 4 == 0:
                cosine += term
            elif k % 4 == 1:
                sine += term
            elif k % 4 == 2:
                cosine -= term
            else:
                sine -= term
            k += 1
            term = term * x / k
        return sine, cosine


def sinusoid_taylor(amplitude, wavelength, shift):
    """Taylor coefficients at u = 0 of a(u + shift), a(u) = A cos(w (u - 1/2)).

    w = 2 pi / wavelength, the wavelength in diamond heights; wavelength 2 is
    the crest-centred half cycle A sin(pi u) of the benchmark.
    """
    from profile_plane_wave import ORDER, PRECISION, decimal_pi

    with localcontext() as context:
        context.prec = PRECISION + 20
        a = Decimal(amplitude)
        omega = 2 * decimal_pi() / Decimal(wavelength)
        sine, cosine = decimal_sin_cos(omega * (Decimal(shift) - Decimal(1) / 2))
        # A cos(w u + theta) = A [cos(theta) cos(w u) - sin(theta) sin(w u)].
        taylor, term = [], Decimal(1)
        for k in range(ORDER):
            sign = -1 if (k // 2) % 2 else 1
            taylor.append(sign * a * term * (cosine if k % 2 == 0 else -sine))
            term = term * omega / (k + 1)
    with localcontext() as context:
        context.prec = PRECISION
        return [+value for value in taylor]


def build_window(kind, amplitude, wavelength, shift):
    from profile_plane_wave import ProfileWave

    if kind == "quadratic":
        from intrinsic_weyl_second_continuum import build

        if shift != 0:
            raise ValueError("The quadratic profile is used at the centre only")
        return build("quadratic", amplitude)
    taylor = sinusoid_taylor(amplitude, wavelength, shift)
    return ProfileWave(taylor, abs(float(amplitude)), f"sinusoid {wavelength} {shift}")


def centre_amplitude(kind, amplitude, wavelength, shift):
    if kind == "quadratic":
        return float(amplitude)
    return float(amplitude) * float(np.cos(2 * np.pi * shift / wavelength))


# ---------------------------------------------------------------------------
# Frames, packed second derivatives, plane-wave operator
# ---------------------------------------------------------------------------


def gauge_frame(points):
    """Proper rotation F with gauge(points)[:, 1:] = points[:, 1:] @ F."""
    first = points[1, 1:] / np.linalg.norm(points[1, 1:])
    second = points[2, 1:] - np.dot(points[2, 1:], first) * first
    second /= np.linalg.norm(second)
    return np.column_stack((first, second, np.cross(first, second)))


def spatial(frame):
    full = np.eye(4)
    full[1:, 1:] = frame
    return full


def weyl_rotation(frame):
    """10x10 map of Weyl coefficients to the coordinates x' = x @ frame."""
    from intrinsic_weyl import WEYL

    full = spatial(frame)
    rotated = np.einsum(
        "ia,jb,kc,ld,nijkl->nabcd", full, full, full, full, WEYL, optimize=True
    )
    matrix = np.linalg.lstsq(
        WEYL.reshape(10, -1).T, rotated.reshape(10, -1).T, rcond=None
    )[0]
    if np.max(abs(np.einsum("mn,mabcd->nabcd", matrix, WEYL) - rotated)) > 1e-12:
        raise ArithmeticError("Rotated Weyl tensor left the Weyl basis")
    return matrix


def unpack_second(packed):
    from intrinsic_weyl import SYMMETRIC_PAIRS

    packed = np.asarray(packed).reshape(10, 10)
    full = np.zeros((4, 4, 10))
    for s, (e, f) in enumerate(SYMMETRIC_PAIRS):
        full[e, f] = full[f, e] = packed[s]
    return full


def pack_second(full):
    from intrinsic_weyl import SYMMETRIC_PAIRS

    return np.array([full[e, f] for e, f in SYMMETRIC_PAIRS]).ravel()


def second_rotation(frame):
    """100x100 map of packed (derivative pair, Weyl) components to x @ frame."""
    full = spatial(frame)
    weyl = weyl_rotation(frame)
    matrix = np.zeros((100, 100))
    for column in range(100):
        unit = np.zeros(100)
        unit[column] = 1
        tensor = np.einsum("ea,fb,efk->abk", full, full, unpack_second(unit))
        matrix[:, column] = pack_second(tensor @ weyl.T)
    return matrix


def plane_wave_operator(k):
    """100x10 map C -> packed k_e k_f C, with k the wave covector."""
    from intrinsic_weyl import SYMMETRIC_PAIRS

    operator = np.zeros((10, 10, 10))
    for s, (e, f) in enumerate(SYMMETRIC_PAIRS):
        operator[s] = k[e] * k[f] * np.eye(10)
    return operator.reshape(100, 10)


def second_metric():
    """Regularity norm on the 42 vacuum second derivatives.

    |D|^2 = sum over ordered (e, f) in the observer frame of the full-tidal
    norm of D_ef; a plane wave C f(phi) with covector (1, n) has
    |D| = 2 |f''| |C|.
    """
    from intrinsic_weyl import SYMMETRIC_PAIRS, second_derivative_basis
    from intrinsic_weyl_joint import TIDAL_GRAM

    blocks = np.zeros((100, 100))
    for s, (e, f) in enumerate(SYMMETRIC_PAIRS):
        blocks[10 * s : 10 * s + 10, 10 * s : 10 * s + 10] = (
            1 if e == f else 2
        ) * TIDAL_GRAM
    basis = second_derivative_basis()
    return basis.T @ blocks @ basis


def inverse_root(matrix):
    values, vectors = np.linalg.eigh(matrix)
    if values[0] <= 0:
        raise ArithmeticError("Metric is not positive definite")
    return vectors @ np.diag(values**-0.5) @ vectors.T


def principal_null_direction(coefficients, grid=4096):
    """Unit spatial n minimizing |C_abcd l^d| for l = (1, n); type N gives 0."""
    from intrinsic_weyl import WEYL
    from scipy.optimize import minimize

    tensor = np.einsum("k,kabcd->abcd", coefficients, WEYL)
    tensor = tensor / np.linalg.norm(tensor)

    def residual(vector):
        n = vector / np.linalg.norm(vector)
        return float(np.sum(np.einsum("abcd,d->abc", tensor, np.r_[1.0, n]) ** 2))

    golden = np.pi * (3 - np.sqrt(5))
    i = np.arange(grid) + 0.5
    z = 1 - 2 * i / grid
    r = np.sqrt(1 - z * z)
    candidates = np.column_stack((r * np.cos(golden * i), r * np.sin(golden * i), z))
    lights = np.column_stack((np.ones(grid), candidates))
    values = (np.einsum("abcd,gd->gabc", tensor, lights) ** 2).sum(axis=(1, 2, 3))
    best = candidates[int(np.argmin(values))]
    result = minimize(
        residual,
        best,
        method="Nelder-Mead",
        options={"xatol": 1e-13, "fatol": 1e-24, "maxiter": 8000},
    )
    n = result.x / np.linalg.norm(result.x)
    return n, residual(n)


def covector(n):
    """Wave covector (1, -n) for propagation direction n: phase t - n.x."""
    return np.r_[1.0, -np.asarray(n)]


def axis_rotation(axis, angle):
    axis = np.asarray(axis) / np.linalg.norm(axis)
    skew = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    return np.eye(3) + np.sin(angle) * skew + (1 - np.cos(angle)) * skew @ skew


def rotation_to(target, source=(0.0, 0.0, -1.0)):
    """Active proper rotation Q with Q @ source = target."""
    source, target = np.asarray(source), np.asarray(target)
    axis = np.cross(source, target)
    sine, cosine = np.linalg.norm(axis), float(source @ target)
    if sine < 1e-12:
        return np.eye(3) if cosine > 0 else axis_rotation([1.0, 0, 0], np.pi)
    return axis_rotation(axis, np.arctan2(sine, cosine))


def unit_truth():
    """Weyl coefficients of unit profile amplitude in the generator frame."""
    from intrinsic_weyl_joint_audit import aligned_tensor_error

    points = np.array(
        [
            [-0.5, 0, 0, 0],
            [0.1, 0.2, 0.03, -0.05],
            [0.0, -0.04, 0.17, 0.02],
            [0.05, 0.01, -0.03, 0.21],
            [0.5, 0, 0, 0],
        ]
    )
    result = aligned_tensor_error(
        {"coefficients": np.zeros(10), "coordinates_inferred_from_counts": points},
        points,
        1.0,
    )
    return np.asarray(result["aligned_truth_coefficients"])


def tidal_norm(coefficients):
    from intrinsic_weyl_joint import TIDAL_ROOT

    return float(np.linalg.norm(TIDAL_ROOT @ coefficients))


def radius_95(covariance):
    from intrinsic_weyl_joint import TIDAL_ROOT
    from scipy.stats import chi2

    values = np.linalg.eigvalsh(TIDAL_ROOT @ covariance @ TIDAL_ROOT.T)
    return float(np.sqrt(chi2.ppf(0.95, 10) * values[-1]))


# ---------------------------------------------------------------------------
# Window fits
# ---------------------------------------------------------------------------


def clock_cache(cache, archive, case):
    """Run clocks of a benchmark case; computed and stored when absent."""
    from intrinsic_weyl_joint import read_order_clocks

    path = Path(cache) / f"run_clocks_case{case}_n{REFERENCE_N}.npz"
    if not path.exists():
        counts, mask, covariance, _, _ = read_order_clocks(
            Path(archive) / f"case_{case}_n{REFERENCE_N}.npz"
        )
        partial = path.with_name(f"{path.stem}.{os.getpid()}.partial.npz")
        np.savez(partial, counts=counts, mask=mask, cov=covariance)
        os.replace(partial, path)
    return path


def check_cache(path, archive, case):
    """The cached clock mask must be the archive's order-only clock mask."""
    from intrinsic_weyl_joint import selection

    with np.load(Path(archive) / f"case_{case}_n{REFERENCE_N}.npz") as saved:
        order = saved["selection_order"]
        ids = saved["selected_ids"]
        design = (
            str(saved["selection_design"].item())
            if "selection_design" in saved
            else "two"
        )
    expected_ids, expected_mask = selection(order, design=design)
    with np.load(path) as saved:
        mask, counts = saved["mask"], saved["counts"]
    return {
        "selected_ids_equal": bool(np.array_equal(ids, expected_ids)),
        "clock_mask_equal": bool(np.array_equal(mask, expected_mask)),
        "tip_count_equals_n": bool(counts[0, -1] == REFERENCE_N),
    }


def window_job(job):
    """One window: exact volumes, normal coordinates, continuum fits."""
    from intrinsic_weyl_joint import gauge
    from intrinsic_weyl_landmark_design import continuum_fit

    begin = time.time()
    with np.load(job["truth"]) as saved:
        points = saved["points"]
    with np.load(job["clocks"]) as saved:
        mask, covariance = saved["mask"], saved["cov"]
    wave = build_window(job["kind"], job["amplitude"], job["wavelength"], job["shift"])
    normal = wave.normal_coordinates(points)
    volume = wave.volumes(points)
    upper = np.triu(mask, 1)
    upper[0, -1] = False
    pairs = np.column_stack(np.nonzero(upper))
    if not np.all(volume[pairs[:, 0], pairs[:, 1]] > 0):
        raise ArithmeticError("A clock has no exact wave volume")
    start = gauge(normal)
    continuum = REFERENCE_N * volume / volume[0, -1]
    result = {
        key: job[key] for key in ("case", "kind", "wavelength", "shift", "orders")
    }
    result["frame"] = gauge_frame(normal).tolist()
    result["start_coordinates"] = start.tolist()
    result["centre_amplitude"] = centre_amplitude(
        job["kind"], job["amplitude"], job["wavelength"], job["shift"]
    )
    result["setup_seconds"] = round(time.time() - begin, 1)
    for order in job["orders"]:
        begin = time.time()
        fit = continuum_fit(continuum, mask, covariance, order, start)
        result[f"order_{order}"] = {
            "coefficients": fit["coefficients"],
            "covariance": fit["nominal_component_covariance"],
            "coordinates": fit["coordinates_inferred_from_counts"],
            "optimizer_success": fit["optimizer_success"],
            "function_evaluations": fit["function_evaluations"],
            "curved_cost": fit["curved_cost"],
            "projected_curvature_rank": fit["projected_curvature_rank"],
            "nuisance_rank": fit["nuisance_rank"],
            "nuisance_parameters": fit["nuisance_parameters"],
            "seconds": round(time.time() - begin, 1),
        }
    return result


# ---------------------------------------------------------------------------
# Alias map
# ---------------------------------------------------------------------------


def whitening(covariance):
    """The whitening of fit_joint and continuum_fit."""
    scales = np.sqrt(np.diag(covariance))
    values, vectors = np.linalg.eigh(covariance / scales[:, None] / scales)
    keep = values > values[-1] * 1e-8
    return (vectors[:, keep] / np.sqrt(values[keep])).T / scales


def alias_map(points, mask, whiten):
    """Weyl rows of J5^+ Jd at the given gauge-fixed coordinates, zero curvature.

    J5: coordinates, ten Weyl and 24 gradient columns of the fifth-degree jet
    with curvature-square terms; Jd: the 100 packed second-derivative
    components (symmetric fill), whose span-of-basis restriction is the
    sixth-degree jet's 42 columns. Both whitened.
    """
    from intrinsic_weyl import PAIRS, SYMMETRIC_PAIRS, second_derivative_basis
    from intrinsic_weyl_joint import BIVECTOR, model

    n = len(points)
    upper = np.triu(mask, 1)
    upper[0, -1] = False
    pairs = np.column_stack(np.nonzero(upper))
    _, dl, dr, dc = model(points, np.zeros(76), pairs, square_correction=True)
    free = np.ones((n, 4), dtype=bool)
    free[[0, -1]] = False
    free[1, 2:] = False
    free[2, 3] = False
    index = np.full((n, 4), -1)
    index[free] = np.arange(free.sum())
    dimension = int(free.sum())
    geometry = np.zeros((len(pairs), dimension))
    for side, derivative in enumerate((dl, dr)):
        for component in range(4):
            column = index[pairs[:, side], component]
            valid = column >= 0
            geometry[np.flatnonzero(valid), column[valid]] = derivative[
                valid, component
            ]
    left, right = points[pairs[:, 0]], points[pairs[:, 1]]
    wedge = np.column_stack(
        [left[:, a] * right[:, b] - left[:, b] * right[:, a] for a, b in PAIRS]
    )
    quartic = np.einsum("pi,kij,pj->pk", wedge, BIVECTOR, wedge, optimize=True)
    middle, step = (left + right) / 2, right - left
    weight = (
        np.einsum("pe,pf->pef", middle, middle) / 20
        + np.einsum("pe,pf->pef", step, step) / 240
    )
    packed = np.zeros((len(pairs), 10, 10))
    for s, (e, f) in enumerate(SYMMETRIC_PAIRS):
        packed[:, s, :] = -quartic * ((1 if e == f else 2) * weight[:, e, f])[:, None]
    packed = packed.reshape(len(pairs), 100)
    consistency = float(
        np.max(abs(packed @ second_derivative_basis() - dc[:, 34:]))
        / np.max(abs(dc[:, 34:]))
    )
    j5 = whiten @ np.column_stack((geometry, dc[:, :34]))
    inverse = np.linalg.pinv(j5)
    alias = (inverse @ (whiten @ packed))[dimension : dimension + 10]
    influence = inverse[dimension : dimension + 10]
    return alias, consistency, influence @ influence.T


# ---------------------------------------------------------------------------
# Corrections
# ---------------------------------------------------------------------------


def combine(windows, weights):
    """Linear combination sum_j W_j c_j and its covariance (independent)."""
    estimate = sum(w @ c for w, (c, _) in zip(weights, windows, strict=True))
    covariance = sum(w @ s @ w.T for w, (_, s) in zip(weights, windows, strict=True))
    return estimate, covariance


def correction_specs(available, alias_n, oracle_k, estimated_k, minimum_norm):
    """(name, h, operator, shifts, weights) of every correction the windows allow.

    The operator maps the worldline second derivative of the Weyl coefficients
    to the first-order bias: L K(k) for a plane wave with covector k (oracle,
    or the principal null direction of the window-0 fit), or the minimum-norm
    map that assumes no wave type. "first": c0 - B D2[c]; "second" adds
    B^2 D4[c], the next term of (I + B d^2/dt^2)^-1.
    """
    specs = [("single_window", 0.0, "none", [0.0], [np.eye(10)])]
    operators = (
        ("plane_wave_oracle_k", alias_n @ plane_wave_operator(oracle_k)),
        ("plane_wave_estimated_k", alias_n @ plane_wave_operator(estimated_k)),
        ("minimum_norm_time_time", minimum_norm),
    )
    for h in (0.25, 0.5, 1.0):
        for name, (offsets, d2, d4) in STENCILS.items():
            chosen = [round(o * h, 10) for o in offsets]
            if not all(s in available for s in chosen):
                continue
            for label, operator in operators:
                first = [
                    (s == 0) * np.eye(10) - operator * w / h**2
                    for s, w in zip(chosen, d2, strict=True)
                ]
                specs.append((f"{name}_first", h, label, chosen, first))
                if d4 is not None:
                    second = [
                        w + operator @ operator * v / h**4
                        for w, v in zip(first, d4, strict=True)
                    ]
                    specs.append((f"{name}_second", h, label, chosen, second))
    return specs


def evaluate_specs(window_set, specs, truth):
    """Bias, nominal radius and event cost of each correction.

    window_set: {shift: (coefficients, covariance)} in the generator frame.
    """
    single_radius = radius_95(window_set[0.0][1])
    truth_norm = tidal_norm(truth)
    rows = []
    for name, h, operator, chosen, weights in specs:
        estimate, covariance = combine([window_set[s] for s in chosen], weights)
        radius = radius_95(covariance)
        rows.append(
            {
                "correction": name,
                "h": h,
                "operator": operator,
                "windows": len(chosen),
                "shifts": list(chosen),
                "relative_bias": tidal_norm(estimate - truth) / truth_norm,
                "relative_radius_95": radius / truth_norm,
                "radius_ratio_to_single": radius / single_radius,
                "events_ratio_to_single": len(chosen) * (radius / single_radius) ** 2,
                "events_for_radius_10pct": len(chosen)
                * REFERENCE_N
                * (radius / truth_norm / TARGET_RADIUS) ** 2,
            }
        )
    return rows


def misalignment(window_set, specs, truth, angles=(1, 3, 10), draws=64):
    """Worst bias when every side window's frame is off by a fixed angle.

    The axis is drawn at random per window and draw; the centre window's frame
    is the reference. Only the h = 1/2 corrections are probed.
    """
    rng = np.random.default_rng(2026112104)
    rows = []
    for name, h, operator, chosen, weights in specs:
        if h != 0.5 or name not in ("three_point_first", "five_point_second"):
            continue
        for angle in angles:
            worst = 0.0
            for _ in range(draws):
                estimate = np.zeros(10)
                for s, w in zip(chosen, weights, strict=True):
                    c = window_set[s][0]
                    if s != 0:
                        turn = axis_rotation(rng.normal(size=3), np.radians(angle))
                        c = weyl_rotation(turn.T) @ c
                    estimate = estimate + w @ c
                worst = max(worst, tidal_norm(estimate - truth) / tidal_norm(truth))
            rows.append(
                {
                    "correction": name,
                    "h": h,
                    "operator": operator,
                    "angle_deg": angle,
                    "worst_relative_bias": worst,
                }
            )
    return rows


# ---------------------------------------------------------------------------
# Bounds and displacement directions
# ---------------------------------------------------------------------------


def fibonacci_sphere(count):
    golden = np.pi * (3 - np.sqrt(5))
    i = np.arange(count) + 0.5
    z = 1 - 2 * i / count
    r = np.sqrt(1 - z * z)
    return np.column_stack((r * np.cos(golden * i), r * np.sin(golden * i), z))


def plane_wave_scan(alias_n, unit, directions=2000):
    """First-order relative bias per unit |f''/f| over direction and polarization.

    For propagation direction n the type-N tensors are the reference tensor
    actively rotated so that its principal null direction becomes n, and its
    45-degree rotation about n; the worst elliptical polarization is the
    largest singular value on that tidal-orthonormal pair.
    """
    from intrinsic_weyl_joint import TIDAL_ROOT

    values, cosines, projections, worst_pnd_error = [], [], [], 0.0
    for n in fibonacci_sphere(directions):
        q = rotation_to(n)
        plus = weyl_rotation(q.T) @ unit
        cross = weyl_rotation(axis_rotation(n, np.pi / 4).T) @ plus
        tensors = np.column_stack((plus, cross))
        pair = TIDAL_ROOT @ tensors
        _, upper = np.linalg.qr(pair)
        response = TIDAL_ROOT @ alias_n @ plane_wave_operator(covector(n))
        biases = response @ tensors
        values.append(
            float(np.linalg.svd(biases @ np.linalg.inv(upper), compute_uv=False)[0])
        )
        # Is the first-order bias a rescaling of the tensor itself?
        for j in range(2):
            projections.append(
                float(biases[:, j] @ pair[:, j] / (pair[:, j] @ pair[:, j]))
            )
            cosines.append(
                float(
                    biases[:, j]
                    @ pair[:, j]
                    / np.linalg.norm(biases[:, j])
                    / np.linalg.norm(pair[:, j])
                )
            )
        if len(values) <= 16:
            estimated, _ = principal_null_direction(plus, grid=1024)
            worst_pnd_error = max(worst_pnd_error, float(np.linalg.norm(estimated - n)))
    values = np.array(values)
    benchmark = TIDAL_ROOT @ alias_n @ plane_wave_operator(covector([0, 0, -1.0]))
    return {
        "directions": directions,
        "worst": float(values.max()),
        "best": float(values.min()),
        "median": float(np.median(values)),
        "quantiles_10_50_90": np.quantile(values, [0.1, 0.5, 0.9]).tolist(),
        "benchmark_direction_benchmark_polarization": float(
            np.linalg.norm(benchmark @ unit) / tidal_norm(unit)
        ),
        "linear_polarizations_bias_to_tensor_cosine_min_max": [
            min(cosines),
            max(cosines),
        ],
        "linear_polarizations_bias_projection_on_tensor_min_max": [
            min(projections),
            max(projections),
        ],
        "rotation_check_pnd_max_error": worst_pnd_error,
    }


def displacement_functional(n):
    """10x42: second derivatives -> n^e n^f D_ef (second difference along n)."""
    from intrinsic_weyl import SYMMETRIC_PAIRS, second_derivative_basis

    weights = np.array([(1 if e == f else 2) * n[e] * n[f] for e, f in SYMMETRIC_PAIRS])
    packed = second_derivative_basis().reshape(10, 10, 42)
    return np.einsum("s,skj->kj", weights, packed)


def displacement_study(alias42_n, metric_root_inverse):
    """Greedy choice of displacement directions until L's first-order bias is
    determined, with no wave type assumed. Residual: the worst-case relative
    part of the first-order bias, in the regularity norm, that the directions
    leave undetermined."""
    from intrinsic_weyl_joint import TIDAL_ROOT

    lhs = TIDAL_ROOT @ alias42_n @ metric_root_inverse
    scale = np.linalg.svd(lhs, compute_uv=False)[0]
    axes = np.eye(4)
    names = {"t": axes[0], "x": axes[1], "y": axes[2], "z": axes[3]}
    for i, a in enumerate("txyz"):
        for j, b in enumerate("txyz"):
            if i < j:
                names[f"{a}+{b}"] = (axes[i] + axes[j]) / np.sqrt(2)
                names[f"{a}-{b}"] = (axes[i] - axes[j]) / np.sqrt(2)

    def residual(chosen):
        stacked = np.vstack([displacement_functional(names[c]) for c in chosen])
        whitened = stacked @ metric_root_inverse
        projector = np.eye(42) - np.linalg.pinv(whitened) @ whitened
        return (
            float(np.linalg.svd(lhs @ projector, compute_uv=False)[0] / scale),
            int(np.linalg.matrix_rank(whitened, tol=1e-9)),
        )

    chosen, path = ["t"], []
    value, rank = residual(chosen)
    path.append({"directions": list(chosen), "residual": value, "rank": rank})
    while value > 1e-9 and len(chosen) < len(names):
        trials = {c: residual(chosen + [c])[0] for c in names if c not in chosen}
        best = min(trials, key=trials.get)
        chosen.append(best)
        value, rank = residual(chosen)
        path.append({"directions": list(chosen), "residual": value, "rank": rank})
    return {
        "greedy_path": path,
        "windows_needed": 1 + 2 * len(chosen) if value <= 1e-9 else None,
    }


# ---------------------------------------------------------------------------
# Selfcheck
# ---------------------------------------------------------------------------


def selfcheck(amplitude):
    from intrinsic_weyl import second_derivative_basis
    from profile_plane_wave import half_cycle

    report = {}
    # 1. The shifted-sinusoid generator reproduces the benchmark half cycle.
    reference = half_cycle(amplitude)
    ours = build_window("sinusoid", amplitude, 2.0, 0.0)
    report["half_cycle_series_max_difference"] = float(
        np.max(abs(ours.coefficients - reference.coefficients))
    )
    assert report["half_cycle_series_max_difference"] < 1e-15
    # 2. Taylor coefficients against the closed form of the shifted profile.
    worst = 0.0
    for wavelength, shift in ((2.0, 0.5), (3.0, -1.0), (6.0, 0.25)):
        taylor = np.array([float(v) for v in sinusoid_taylor(1, wavelength, shift)])
        for u in (0.0, 0.3, 0.77, 1.0):
            exact = np.cos(2 * np.pi * (u + shift - 0.5) / wavelength)
            worst = max(worst, abs(np.polynomial.polynomial.polyval(u, taylor) - exact))
    report["shifted_profile_max_error"] = float(worst)
    assert worst < 1e-13
    # 3. Window isometry: actions of the shifted window equal the base wave's
    # actions between the shifted events (u + 1/4, v - 1/4).
    rng = np.random.default_rng(2026112103)
    left = np.column_stack(
        (
            rng.uniform(0, 0.5, 64),
            rng.uniform(-0.3, 0.3, 64),
            rng.normal(0, 0.1, (64, 2)),
        )
    )
    right = left.copy()
    right[:, 0] = left[:, 0] + rng.uniform(0.05, 0.25, 64)
    right[:, 2:] += rng.normal(0, 0.05, (64, 2))
    shifted = build_window("sinusoid", amplitude, 2.0, 0.25)
    move = np.array([0.25, -0.25, 0, 0])
    mine = shifted._pair(left, right)[0]
    theirs = reference._pair(left + move, right + move)[0]
    report["window_isometry_action_max_difference"] = float(np.max(abs(mine - theirs)))
    assert report["window_isometry_action_max_difference"] < 1e-12
    # 4. Rotations: tidal norm preserved, inverse, 42-space invariant.
    frame = rotation_to(np.array([0.3, -0.5, 0.81]) / np.linalg.norm([0.3, -0.5, 0.81]))
    forward, backward = weyl_rotation(frame), weyl_rotation(frame.T)
    report["weyl_rotation_inverse_error"] = float(
        np.max(abs(forward @ backward - np.eye(10)))
    )
    vector = rng.normal(size=10)
    report["weyl_rotation_norm_change"] = abs(
        tidal_norm(forward @ vector) - tidal_norm(vector)
    )
    basis = second_derivative_basis()
    rotated = second_rotation(frame) @ basis
    report["second_rotation_leaves_vacuum_space"] = float(
        np.max(abs(rotated - basis @ np.linalg.lstsq(basis, rotated, rcond=None)[0]))
    )
    assert report["weyl_rotation_inverse_error"] < 1e-12
    assert report["weyl_rotation_norm_change"] < 1e-12
    assert report["second_rotation_leaves_vacuum_space"] < 1e-12
    # 5. The benchmark truth: unit amplitude has full-tidal norm 2, principal
    # null direction (0, 0, -1), and k k C with k = (1, 0, 0, 1) is a vacuum
    # second derivative (differential Bianchi identity) while k = (1, 0, 0, -1)
    # is not.
    unit = unit_truth()
    report["unit_truth_tidal_norm"] = tidal_norm(unit)
    assert abs(report["unit_truth_tidal_norm"] - 2) < 1e-12
    n, value = principal_null_direction(unit)
    report["unit_truth_pnd"] = n.tolist()
    report["unit_truth_pnd_residual"] = value
    assert np.linalg.norm(n - [0, 0, -1]) < 1e-6
    projector = basis @ np.linalg.pinv(basis)
    for label, k in (("matched", [1.0, 0, 0, 1]), ("mismatched", [1.0, 0, 0, -1])):
        packed = plane_wave_operator(np.array(k)) @ unit
        report[f"plane_wave_bianchi_residual_{label}"] = float(
            np.linalg.norm(packed - projector @ packed) / np.linalg.norm(packed)
        )
    assert report["plane_wave_bianchi_residual_matched"] < 1e-12
    assert report["plane_wave_bianchi_residual_mismatched"] > 0.1
    # 6. Regularity metric: a plane wave has |D| = 2 |C|.
    metric = second_metric()
    coefficients = np.linalg.lstsq(
        basis, plane_wave_operator(np.array([1.0, 0, 0, 1])) @ unit, rcond=None
    )[0]
    report["plane_wave_regularity_norm_over_curvature"] = float(
        np.sqrt(coefficients @ metric @ coefficients) / tidal_norm(unit)
    )
    assert abs(report["plane_wave_regularity_norm_over_curvature"] - 2) < 1e-12
    return report


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def jobs_for(case, archive, cache, amplitude, long_wavelengths):
    clocks = str(Path(cache) / f"run_clocks_case{case}_n{REFERENCE_N}.npz")
    truth = str(Path(archive) / f"truth_{case}.npz")
    common = {
        "case": case,
        "amplitude": str(amplitude),
        "clocks": clocks,
        "truth": truth,
    }
    jobs = [
        dict(
            common,
            kind="sinusoid",
            wavelength=2.0,
            shift=s,
            orders=(5, 6) if s == 0 else (5,),
        )
        for s in MAIN_SHIFTS
    ]
    jobs.append(dict(common, kind="quadratic", wavelength=None, shift=0.0, orders=(5,)))
    for wavelength in long_wavelengths:
        jobs += [
            dict(
                common,
                kind="sinusoid",
                wavelength=wavelength,
                shift=s,
                orders=(5, 6) if s == 0 else (5,),
            )
            for s in LONG_SHIFTS
        ]
    return jobs


def analyse_case(case, windows, cache, amplitude):
    """Alias map, validation, corrections and bounds for one benchmark case."""
    from intrinsic_weyl import second_derivative_basis
    from intrinsic_weyl_joint import TIDAL_ROOT

    with np.load(Path(cache) / f"run_clocks_case{case}_n{REFERENCE_N}.npz") as saved:
        mask, covariance = saved["mask"], saved["cov"]
    whiten = whitening(covariance)
    unit = unit_truth()
    centre = next(
        w
        for w in windows
        if w["kind"] == "sinusoid" and w["wavelength"] == 2.0 and w["shift"] == 0.0
    )
    frame0 = np.asarray(centre["frame"])
    inferred = np.asarray(centre["order_5"]["coordinates"])
    alias_g, consistency, influence_cov = alias_map(inferred, mask, whiten)
    alias_truth, _, _ = alias_map(np.asarray(centre["start_coordinates"]), mask, whiten)
    to_normal = weyl_rotation(frame0.T)
    alias_n = to_normal @ alias_g @ second_rotation(frame0)
    basis = second_derivative_basis()
    alias42_n = alias_n @ basis
    metric = second_metric()
    root_inverse = inverse_root(metric)
    oracle_k = np.array([1.0, 0, 0, 1])
    result = {
        "case": case,
        "landmarks": len(inferred) - 2,
        "clocks": int(np.triu(mask, 1).sum()),
        "alias": {
            "evaluated_at": "window-0 fifth-degree inferred coordinates",
            "packed_columns_vs_sixth_degree_model_max_relative": consistency,
            "inferred_vs_truth_coordinates_relative_difference": float(
                np.linalg.norm(alias_g - alias_truth) / np.linalg.norm(alias_truth)
            ),
            "linearized_covariance_vs_fit_relative_difference": float(
                np.linalg.norm(
                    influence_cov - np.asarray(centre["order_5"]["covariance"])
                )
                / np.linalg.norm(influence_cov)
            ),
        },
    }

    def in_normal(window, order=5):
        rotation = weyl_rotation(np.asarray(window["frame"]).T)
        fit = window[f"order_{order}"]
        return (
            rotation @ np.asarray(fit["coefficients"]),
            rotation @ np.asarray(fit["covariance"]) @ rotation.T,
        )

    # Validation: the quadratic profile's continuum bias against L K D.
    quadratic = next(w for w in windows if w["kind"] == "quadratic")
    estimate, _ = in_normal(quadratic)
    truth = amplitude * unit
    bias = estimate - truth
    second = -(np.pi**2) * amplitude * unit
    prediction = alias_n @ plane_wave_operator(oracle_k) @ second
    result["alias"]["quadratic_validation"] = {
        "continuum_relative_bias": tidal_norm(bias) / tidal_norm(truth),
        "first_order_prediction": tidal_norm(prediction) / tidal_norm(truth),
        "difference_relative_to_truth": tidal_norm(bias - prediction)
        / tidal_norm(truth),
        "cosine": float(
            (TIDAL_ROOT @ bias)
            @ (TIDAL_ROOT @ prediction)
            / tidal_norm(bias)
            / tidal_norm(prediction)
        ),
    }
    # Minimum-norm second derivative with a given time-time slice.
    tt = displacement_functional(np.array([1.0, 0, 0, 0]))
    metric_inverse = np.linalg.inv(metric)
    gram = tt @ metric_inverse @ tt.T
    result["time_time_slice_rank"] = int(np.linalg.matrix_rank(gram))
    minimum_norm = alias42_n @ metric_inverse @ tt.T @ np.linalg.pinv(gram)
    rows, sixth, per_window, pnd, turned = [], [], [], [], []
    context = {"alias_n": alias_n, "minimum_norm": minimum_norm, "truth": truth}
    context["window_sets"], context["estimated_k"] = {}, {}
    for wavelength in sorted(
        {w["wavelength"] for w in windows if w["kind"] == "sinusoid"}
    ):
        chosen = [
            w
            for w in windows
            if w["kind"] == "sinusoid" and w["wavelength"] == wavelength
        ]
        window_set = {w["shift"]: in_normal(w) for w in chosen}
        for w in chosen:
            c, s = window_set[w["shift"]]
            per_window.append(
                {
                    "wavelength": wavelength,
                    "shift": w["shift"],
                    "centre_amplitude_over_benchmark": w["centre_amplitude"]
                    / amplitude,
                    "relative_bias_order_5": tidal_norm(
                        c - w["centre_amplitude"] * unit
                    )
                    / tidal_norm(truth),
                    "relative_radius_95_order_5": radius_95(s) / tidal_norm(truth),
                    "optimizer_success": w["order_5"]["optimizer_success"],
                    "function_evaluations": w["order_5"]["function_evaluations"],
                    "seconds": w["order_5"]["seconds"],
                }
            )
        centre_window = next(w for w in chosen if w["shift"] == 0.0)
        estimate6, covariance6 = in_normal(centre_window, 6)
        radius5 = radius_95(window_set[0.0][1])
        radius6 = radius_95(covariance6)
        sixth.append(
            {
                "wavelength": wavelength,
                "relative_bias": tidal_norm(estimate6 - truth) / tidal_norm(truth),
                "relative_radius_95": radius6 / tidal_norm(truth),
                "radius_ratio_to_single": radius6 / radius5,
                "events_ratio_to_single": (radius6 / radius5) ** 2,
                "optimizer_success": centre_window["order_6"]["optimizer_success"],
            }
        )
        n, value = principal_null_direction(window_set[0.0][0])
        angle = float(np.degrees(np.arccos(np.clip(-n[2], -1, 1))))
        pnd.append(
            {
                "wavelength": wavelength,
                "direction": n.tolist(),
                "residual": value,
                "angle_to_truth_deg": angle,
            }
        )
        first_order = (
            alias_n
            @ plane_wave_operator(oracle_k)
            @ (-((2 * np.pi / wavelength) ** 2) * truth)
        )
        specs = correction_specs(
            sorted(window_set), alias_n, oracle_k, covector(n), minimum_norm
        )
        for item in evaluate_specs(window_set, specs, truth):
            item["wavelength"] = wavelength
            item["first_order_prediction_single"] = tidal_norm(
                first_order
            ) / tidal_norm(truth)
            rows.append(item)
        for item in misalignment(window_set, specs, truth):
            item["wavelength"] = wavelength
            turned.append(item)
        context["window_sets"][wavelength] = window_set
        context["estimated_k"][wavelength] = covector(n)
    result["windows"] = per_window
    result["sixth_degree"] = sixth
    result["principal_null_direction"] = pnd
    result["corrections"] = rows
    result["misaligned_side_windows"] = turned
    # Bounds.
    scan = plane_wave_scan(alias_n, unit)
    operator_norm = float(
        np.linalg.svd(TIDAL_ROOT @ alias42_n @ root_inverse, compute_uv=False)[0]
    )
    result["bounds"] = {
        "plane_wave_first_order_per_unit_f2_over_f": scan,
        "plane_wave_worst_at_wavelength_2": scan["worst"] * np.pi**2,
        "plane_wave_benchmark_at_wavelength_2": scan[
            "benchmark_direction_benchmark_polarization"
        ]
        * np.pi**2,
        "general_operator_norm": operator_norm,
        "general_bound_at_plane_wave_regularity_wavelength_2": 2
        * operator_norm
        * np.pi**2,
        "regularity_norm": (
            "sum over ordered (e, f) in the observer frame of the full-tidal norm "
            "squared of D_ef; a plane wave has |D| = 2 |f''/f| |C|"
        ),
    }
    result["displacements"] = displacement_study(alias42_n, root_inverse)
    # Worldline-only: the part of L that the time-time slice cannot fix.
    lhs = TIDAL_ROOT @ alias42_n @ root_inverse
    whitened = tt @ root_inverse
    projector = np.eye(42) - np.linalg.pinv(whitened) @ whitened
    result["bounds"]["worldline_undetermined_operator_norm"] = float(
        np.linalg.svd(lhs @ projector, compute_uv=False)[0]
    )
    return result, context


def mixed_geometry(own, other):
    """Centre window from one landmark layout, side windows from the other.

    Each window is fitted with its own landmarks, clocks and covariance and
    mapped to the generator frame with its own (oracle) frame; the correction
    uses the centre layout's alias map. Wavelength 2 only.
    """
    window_set = dict(other["window_sets"][2.0])
    window_set[0.0] = own["window_sets"][2.0][0.0]
    specs = correction_specs(
        sorted(window_set),
        own["alias_n"],
        np.array([1.0, 0, 0, 1]),
        own["estimated_k"][2.0],
        own["minimum_norm"],
    )
    return evaluate_specs(window_set, specs, own["truth"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--long-cases", type=int, nargs="*", default=[0])
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--selfcheck-only", action="store_true")
    args = parser.parse_args()
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_landmark_design import code_state

    amplitude = benchmark_amplitude()
    commit, clean = code_state()
    report = selfcheck(amplitude)
    print(json.dumps(report, indent=1), flush=True)
    if args.selfcheck_only:
        return
    args.cache.mkdir(parents=True, exist_ok=True)
    cache_checks = {}
    jobs = []
    for case in args.cases:
        path = clock_cache(args.cache, args.archive, case)
        cache_checks[case] = check_cache(path, args.archive, case)
        if not all(cache_checks[case].values()):
            raise ValueError(f"Clock cache of case {case} does not match the archive")
        jobs += jobs_for(
            case,
            args.archive,
            args.cache,
            amplitude,
            LONG_WAVELENGTHS if case in args.long_cases else (),
        )
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[name] = str(args.threads)
    begin = time.time()
    windows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for window in pool.map(window_job, jobs):
            windows.append(window)
            print(
                window["case"],
                window["kind"],
                window["wavelength"],
                window["shift"],
                *(f"o{o} {window[f'order_{o}']['seconds']}s" for o in window["orders"]),
                flush=True,
            )
    fit_seconds = round(time.time() - begin, 1)
    analysed = {
        case: analyse_case(
            case,
            [w for w in windows if w["case"] == case],
            args.cache,
            float(amplitude),
        )
        for case in args.cases
    }
    cases = []
    for case in args.cases:
        result, context = analysed[case]
        others = [c for c in args.cases if c != case]
        if others:
            result["mixed_geometry"] = {
                "side_windows_from_case": others[0],
                "rows": mixed_geometry(context, analysed[others[0]][1]),
            }
        cases.append(result)
    summary = {
        "diagnostic": (
            "multi-window bias correction and bias bounds; deterministic continuum; "
            "oracle frame alignment; no criterion"
        ),
        "code_version": commit,
        "worktree_clean": clean,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "benchmark_amplitude": str(amplitude),
        "reference_n": REFERENCE_N,
        "selfcheck": report,
        "cache_checks": {str(k): v for k, v in cache_checks.items()},
        "fit_seconds": fit_seconds,
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    raw = args.output.with_name(args.output.stem + "_windows.json")
    raw.write_text(json.dumps(windows) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
