"""Local vacuum geometry and Weyl fit using only interval cardinalities.

The two outer tips set the observer and the dimensionless length unit.
This development fitter has no access to generator coordinates or profiles.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from intrinsic_weyl import ETA, PAIRS, WEYL, gradient_basis
from scipy.optimize import least_squares
from scipy.stats import chi2

BIVECTOR = np.array(
    [[[c[a, b, d, e] for d, e in PAIRS] for a, b in PAIRS] for c in WEYL]
)
_TIDAL_BASIS = WEYL[:5, 0, 1:, 0, 1:].reshape(5, -1)
TIDAL_GRAM = np.kron(np.eye(2), _TIDAL_BASIS @ _TIDAL_BASIS.T)
TIDAL_ROOT = np.linalg.cholesky(TIDAL_GRAM).T


def axis_signs(signs):
    """Coefficient transform under diagonal O(3), including magnetic parity."""
    x, y, z = signs
    electric = np.array([1, 1, x * y, x * z, y * z])
    return np.r_[electric, x * y * z * electric]


def confidence_ball(fit, coverage=0.95, relative_bias=0.02):
    """Nominal Gaussian full-tidal ball modulo the remaining z reflection.

    The first two landmark vectors fix x and the xy half-plane. No third
    landmark is used to choose handedness. The bias allowance is an explicit
    model assumption, not a bound inferred from an unknown wave profile.
    """
    if not 0 < coverage < 1 or not 0 <= relative_bias < 1:
        raise ValueError("Invalid confidence or relative-bias level")
    if (
        not fit["optimizer_success"]
        or fit["boundary_excess"] > 0
        or fit["coordinate_rank"] != fit["coordinate_parameters"]
        or fit["projected_curvature_rank"] != 10
        or fit.get("nuisance_rank", fit["coordinate_rank"])
        != fit.get("nuisance_parameters", fit["coordinate_parameters"])
    ):
        raise ValueError("Incomplete or rank-deficient fit")
    points = np.asarray(fit["coordinates_inferred_from_counts"])
    if min(abs(points[1, 1]), abs(points[2, 2])) <= 1e-12:
        raise ValueError("Degenerate internal frame")
    transform = axis_signs([np.sign(points[1, 1]), np.sign(points[2, 2]), 1])
    coefficients = transform * np.asarray(fit["coefficients"])
    covariance = (
        transform[:, None] * np.asarray(fit["nominal_component_covariance"]) * transform
    )
    if not np.isfinite(covariance).all() or not np.isfinite(coefficients).all():
        raise ValueError("Nonfinite tidal estimate or covariance")
    eigenvalues = np.linalg.eigvalsh(TIDAL_ROOT @ covariance @ TIDAL_ROOT.T)
    if eigenvalues[0] <= 0:
        raise ValueError("Positive definite full-tidal covariance required")
    norm = float(np.linalg.norm(TIDAL_ROOT @ coefficients))
    statistical = float(np.sqrt(chi2.ppf(coverage, 10) * eigenvalues[-1]))
    return {
        "coefficients_xy_positive": coefficients.tolist(),
        "nominal_covariance_xy_positive": covariance.tolist(),
        "tidal_norm": norm,
        "largest_tidal_variance": float(eigenvalues[-1]),
        "nominal_coverage": coverage,
        "assumed_relative_bias_bound": relative_bias,
        "statistical_radius": statistical,
        "radius": (statistical + relative_bias * norm) / (1 - relative_bias),
        "scope": (
            "Nominal Gaussian ball with assumed bias bound; z reflection quotiented"
        ),
    }


def tidal_distance(left, right):
    """Full-tidal distance modulo reflection, with x and y signs fixed."""
    left, right = np.asarray(left), np.asarray(right)
    reflection = axis_signs([1, 1, -1])
    return min(
        float(np.linalg.norm(TIDAL_ROOT @ (left - candidate)))
        for candidate in (right, reflection * right)
    )


def select_landmarks(selection_counts, maximum=64, minimum_count=2, design="two"):
    """Peel poorly connected candidates, using an independent count sample.

    The first and last events are the outer diamond tips. Candidate IDs carry
    no geometric information. Exact score ties retain their input ordering.
    """
    counts = np.asarray(selection_counts, dtype=float)
    if counts.ndim != 2 or counts.shape[0] != counts.shape[1]:
        raise ValueError("Expected a square matrix of interval counts")
    if not np.isfinite(counts).all() or np.any(counts < 0):
        raise ValueError("Counts must be finite and nonnegative")
    if not np.allclose(counts, counts.T) or np.any(np.diag(counts)):
        raise ValueError("Counts must be symmetric with zero diagonal")
    total = counts[0, -1]
    if total <= 0:
        raise ValueError("The outer interval must contain events")
    if design == "nested":
        base = select_landmarks(counts, maximum, minimum_count, "connected")
        expanded = select_landmarks(counts, 2 * maximum, minimum_count, "connected")
        extra = expanded[~np.isin(expanded, base)]
        # Preserve every old clock and the two vectors defining the frame.
        # The union has at most 3*maximum interior landmarks, not exactly 2*maximum.
        return np.r_[base[:-1], extra, base[-1]]
    tau_p, tau_q = np.sqrt(counts[0] / total), np.sqrt(counts[-1] / total)
    times = (tau_p - tau_q) / 2
    radius2 = 0.25 + times**2 - (tau_p + tau_q) / 2
    if design == "wide":
        # Keep the original outer landmarks, including the two frame vectors.
        # Only central spatial spread changes in the prospective paired pilot.
        base = select_landmarks(counts, maximum, minimum_count, "connected")[1:-1]
        past, future = base[times[base] <= -0.24], base[times[base] >= 0.24]
        middle = np.flatnonzero(
            (abs(times) <= 0.075)
            & (radius2 >= 0.075**2)
            & (radius2 <= 0.3**2)
            & (counts[0] >= minimum_count)
            & (counts[-1] >= minimum_count)
        )
        degrees = np.column_stack(
            [
                (counts[middle[:, None], end] >= minimum_count).sum(axis=1)
                for end in (past, future)
            ]
        )
        keep = np.all(degrees >= 4, axis=1)
        middle, degrees = middle[keep], degrees[keep]
        score = degrees.min(axis=1) * radius2[middle]
        ranked = np.argsort(-score, kind="stable")
        middle = middle[ranked[: maximum - len(past) - len(future)]]
        ids = np.r_[past, middle, future]
        if len(middle) < 8 or np.any(
            (counts[ids[:, None], ids] >= minimum_count).sum(axis=1) < 8
        ):
            raise ValueError("Wider centers require eight centers and degree eight")
        return np.r_[0, ids, len(counts) - 1]
    if design not in ("two", "three", "connected"):
        raise ValueError("Unknown prospective selection design")
    outer_floor = 0.24 if design == "connected" else 0.15
    outer_ceiling = 0.36 if design == "connected" else 0.32
    outer = (abs(times) >= outer_floor) & (abs(times) <= outer_ceiling)
    eligible_time = outer if design == "two" else outer | (abs(times) <= 0.075)
    eligible_radius = (
        (radius2 >= 0.075**2) & (radius2 <= 0.2**2)
        if design == "connected"
        else radius2 >= 0.01
    )
    ids = np.flatnonzero(
        eligible_time
        & eligible_radius
        & (counts[0] >= minimum_count)
        & (counts[-1] >= minimum_count)
    )
    related = counts[ids[:, None], ids] >= minimum_count
    score = related.sum(axis=1) * np.minimum(radius2[ids], 0.04)
    if design == "connected":
        third = maximum // 3
        ends = []
        for condition in (times[ids] <= -outer_floor, times[ids] >= outer_floor):
            choices = np.flatnonzero(condition)
            ranked = choices[np.argsort(-score[choices], kind="stable")]
            ends.append(ids[ranked[:third]])
        middle = ids[abs(times[ids]) <= 0.075]
        degrees = np.column_stack(
            [
                (counts[middle[:, None], end] >= minimum_count).sum(axis=1)
                for end in ends
            ]
        )
        connected = np.all(degrees >= 4, axis=1)
        middle, degrees = middle[connected], degrees[connected]
        score = degrees.min(axis=1) * np.minimum(radius2[middle], 0.04)
        ranked = np.argsort(-score, kind="stable")
        ids = np.r_[ends[0], middle[ranked[: maximum - 2 * third]], ends[1]]
        strata = []
    elif design == "two":
        strata = [(times[ids] < 0, maximum // 2), (times[ids] >= 0, maximum // 2)]
    else:
        third = maximum // 3
        strata = [
            (times[ids] <= -0.15, third),
            (abs(times[ids]) <= 0.075, maximum - 2 * third),
            (times[ids] >= 0.15, third),
        ]
    selected = []
    for condition, limit in strata:
        choices = np.flatnonzero(condition)
        ranked = choices[np.argsort(-score[choices], kind="stable")]
        selected.append(ids[ranked[:limit]])
    if selected:
        ids = np.concatenate(selected)
    while len(ids):
        graph = counts[ids[:, None], ids] >= minimum_count
        degrees = graph.sum(axis=1)
        if degrees.min() >= 8 and len(ids) <= maximum:
            break
        score = degrees + np.log1p(counts[ids[:, None], ids]).sum(axis=1) * 1e-8
        ids = np.delete(ids, int(np.argmin(score)))
    if len(ids) < 16:
        raise ValueError("Insufficient mutually informative landmarks")
    if design == "connected":
        past, middle, future = (
            ids[times[ids] <= -outer_floor],
            ids[abs(times[ids]) <= 0.075],
            ids[times[ids] >= outer_floor],
        )
        if len(past) < 16 or len(middle) < 8 or len(future) < 16:
            raise ValueError(
                "Connected design requires at least 16/8/16 retained landmarks"
            )
        if any(
            np.any((counts[middle[:, None], end] >= minimum_count).sum(axis=1) < 4)
            for end in (past, future)
        ):
            raise ValueError(
                "Every central landmark needs four retained clocks to each side"
            )
    return np.r_[0, ids, len(counts) - 1]


def gauge(points):
    """Rotate an inferred spatial embedding; no external frame is supplied."""
    first = points[1, 1:]
    first = first / np.linalg.norm(first)
    second = points[2, 1:] - np.dot(points[2, 1:], first) * first
    if np.linalg.norm(second) < 1e-10:
        raise ValueError("Degenerate initial spatial frame")
    second /= np.linalg.norm(second)
    frame = np.column_stack((first, second, np.cross(first, second)))
    result = points.copy()
    result[:, 1:] = result[:, 1:] @ frame
    result[1, 2:] = 0
    result[2, 3] = 0
    return result


def model(points, coefficients, pairs, *, square_correction=False):
    """Quartic or fifth-degree world function with analytic derivatives."""
    if square_correction:
        from intrinsic_weyl_square import correction

        return tuple(
            a + b
            for a, b in zip(
                model(points, coefficients, pairs),
                correction(points, coefficients, pairs),
                strict=True,
            )
        )
    left, right = points[pairs[:, 0]], points[pairs[:, 1]]
    delta = left - right
    flat = np.einsum("pa,ab,pb->p", delta, ETA, delta)
    wedge = np.column_stack(
        [left[:, a] * right[:, b] - left[:, b] * right[:, a] for a, b in PAIRS]
    )
    basis = np.einsum("pi,kij,pj->pk", wedge, BIVECTOR, wedge, optimize=True)
    gradient = None
    if len(coefficients) == 34:
        gradient = (gradient_basis() @ coefficients[10:]).reshape(4, 10)
        effective = coefficients[:10] + (left + right) @ gradient / 4
        curvature = np.einsum("pk,kij->pij", effective, BIVECTOR)
        cw = np.einsum("pi,pij->pj", wedge, curvature)
    elif len(coefficients) == 10:
        curvature = np.einsum("k,kij->ij", coefficients, BIVECTOR)
        cw = wedge @ curvature
    else:
        raise ValueError("Ten Weyl or ten Weyl plus 24 gradient coefficients required")
    dl = 2 * delta @ ETA
    dr = -dl.copy()
    for k, (a, b) in enumerate(PAIRS):
        dl[:, a] -= 2 * cw[:, k] * right[:, b] / 3
        dl[:, b] += 2 * cw[:, k] * right[:, a] / 3
        dr[:, a] += 2 * cw[:, k] * left[:, b] / 3
        dr[:, b] -= 2 * cw[:, k] * left[:, a] / 3
    dc = -basis / 3
    if gradient is None:
        return flat - basis @ coefficients / 3, dl, dr, dc
    slope = -basis @ gradient.T / 12
    dl += slope
    dr += slope
    dg = (
        np.einsum("pe,pk->pek", (left + right) / 4, dc).reshape(len(pairs), 40)
        @ gradient_basis()
    )
    value = flat - np.einsum("pk,pk->p", basis, effective) / 3
    return value, dl, dr, np.column_stack((dc, dg))


def fit_joint(
    counts,
    minimum_count=32,
    starts=4,
    clock_mask=None,
    clock_covariance=None,
    model_order=4,
    square_correction=False,
):
    """Jointly infer positions and all ten Weyl components from counts.

    Input landmarks are selected separately. A returned point estimate is not
    a confidence statement. The local truncation assumes a small vacuum
    diamond and the outer-tip observer; nonlocal and matter cases are excluded.
    """
    if model_order not in (4, 5):
        raise ValueError("World-function order must be four or five")
    if square_correction and model_order != 5:
        raise ValueError("The curvature-square variant requires the fifth-degree fit")
    curvature_parameters = 10 if model_order == 4 else 34
    counts = np.asarray(counts, dtype=float)
    if (
        counts.ndim != 2
        or counts.shape[0] != counts.shape[1]
        or len(counts) < 18
        or not np.isfinite(counts).all()
        or np.any(counts < 0)
        or not np.allclose(counts, counts.T)
        or np.any(np.diag(counts))
    ):
        raise ValueError("Finite symmetric nonnegative interval matrix required")
    n = len(counts)
    if (
        counts[0, -1] <= 0
        or np.any(counts[0, 1:-1] <= 0)
        or np.any(counts[-1, 1:-1] <= 0)
    ):
        raise ValueError("Every landmark requires clocks to both outer tips")
    available = np.triu(
        counts >= minimum_count if clock_mask is None else clock_mask, 1
    )
    available[0, -1] = False
    pairs = np.column_stack(np.nonzero(available))
    if np.any(counts[pairs[:, 0], pairs[:, 1]] <= 0):
        raise ValueError("A preselected clock has zero measurement events")
    target = -np.sqrt(counts[pairs[:, 0], pairs[:, 1]] / counts[0, -1])
    whiten = None
    if clock_covariance is not None:
        covariance = np.asarray(clock_covariance)
        if (
            covariance.shape != (len(pairs), len(pairs))
            or not np.isfinite(covariance).all()
        ):
            raise ValueError("Clock covariance has the wrong shape or nonfinite values")
        scales = np.sqrt(np.diag(covariance))
        if np.any(scales <= 0):
            raise ValueError("Clock variances must be positive")
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
    if len(pairs) <= dimension + curvature_parameters:
        raise ValueError("Not enough measured intervals for a curvature fit")

    def unpack(parameters):
        points = np.zeros((n, 4))
        points[0, 0], points[-1, 0] = -0.5, 0.5
        points[mask] = parameters[:dimension]
        c = parameters[dimension:] if len(parameters) > dimension else np.zeros(10)
        return points, c

    def residual(parameters):
        value = (
            model(*unpack(parameters), pairs, square_correction=square_correction)[0]
            - target
        )
        return value if whiten is None else whiten @ value

    def jacobian(parameters):
        _, dl, dr, dc = model(
            *unpack(parameters), pairs, square_correction=square_correction
        )
        result = np.zeros((len(pairs), len(parameters)))
        for side, derivative in enumerate((dl, dr)):
            for component in range(4):
                index = parameter_index[pairs[:, side], component]
                valid = index >= 0
                result[np.flatnonzero(valid), index[valid]] = derivative[
                    valid, component
                ]
        if len(parameters) > dimension:
            result[:, dimension:] = dc
        return result if whiten is None else whiten @ result

    tau_p = np.sqrt(counts[0, 1:-1] / counts[0, -1])
    tau_q = np.sqrt(counts[-1, 1:-1] / counts[0, -1])
    times = (tau_p - tau_q) / 2
    radii = np.sqrt(np.maximum(0.25 + times**2 - (tau_p + tau_q) / 2, 1e-8))
    candidates = []
    # ponytail: four deterministic starts expose local minima; no optimizer framework.
    for seed in range(starts):
        rng = np.random.default_rng(seed)
        directions = rng.normal(size=(n - 2, 3))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        points = np.zeros((n, 4))
        points[0, 0], points[-1, 0] = -0.5, 0.5
        points[1:-1, 0] = times
        points[1:-1, 1:] = radii[:, None] * directions
        initial = gauge(points)[mask]
        fit = least_squares(
            residual,
            initial,
            jac=jacobian,
            method="lm",
            max_nfev=600,
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-10,
        )
        candidates.append(fit)
    initial = min(candidates, key=lambda fit: fit.cost)
    curved = least_squares(
        residual,
        np.r_[initial.x, np.zeros(curvature_parameters)],
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
    nuisance = (
        np.column_stack((geometry, jac[:, dimension + 10 :]))
        if model_order == 5
        else geometry
    )
    projected = (
        curvature - nuisance @ np.linalg.lstsq(nuisance, curvature, rcond=None)[0]
    )
    singular = np.linalg.svd(projected, compute_uv=False)
    tensor = np.einsum("k,kabcd->abcd", c[:10], WEYL)
    result = {
        "coefficients": c[:10].tolist(),
        "electric_eigenvalues": np.linalg.eigvalsh(tensor[0, 1:, 0, 1:]).tolist(),
        "electric_norm_squared": float(np.sum(tensor[0, 1:, 0, 1:] ** 2)),
        "coordinates_inferred_from_counts": points.tolist(),
        "landmarks": n - 2,
        "measured_intervals": len(pairs),
        "coordinate_rank": int(np.linalg.matrix_rank(geometry)),
        "coordinate_parameters": dimension,
        "projected_curvature_rank": int(np.linalg.matrix_rank(projected)),
        "projected_curvature_singular_values": singular.tolist(),
        "flat_costs": [float(fit.cost) for fit in candidates],
        "curved_cost": float(curved.cost),
        "optimizer_success": bool(curved.success),
        "optimizer_message": curved.message,
        "function_evaluations": int(curved.nfev),
        "boundary_excess": float(max(0, np.max(abs(points[:, 0])) - 0.5)),
    }
    if model_order == 5:
        result.update(
            {
                "model_order": 5,
                "gradient_coefficients": c[10:].tolist(),
                "covariant_gradient_components": (gradient_basis() @ c[10:])
                .reshape(4, 10)
                .tolist(),
                "nuisance_rank": int(np.linalg.matrix_rank(nuisance)),
                "nuisance_parameters": nuisance.shape[1],
            }
        )
    if whiten is not None:
        inverse = np.linalg.pinv(jac)
        influence = inverse[dimension : dimension + 10]
        component_covariance = influence @ influence.T
        result["nominal_component_covariance"] = component_covariance.tolist()
        result["uncertainty_scope"] = (
            "Local linearized covariance under the quartic model; "
            if model_order == 4
            else "Local linearized covariance under the fifth-degree model; "
        ) + "coverage and truncation bias not certified"
        result["retained_clock_covariance_rank"] = len(whiten)
    if square_correction:
        result["curvature_square_correction"] = True
        result["model_scope"] = (
            "Fifth-degree jet plus curvature-squared world-function and volume-clock "
            "terms; second covariant curvature derivatives still omitted"
        )
    return result


def selection(order, design="two"):
    """Order-only choice of landmarks and of informative pair clocks."""
    order = np.asarray(order, dtype=bool)
    if order.ndim != 2 or order.shape[0] != order.shape[1] or np.any(np.diag(order)):
        raise ValueError("Strict selection order required")
    if not np.all(order[0, 1:]) or not np.all(order[:-1, -1]):
        raise ValueError("Selection events must be in the supplied outer interval")
    counts = order.astype(float) @ order.astype(float)
    counts += counts.T.copy()
    ids = select_landmarks(counts, design=design)
    restricted = counts[ids[:, None], ids]
    clock_mask = restricted >= 2
    clock_mask[0, -1] = clock_mask[-1, 0] = False
    return ids, clock_mask


def _clock_second_moment(after, before, a, b):
    """Sum repeated causal-membership patterns with exact integer weights."""
    landmarks, bytes_per_landmark = after.shape
    n = 8 * bytes_per_landmark
    # ponytail: keep signatures in memory, like the packed input; external
    # grouping is needed only if this linear memory cost exceeds the run budget.
    signatures = np.empty((n, (2 * landmarks + 7) // 8), dtype=np.uint8)
    for start in range(0, bytes_per_landmark, 4096):
        end = min(start + 4096, bytes_per_landmark)
        bits = np.unpackbits(
            np.vstack((after[:, start:end], before[:, start:end])), axis=1
        )
        signatures[8 * start : 8 * end] = np.packbits(bits, axis=0).T
    width = signatures.shape[1]
    unique, weights = np.unique(
        signatures.view(np.dtype(f"V{width}")).ravel(), return_counts=True
    )
    unique = unique.view(np.uint8).reshape(-1, width)
    second = np.zeros((len(a), len(a)))
    for start in range(0, len(unique), 4096):
        bits = np.unpackbits(unique[start : start + 4096], axis=1)
        membership = (bits[:, a] & bits[:, landmarks + b]).T.astype(float)
        # No sqrt(weights): every product and partial sum stays an exactly
        # representable integer when the event count is below 2**53.
        second += (membership * weights[start : start + 4096]) @ membership.T
    return second


def read_order_clocks(path):
    """Derive clocks and their covariance from the saved order bits."""
    with np.load(path, allow_pickle=False) as source:
        allowed = {
            "selection_order",
            "selected_ids",
            "after",
            "before",
            "measurement_n",
        }
        if set(source.files) not in (allowed, allowed | {"selection_design"}):
            raise ValueError("Archive contains unexpected fields")
        order, ids = source["selection_order"], source["selected_ids"]
        after, before = source["after"], source["before"]
        n = int(source["measurement_n"])
        design = (
            str(source["selection_design"].item())
            if "selection_design" in source
            else "two"
        )
    expected, clock_mask = selection(order, design=design)
    if not np.array_equal(ids, expected):
        raise ValueError("Saved landmarks differ from the order-only selection")
    if after.shape != (len(ids), n // 8) or before.shape != after.shape or n % 8:
        raise ValueError("Packed order dimensions are inconsistent")
    if not (
        np.all(after[0] == 255)
        and np.all(before[-1] == 255)
        and np.all(before[0] == 0)
        and np.all(after[-1] == 0)
    ):
        raise ValueError("Measurement order is not contained in the outer interval")
    counts = np.zeros((len(ids), len(ids)))
    restricted_order = order[ids[:, None], ids]
    for i, j in itertools.combinations(range(len(ids)), 2):
        if restricted_order[i, j]:
            a, b = i, j
        elif restricted_order[j, i]:
            a, b = j, i
        else:
            continue
        counts[i, j] = counts[j, i] = np.bitwise_count(after[a] & before[b]).sum()
    pairs = np.column_stack(np.nonzero(np.triu(clock_mask, 1)))
    forward = restricted_order[pairs[:, 0], pairs[:, 1]]
    a = np.where(forward, pairs[:, 0], pairs[:, 1])
    b = np.where(forward, pairs[:, 1], pairs[:, 0])
    p = counts[pairs[:, 0], pairs[:, 1]] / n
    if np.any(p <= 0):
        raise ValueError("A selected clock has zero measurement count")
    second = _clock_second_moment(after, before, a, b)
    indicator_covariance = (second - n * np.outer(p, p)) / (n - 1)
    derivative = -0.5 / np.sqrt(p)
    clock_covariance = indicator_covariance * np.outer(derivative, derivative) / n
    return counts, clock_mask, clock_covariance, ids, n


def fit_order_archive(
    path, *, with_clocks=False, model_order=4, square_correction=False
):
    """Blind boundary: accept only saved order bits; infer all other inputs."""
    clocks = read_order_clocks(path)
    counts, clock_mask, clock_covariance, ids, n = clocks
    result = fit_joint(
        counts,
        clock_mask=clock_mask,
        clock_covariance=clock_covariance,
        model_order=model_order,
        square_correction=square_correction,
    )
    result["selected_ids"] = ids.tolist()
    result["measurement_n"] = n
    result["smallest_measurement_clock"] = int(counts[clock_mask].min())
    result["order_archive_sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return (result, clocks) if with_clocks else result


def selfcheck_jet():
    """One deterministic check of the full fifth-degree inverse path."""
    rng = np.random.default_rng(719)
    for landmarks, n, distinct in ((13, 8192, 7), (9, 512, 512), (11, 256, 1)):
        patterns = rng.integers(0, 2, (2 * landmarks, distinct), dtype=np.uint8)
        bits = patterns[:, np.arange(n) % distinct]
        a, b = np.arange(landmarks - 1), np.arange(1, landmarks)
        membership = (bits[a] & bits[landmarks + b]).astype(float)
        grouped = _clock_second_moment(
            np.packbits(bits[:landmarks], axis=1),
            np.packbits(bits[landmarks:], axis=1),
            a,
            b,
        )
        assert np.array_equal(grouped, membership @ membership.T)
    points = 0.1 * np.sin(np.arange(32).reshape(8, 4))
    coefficients = np.linspace(-0.2, 0.3, 34)
    pairs = np.array(list(itertools.combinations(range(8), 2)))
    value, dl, dr, dc = model(points, coefficients, pairs)
    maximum = 0.0
    for i in range(34):
        step = np.zeros(34)
        step[i] = 1e-6
        numeric = (
            model(points, coefficients + step, pairs)[0]
            - model(points, coefficients - step, pairs)[0]
        ) / 2e-6
        maximum = max(maximum, float(np.max(abs(numeric - dc[:, i]))))
    for i, a in ((1, 0), (1, 2), (4, 1), (6, 3)):
        step = np.zeros_like(points)
        step[i, a] = 1e-6
        numeric = (
            model(points + step, coefficients, pairs)[0]
            - model(points - step, coefficients, pairs)[0]
        ) / 2e-6
        expected = (pairs[:, 0] == i) * dl[:, a] + (pairs[:, 1] == i) * dr[:, a]
        maximum = max(maximum, float(np.max(abs(numeric - expected))))
    assert maximum < 1e-9 and np.isfinite(value).all()
    reduced = model(points, np.r_[coefficients[:10], np.zeros(24)], pairs)
    original = model(points, coefficients[:10], pairs)
    assert all(
        np.max(abs(a - b)) < 1e-13
        for a, b in zip(reduced[:3], original[:3], strict=True)
    )
    k = np.arange(48)
    z = 1 - 2 * (k + 0.5) / len(k)
    phi = k * np.pi * (3 - np.sqrt(5))
    directions = np.column_stack(
        (np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z)
    )
    points = np.zeros((50, 4))
    points[0, 0], points[-1, 0] = -0.5, 0.5
    points[1:-1, 0] = 0.27 * ((k % 3) - 1) + 0.015 * np.sin(k)
    points[1:-1, 1:] = (0.12 + 0.08 * np.sin(k * 2.1) ** 2)[:, None] * directions
    delta = points[:, None] - points[None]
    tau2 = delta[..., 0] ** 2 - np.sum(delta[..., 1:] ** 2, axis=-1)
    fitted = fit_joint(1e8 * np.maximum(tau2, 0) ** 2, starts=2, model_order=5)
    assert fitted["projected_curvature_rank"] == 10
    assert fitted["nuisance_rank"] == fitted["nuisance_parameters"]
    assert np.linalg.norm(fitted["coefficients"]) < 1e-7
    assert np.linalg.norm(fitted["gradient_coefficients"]) < 1e-6
    return {
        "shared_interval_moments_exact": True,
        "maximum_jacobian_error": maximum,
        "exact_flat_weyl_rank": 10,
        "exact_flat_full_nuisance_rank": fitted["nuisance_rank"],
        "exact_flat_coefficient_norm": float(np.linalg.norm(fitted["coefficients"])),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    print(json.dumps(fit_order_archive(args.archive), indent=2))
