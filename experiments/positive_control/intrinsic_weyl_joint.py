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
from intrinsic_weyl import ETA, PAIRS, WEYL
from scipy.optimize import least_squares

BIVECTOR = np.array(
    [[[c[a, b, d, e] for d, e in PAIRS] for a, b in PAIRS] for c in WEYL]
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
    tau_p, tau_q = np.sqrt(counts[0] / total), np.sqrt(counts[-1] / total)
    times = (tau_p - tau_q) / 2
    radius2 = 0.25 + times**2 - (tau_p + tau_q) / 2
    if design not in ("two", "three"):
        raise ValueError("Unknown prospective selection design")
    outer = (abs(times) >= 0.15) & (abs(times) <= 0.32)
    eligible_time = outer if design == "two" else outer | (abs(times) <= 0.075)
    ids = np.flatnonzero(
        eligible_time
        & (radius2 >= 0.01)
        & (counts[0] >= minimum_count)
        & (counts[-1] >= minimum_count)
    )
    related = counts[ids[:, None], ids] >= minimum_count
    score = related.sum(axis=1) * np.minimum(radius2[ids], 0.04)
    if design == "two":
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


def model(points, coefficients, pairs):
    """Quartic world function and analytic position/curvature derivatives."""
    left, right = points[pairs[:, 0]], points[pairs[:, 1]]
    delta = left - right
    flat = np.einsum("pa,ab,pb->p", delta, ETA, delta)
    wedge = np.column_stack(
        [left[:, a] * right[:, b] - left[:, b] * right[:, a] for a, b in PAIRS]
    )
    basis = np.einsum("pi,kij,pj->pk", wedge, BIVECTOR, wedge, optimize=True)
    curvature = np.einsum("k,kij->ij", coefficients, BIVECTOR)
    cw = wedge @ curvature
    dl = 2 * delta @ ETA
    dr = -dl.copy()
    for k, (a, b) in enumerate(PAIRS):
        dl[:, a] -= 2 * cw[:, k] * right[:, b] / 3
        dl[:, b] += 2 * cw[:, k] * right[:, a] / 3
        dr[:, a] += 2 * cw[:, k] * left[:, b] / 3
        dr[:, b] -= 2 * cw[:, k] * left[:, a] / 3
    return flat - basis @ coefficients / 3, dl, dr, -basis / 3


def fit_joint(
    counts, minimum_count=32, starts=4, clock_mask=None, clock_covariance=None
):
    """Jointly infer positions and all ten Weyl components from counts.

    Input landmarks are selected separately. A returned point estimate is not
    a confidence statement. The quartic truncation assumes a small vacuum
    diamond and the outer-tip observer; nonlocal and matter cases are excluded.
    """
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
    if len(pairs) <= dimension + 10:
        raise ValueError("Not enough measured intervals for a curvature fit")

    def unpack(parameters):
        points = np.zeros((n, 4))
        points[0, 0], points[-1, 0] = -0.5, 0.5
        points[mask] = parameters[:dimension]
        c = parameters[dimension:] if len(parameters) > dimension else np.zeros(10)
        return points, c

    def residual(parameters):
        value = model(*unpack(parameters), pairs)[0] - target
        return value if whiten is None else whiten @ value

    def jacobian(parameters):
        _, dl, dr, dc = model(*unpack(parameters), pairs)
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
        np.r_[initial.x, np.zeros(10)],
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
    curvature = jac[:, dimension:]
    projected = (
        curvature - geometry @ np.linalg.lstsq(geometry, curvature, rcond=None)[0]
    )
    singular = np.linalg.svd(projected, compute_uv=False)
    tensor = np.einsum("k,kabcd->abcd", c, WEYL)
    result = {
        "coefficients": c.tolist(),
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
    if whiten is not None:
        inverse = np.linalg.pinv(jac)
        component_covariance = inverse[-10:] @ inverse[-10:].T
        result["nominal_component_covariance"] = component_covariance.tolist()
        result["uncertainty_scope"] = (
            "Local linearized covariance under the quartic model; "
            "coverage and truncation bias not certified"
        )
        result["retained_clock_covariance_rank"] = len(whiten)
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
    second = np.zeros((len(pairs), len(pairs)))
    for start in range(0, n // 8, 4096 // 8):
        membership = np.unpackbits(
            after[a, start : start + 4096 // 8] & before[b, start : start + 4096 // 8],
            axis=1,
        ).astype(float)
        second += membership @ membership.T
    indicator_covariance = (second - n * np.outer(p, p)) / (n - 1)
    derivative = -0.5 / np.sqrt(p)
    clock_covariance = indicator_covariance * np.outer(derivative, derivative) / n
    return counts, clock_mask, clock_covariance, ids, n


def fit_order_archive(path):
    """Blind boundary: accept only saved order bits; infer all other inputs."""
    counts, clock_mask, clock_covariance, ids, n = read_order_clocks(path)
    result = fit_joint(counts, clock_mask=clock_mask, clock_covariance=clock_covariance)
    result["selected_ids"] = ids.tolist()
    result["measurement_n"] = n
    result["smallest_measurement_clock"] = int(counts[clock_mask].min())
    result["order_archive_sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    print(json.dumps(fit_order_archive(args.archive), indent=2))
