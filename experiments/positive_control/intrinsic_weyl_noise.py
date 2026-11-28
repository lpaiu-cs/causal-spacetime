"""Prospective flat-limit cost diagnostic; NOT a curved recovery campaign.

Protocol: docs/paper/paper_a/intrinsic_weyl/README.md section 5.
The flat generator is a design tool, not an input to a curved-data estimator.
"""

import hashlib
import itertools
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path
from statistics import NormalDist

import numpy as np
from intrinsic_weyl import ETA, WEYL, fit_from_interval_counts, residual_design
from intrinsic_weyl_development import TANGENTS, physical_profile, target_tensor

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from causal_spacetime_lab.sprinkling import (  # noqa: E402
    sprinkle_minkowski_causal_diamond,
)

OUTER_T = 11.0
PER_SEED = 262144
CHUNK = 16384
SEEDS = range(2026092200, 2026092204)


def linearized_fit():
    """Analytic Jacobian with respect to interval occupancy probabilities."""
    delta = TANGENTS[:, None] - TANGENTS[None, :]
    tau2 = -np.einsum("ija,ab,ijb->ij", delta, ETA, delta)
    reference_tau2 = tau2[0, 4]
    points = TANGENTS / math.sqrt(reference_tau2)
    references = points[1:5]
    coefficients = points[5:] @ np.linalg.inv(references)
    design, residual_pairs = residual_design(references, coefficients)
    keep = [i == j or tau2[i + 5, j + 5] > 0 for i, j in residual_pairs]
    design = design[keep]
    residual_pairs = [p for p, k in zip(residual_pairs, keep, strict=True) if k]
    clocks = [
        (i, j)
        for i, j in itertools.combinations(range(len(points)), 2)
        if tau2[i, j] > 0
    ]
    probabilities = np.array([(tau2[i, j] / OUTER_T**2) ** 2 for i, j in clocks])
    # At the flat limit D=0, so derivatives of the fitted design do not enter.
    projection = np.column_stack((-coefficients, np.eye(len(points) - 5)))
    sums = projection.sum(axis=1)
    derivatives = []
    for (a, b), probability in zip(clocks, probabilities, strict=True):
        wb = projection[:, b - 1]
        if a == 0:
            response = (np.outer(wb, sums) + np.outer(sums, wb)) / 2
        else:
            wa = projection[:, a - 1]
            response = -(np.outer(wa, wb) + np.outer(wb, wa)) / 2
        reference_probability = (reference_tau2 / OUTER_T**2) ** 2
        response *= -1 / (2 * math.sqrt(probability * reference_probability))
        derivatives.append([response[i, j] for i, j in residual_pairs])
    response = np.array(derivatives).T
    jacobian = np.linalg.pinv(design) @ response
    # The omitted derivative of the normalization multiplies the zero flat D.
    assert np.linalg.norm(jacobian @ probabilities) < 1e-8
    base = np.zeros((len(points), len(points)))
    direction = np.zeros_like(base)
    relative_direction = np.cos(np.arange(len(clocks)))
    for k, (a, b) in enumerate(clocks):
        base[a, b] = base[b, a] = probabilities[k]
        direction[a, b] = direction[b, a] = probabilities[k] * relative_direction[k]
    step = 1e-6
    plus = np.array(fit_from_interval_counts(base + step * direction)["coefficients"])
    minus = np.array(fit_from_interval_counts(base - step * direction)["coefficients"])
    numeric = (plus - minus) / (2 * step)
    analytic = jacobian @ (probabilities * relative_direction)
    discrepancy = np.linalg.norm(numeric - analytic) / np.linalg.norm(analytic)
    assert discrepancy < 2e-6, discrepancy
    return (
        clocks,
        probabilities,
        jacobian,
        reference_tau2,
        discrepancy,
        design,
        response,
    )


def memberships(rng, clocks):
    """Shared event membership, retaining correlations between interval counts."""
    a, b = np.array(clocks).T
    events = sprinkle_minkowski_causal_diamond(CHUNK, 4, OUTER_T, rng)
    events[:, 0] += OUTER_T / 2
    dt = events[:, None, 0] - TANGENTS[None, :, 0]
    distance = np.linalg.norm(events[:, None, 1:] - TANGENTS[None, :, 1:], axis=2)
    return (dt > distance)[:, a] & (-dt > distance)[:, b]


def price(covariance, reference_tau2):
    tensor_metric = np.einsum("kabcd,labcd->kl", WEYL, WEYL)
    normal = NormalDist()
    factor = ((normal.inv_cdf(0.975) + normal.inv_cdf(0.95)) / 0.25) ** 2
    rows = []
    for scale in (0.04, 0.02):
        truth = target_tensor(physical_profile("constant")) * reference_tau2 * scale**2
        norm2 = float(np.sum(truth**2))
        weights = np.einsum("kabcd,abcd->k", WEYL, truth) / norm2
        projection_variance = float(weights @ covariance @ weights)
        rows.append(
            {
                "scale": scale,
                "predicted_relative_tensor_rms_at_1e6": math.sqrt(
                    np.sum(tensor_metric * covariance) / (1e6 * norm2)
                ),
                "projection_relative_sd_at_1e6": math.sqrt(projection_variance / 1e6),
                "optimistic_single_projection_n_95ci_90power_25pct": math.ceil(
                    factor * projection_variance
                ),
            }
        )
    return rows


def run():
    # Freeze precedes random data; preserve a resolvable commit and a clean tree.
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True
    )
    if dirty.strip():
        raise RuntimeError("Commit the prospective protocol and implementation first")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    clocks, probabilities, jacobian, reference_tau2, discrepancy, _, _ = (
        linearized_fit()
    )
    a, b = np.array(clocks).T
    rows = []
    for seed in SEEDS:
        start = time.perf_counter()
        rng = np.random.default_rng(seed)
        count = np.zeros(len(clocks), dtype=np.int64)
        influence_sum = np.zeros(10)
        second = np.zeros((10, 10))
        for _ in range(PER_SEED // CHUNK):
            membership = memberships(rng, clocks)
            count += membership.sum(axis=0)
            influence = (membership - probabilities) @ jacobian.T
            influence_sum += influence.sum(axis=0)
            second += influence.T @ influence
        covariance = (second - np.outer(influence_sum, influence_sum) / PER_SEED) / (
            PER_SEED - 1
        )
        observed = np.zeros((len(TANGENTS), len(TANGENTS)))
        observed[a, b] = observed[b, a] = count
        try:
            fitted = fit_from_interval_counts(observed)
        except ValueError as error:
            fitted = {"failure": str(error)}
        row = {
            "seed": seed,
            "n": PER_SEED,
            "seconds": time.perf_counter() - start,
            "minimum_interval_count": int(count.min()),
            "influence_covariance": covariance.tolist(),
            "price": price(covariance, reference_tau2),
            "flat_fit_from_actual_counts": fitted,
        }
        rows.append(row)
        print(json.dumps({k: row[k] for k in ("seed", "seconds", "price")}), flush=True)
    sources = [
        Path(__file__),
        Path(__file__).with_name("intrinsic_weyl.py"),
        Path(__file__).with_name("intrinsic_weyl_development.py"),
        ROOT / "src/causal_spacetime_lab/sprinkling.py",
        ROOT / "docs/paper/paper_a/intrinsic_weyl/README.md",
    ]
    return {
        "kind": (
            "Flat-limit development cost projection; not curved recovery or coverage"
        ),
        "code_version": stamp,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "jacobian_directional_check_relative_error": discrepancy,
        "clock_count": len(clocks),
        "replicates": rows,
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }


if __name__ == "__main__":
    result = run()
    destination = ROOT / "docs/paper/paper_a/intrinsic_weyl/noise_cost_diagnostic.json"
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Written {destination}")
