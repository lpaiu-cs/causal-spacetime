"""Deterministic checks for plane_wave_chain_literature_and_limits.md.

Reuses the earlier Jacobi/Volterra implementation without changing its outputs.
No sampled posets, new dependency, or statistical accuracy claim.
"""

import hashlib
import itertools
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from check_plane_wave_chains import (
    chain_statistics,
    constant_kernel,
    convolution,
    log_flat_probability,
    matrix_kernel,
    rotating,
)
from scipy.integrate import simpson
from scipy.special import gammaln


def law_from_moments(p, chain, past_second, future_second):
    """Five unlabelled three-point orders; time orientation is retained."""
    stars = past_second + future_second
    return np.array(
        [
            1 - 3 * p + 3 * stars,
            3 * p + chain - 6 * stars,
            3 * future_second - chain,
            3 * past_second - chain,
            chain,
        ]
    )


def three_point_law(kernel, stats):
    left, right = kernel[:, 0], kernel[-1]
    c2left = convolution(kernel, left)
    c2right = convolution(kernel[::-1, ::-1].T, right[::-1])[::-1]
    alpha = np.divide(c2left, left, out=np.zeros_like(left), where=left > 0)
    beta = np.divide(c2right, right, out=np.zeros_like(right), where=right > 0)
    gamma = stats["C_r"][1] / stats["C_r"][0]
    weight = left * right / stats["C_r"][0]
    second = [
        float(simpson(weight * value**2, dx=1 / (len(left) - 1)) / (15 * gamma**3))
        for value in (alpha, beta)
    ]
    p, chain = stats["P_k"][1:3]
    law = law_from_moments(p, chain, *second)
    assert np.min(law) >= 0 and abs(sum(law) - 1) < 1e-12
    # C4 convolution and the independent star integrals must agree on E[q^2].
    assert abs(sum(second) + chain / 3 - stats["q_second"]) < 2e-10
    assert abs(law @ np.array([0, 1, 2, 2, 3]) - 3 * p) < 1e-12
    return {"past_second": second[0], "future_second": second[1], "law": law.tolist()}


def combinatorial_check():
    """Exhaust all 19 labelled strict posets on three vertices."""
    edges = [(i, j) for i in range(3) for j in range(3) if i != j]
    checked = 0
    for bits in itertools.product((0, 1), repeat=6):
        relation = np.zeros((3, 3), dtype=int)
        for edge, bit in zip(edges, bits, strict=True):
            relation[edge] = bit
        if np.any((relation @ relation > 0) & (relation == 0)):
            continue
        number = int(relation.sum())
        past = int(np.count_nonzero(relation.sum(axis=0) == 2))
        future = int(np.count_nonzero(relation.sum(axis=1) == 2))
        index = {0: 0, 1: 1, 3: 4}.get(number, 2 if future else 3)
        predicted = law_from_moments(number / 3, int(number == 3), past / 3, future / 3)
        assert np.max(np.abs(predicted - np.eye(5)[index])) < 1e-14
        checked += 1
    assert checked == 19
    return checked


def literature_normalization():
    rows = []
    for d in (2, 3, 4, 6, 10):
        for k in range(1, 7):
            n = 37
            falling = gammaln(n + 1) - gammaln(n - k + 1)
            # Aghili--Bombelli--Pilgrim Eq. (6): j=k+1 edges, k interior points.
            log_abp = (
                falling
                + (k - 1) * (gammaln(d + 1) - math.log(2))
                + gammaln(d / 2 + 1)
                + gammaln(d)
                - gammaln(k * d / 2 + 1)
                - gammaln((k + 1) * d / 2)
            )
            log_here = falling - gammaln(k + 1) + log_flat_probability(d, k)
            rows.append(abs(log_here - log_abp))
    error = float(max(rows))
    assert error < 1e-12
    return {"cases": len(rows), "max_log_difference": error, "fixed_N": 37}


def flat_ratio_asymptotics():
    def rate(k):
        return (k * math.log(k) + (k + 1) * math.log(k + 1)) / 2 - k * math.log(2)

    rows = []
    for k in range(3, 7):
        lam = rate(k) - (k - 1) * rate(2)
        log_prefactor = gammaln(k) + math.log(k * (k + 1) / 2) / 2
        log_prefactor -= (k - 1) * math.log(3) / 2
        residuals = []
        for d in (100, 1000):
            log_ratio = log_flat_probability(d, k) - (k - 1) * log_flat_probability(
                d, 2
            )
            residuals.append(float(log_ratio + lam * d - log_prefactor))
        assert lam > 0 and abs(residuals[-1]) < 1e-3
        assert abs(residuals[-1]) < abs(residuals[0]) / 5
        rows.append(
            {
                "k": k,
                "lambda_k": lam,
                "log_prefactor": log_prefactor,
                "log_asymptotic_residual_d100_d1000": residuals,
            }
        )
    return rows


def ratio_from_powers(stats, k):
    c = stats["C_r"]
    return (
        math.factorial(k)
        * 2**k
        * 3 ** (k - 1)
        / math.factorial(2 * k)
        * (c[k] * c[1] ** (k - 2) / c[2] ** (k - 1))
    )


def profile_checks():
    rows = []
    flat_law = np.array([265, 66, 9, 9, 1]) / 350
    for name, delta in (
        ("flat", None),
        ("rotating", None),
        ("caustic_1e-3", 1e-3),
        ("caustic_1e-6", 1e-6),
    ):
        refinements = []
        for panels in (1024, 2048):
            kernel = (
                matrix_kernel(rotating, panels)[0]
                if name == "rotating"
                else constant_kernel(math.pi - delta if delta else 0, panels)
            )
            stats = chain_statistics(kernel)
            law = three_point_law(kernel, stats)
            ratios = [ratio_from_powers(stats, k) for k in range(3, 7)]
            for k, ratio in enumerate(ratios, 3):
                assert (
                    abs(ratio / (stats["P_k"][k - 1] / stats["P_k"][1] ** (k - 1)) - 1)
                    < 1e-12
                )
            refinements.append(
                {
                    "panels": panels,
                    "P_k": stats["P_k"],
                    "R_3_through_R_6": ratios,
                    **law,
                }
            )
        last, previous = refinements[-1], refinements[-2]
        change = float(np.max(np.abs(np.array(last["P_k"]) / previous["P_k"] - 1)))
        assert change < 1e-7
        if name == "flat":
            assert np.max(np.abs(np.array(last["law"]) - flat_law)) < 1e-9
        # Actual three-point laws at two geometries, not a fractional-dimensional model.
        tv = float(np.abs(np.array(last["law"]) - flat_law).sum() / 2)
        antichain_bound = min(1, 3 * max(last["P_k"][1], 0.1))
        assert tv <= antichain_bound + 1e-9
        rows.append(
            {
                "profile": name,
                "delta": delta,
                "refinements": refinements,
                "last_relative_P_change": change,
                "TV_to_flat4_at_N3": tv,
                "antichain_TV_upper_bound": antichain_bound,
            }
        )
    near, nearer = rows[-2:]
    for k in range(3, 7):
        r1, r2 = (
            row["refinements"][-1]["R_3_through_R_6"][k - 3] for row in (near, nearer)
        )
        assert r2 > 0 and abs(r1 / r2 - 1) < 0.002
    return rows


if __name__ == "__main__":
    started = time.perf_counter()
    here = Path(__file__)
    sources = (here, here.with_name("check_plane_wave_chains.py"))
    result = {
        "kind": "Post-study deterministic normalization and corollary checks",
        "date": "2026-09-22",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256_lf": {
            p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
            for p in sources
        },
        "three_point_order": [
            "antichain",
            "one_relation",
            "one_minimum_fork",
            "one_maximum_fork",
            "chain",
        ],
        "labelled_posets_checked": combinatorial_check(),
        "ABP_finite_N_normalization": literature_normalization(),
        "flat_ratio_asymptotics": flat_ratio_asymptotics(),
        "profiles": profile_checks(),
        "numerical_scope": (
            "Floating-point refinement diagnostics, "
            "not certified enclosures or a sampling campaign."
        ),
    }
    result["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    destination = here.with_name("chain_limits_audit.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "written": str(destination),
                "all_checks_passed": True,
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
        )
    )
