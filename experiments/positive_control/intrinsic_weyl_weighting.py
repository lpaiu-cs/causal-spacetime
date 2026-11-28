"""Independent noise-cost validation for GLS; no curved recovery claim."""

import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from intrinsic_weyl_noise import (
    CHUNK,
    PER_SEED,
    ROOT,
    SEEDS,
    linearized_fit,
    memberships,
    price,
)


def run():
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True
    )
    if dirty.strip():
        raise RuntimeError("Commit this diagnostic's protocol and implementation first")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    clocks, p, ordinary, tau2, _, design, response = linearized_fit()
    n = PER_SEED * len(SEEDS)
    total = np.zeros(len(design))
    second = np.zeros((len(design), len(design)))
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        for _ in range(PER_SEED // CHUNK):
            values = (memberships(rng, clocks) - p) @ response.T
            total += values.sum(axis=0)
            second += values.T @ values
    covariance = (second - np.outer(total, total) / n) / (n - 1)
    scales = np.sqrt(np.diag(covariance))
    standardized = covariance / scales[:, None] / scales[None, :]
    eigenvalues, vectors = np.linalg.eigh(standardized)
    keep = eigenvalues > eigenvalues[-1] * 1e-9
    whiten = (vectors[:, keep] / np.sqrt(eigenvalues[keep])).T / scales
    whitened_design = whiten @ design
    assert np.linalg.matrix_rank(whitened_design) == 10
    operator = np.linalg.pinv(whitened_design) @ whiten
    assert np.max(abs(operator @ design - np.eye(10))) < 1e-7
    assert np.max(abs(np.linalg.pinv(design) @ design - np.eye(10))) < 1e-7
    weighted = operator @ response
    assert np.linalg.norm(weighted @ p) < 1e-8
    print(
        json.dumps(
            {
                "training_complete": True,
                "covariance_rank": int(keep.sum()),
                "residual_rows": len(design),
            }
        ),
        flush=True,
    )
    rows = []
    for seed in range(2026092210, 2026092214):
        start = time.perf_counter()
        rng = np.random.default_rng(seed)
        summed = np.zeros((2, 10))
        moments = np.zeros((2, 10, 10))
        both = np.vstack((ordinary, weighted))
        for _ in range(PER_SEED // CHUNK):
            values = (memberships(rng, clocks) - p) @ both.T
            for index in range(2):
                value = values[:, index * 10 : (index + 1) * 10]
                summed[index] += value.sum(axis=0)
                moments[index] += value.T @ value
        estimates = [
            (moments[i] - np.outer(summed[i], summed[i]) / PER_SEED) / (PER_SEED - 1)
            for i in range(2)
        ]
        prices = [price(value, tau2) for value in estimates]
        price_key = "optimistic_single_projection_n_95ci_90power_25pct"
        ratios = [
            a[price_key] / b[price_key]
            for a, b in zip(prices[0], prices[1], strict=True)
        ]
        row = {
            "seed": seed,
            "n": PER_SEED,
            "seconds": time.perf_counter() - start,
            "ordinary": prices[0],
            "weighted": prices[1],
            "ordinary_to_weighted_price_ratio": ratios,
            "influence_covariances": [value.tolist() for value in estimates],
        }
        rows.append(row)
        print(
            json.dumps({k: v for k, v in row.items() if k != "influence_covariances"}),
            flush=True,
        )
    sources = [
        Path(__file__),
        Path(__file__).with_name("intrinsic_weyl.py"),
        Path(__file__).with_name("intrinsic_weyl_development.py"),
        Path(__file__).with_name("intrinsic_weyl_noise.py"),
        ROOT / "src/causal_spacetime_lab/sprinkling.py",
        ROOT / "docs/paper/paper_a/intrinsic_weyl/noise_design.md",
    ]
    return {
        "kind": "Independent validation of linearized noise price; not curved recovery",
        "code_version": stamp,
        "training_seeds": list(SEEDS),
        "training_events": n,
        "standardized_covariance_eigenvalues": eigenvalues.tolist(),
        "retained_covariance_rank": int(keep.sum()),
        "validation": rows,
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }


if __name__ == "__main__":
    result = run()
    destination = (
        ROOT / "docs/paper/paper_a/intrinsic_weyl/weighted_noise_diagnostic.json"
    )
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Written {destination}")
