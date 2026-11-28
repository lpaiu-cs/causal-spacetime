"""Deterministic control for the finite-component smoothing lemma."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

ROOT = Path(__file__).resolve().parents[2]
PROOF = ROOT / "docs/paper/paper_a/order_wave_energy/quantized_order_bound.md"
EDGES = np.array([[-1.0, 1, 0], [0, -1.0, 1], [1, 0, -1.0]])


def pattern_probabilities(eps):
    """Volumes in [0,1]^3 cut by all cyclic pair inequalities."""
    probabilities = []
    for bits in itertools.product((False, True), repeat=3):
        sign = np.where(bits, 1.0, -1.0)
        matrix = np.vstack((np.eye(3), -np.eye(3), EDGES * sign[:, None]))
        threshold = np.concatenate((np.zeros(3), -np.ones(3), -eps * sign))
        vertices = []
        for active in itertools.combinations(range(len(threshold)), 3):
            indices = np.array(active)
            face = matrix[indices]
            if abs(np.linalg.det(face)) < 0.5:
                continue
            point = np.linalg.solve(face, threshold[indices])
            if np.all(matrix @ point >= threshold - 2e-14):
                vertices.append(point)
        vertices = np.unique(np.round(vertices, 14), axis=0)
        if len(vertices) < 4 or np.linalg.matrix_rank(vertices - vertices[0]) < 3:
            probabilities.append(0.0)
        else:
            probabilities.append(float(ConvexHull(vertices).volume))
    return np.array(probabilities)


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Freeze source and proof; run from a clean worktree")
    if output.exists():
        raise FileExistsError(output)
    baseline = pattern_probabilities(0)
    rows = []
    for eps in (0.0, 0.001, 0.01, 0.1, 0.4):
        probabilities = pattern_probabilities(eps)
        expected = 4.5 * eps**2 - 4 * eps**3
        assert abs(probabilities.sum() - 1) < 2e-12
        assert abs(probabilities[-1] - expected) < 2e-12
        hellinger = float(np.sum((np.sqrt(probabilities) - np.sqrt(baseline)) ** 2))
        rows.append(
            {
                "eps": eps,
                "probabilities": probabilities.tolist(),
                "cycle_probability": float(probabilities[-1]),
                "cycle_polynomial": expected,
                "hellinger_squared": hellinger,
                "hellinger_squared_over_eps_squared": hellinger / eps**2
                if eps
                else None,
            }
        )
    record = {
        "scope": "unit-cube three-time polytope check; not wave detection power",
        "code_version": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
        .decode()
        .strip(),
        "files": {
            path.relative_to(ROOT).as_posix(): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in (Path(__file__), PROOF)
        },
        "rows": rows,
        "frozen_absolute_tolerance": 2e-12,
        "all_controls_passed": True,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    run(parser.parse_args().output)
