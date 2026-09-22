"""Deterministic same-landmark continuum bias scaling, not finite-order recovery."""

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from intrinsic_weyl_joint import (
    confidence_ball,
    fit_joint,
    gauge,
    read_order_clocks,
    tidal_distance,
)
from intrinsic_weyl_joint_audit import aligned_tensor_error
from intrinsic_weyl_joint_pilot import ROOT
from rotating_plane_wave import RotatingWave


def run():
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit the deterministic diagnostic before running it")
    output = ROOT / "outputs/intrinsic_weyl_jet_scaling_v1.json"
    if output.exists():
        raise FileExistsError("Preserve the earlier deterministic result")
    base = ROOT / "outputs/intrinsic_weyl_jet_pilot_v1"
    order, truth = base / "case_7.npz", base / "truth_7.npz"
    counts, mask, covariance, _, _ = read_order_clocks(order)
    with np.load(truth, allow_pickle=False) as source:
        points = source["points"]
        amplitude, rotation = float(source["amplitude"]), float(source["rotation_rate"])
    assert (amplitude, rotation) == (0.75, -0.5)
    rows = []
    for scale in (1.0, 0.5, 0.25, 0.125):
        wave = RotatingWave(amplitude * scale**2, rotation * scale)
        normal = wave.normal_coordinates(points)
        expected = aligned_tensor_error(
            {
                "coefficients": np.zeros(10),
                "coordinates_inferred_from_counts": gauge(normal),
            },
            normal,
            wave.amplitude,
        )["aligned_truth_coefficients"]
        volumes = wave.volumes(points)
        assert np.all(volumes[mask] > 0)
        refinement = float(
            np.max(abs(wave.volumes(points, 32)[mask] / volumes[mask] - 1))
        )
        assert refinement < 1e-10
        row = {"scale": scale, "quadrature_relative_change": refinement}
        for degree in (4, 5):
            fitted = fit_joint(
                counts[0, -1] * volumes / volumes[0, -1],
                clock_mask=mask,
                clock_covariance=covariance,
                model_order=degree,
            )
            inferred = confidence_ball(fitted)
            row[str(degree)] = {
                "relative_continuum_bias": tidal_distance(
                    inferred["coefficients_xy_positive"], expected
                )
                / (2 * wave.amplitude),
                "fit": fitted,
            }
        rows.append(row)
        print(
            json.dumps(
                {
                    "scale": scale,
                    **{str(k): row[str(k)]["relative_continuum_bias"] for k in (4, 5)},
                }
            ),
            flush=True,
        )
    sources = [
        Path(__file__),
        Path(__file__).with_name("intrinsic_weyl_joint.py"),
        Path(__file__).with_name("intrinsic_weyl.py"),
        Path(__file__).with_name("intrinsic_weyl_joint_audit.py"),
        Path(__file__).with_name("rotating_plane_wave.py"),
        ROOT / "docs/paper/paper_a/intrinsic_weyl/local_consistency_scope.md",
    ]
    result = {
        "code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256_lf": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
        "input_sha256": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (order, truth)
        },
        "scope": (
            "Fixed existing landmarks and weights, continuum-only oracle audit; "
            "no new random orders, accuracy or coverage verdict"
        ),
        "cases": rows,
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Written {output}", flush=True)


if __name__ == "__main__":
    run()
