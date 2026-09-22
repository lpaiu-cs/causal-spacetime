"""Frozen independent coverage check for the blind fifth-degree inverse."""

import argparse
import hashlib
import json
import platform
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy
from intrinsic_weyl_coverage import selfcheck as check_confidence
from intrinsic_weyl_joint import selfcheck_jet
from intrinsic_weyl_joint_pilot import ROOT
from intrinsic_weyl_rotating_pilot import evaluate, generate
from rotating_plane_wave import selfcheck as check_generator
from scipy.stats import beta


def verdict(rows):
    """Failed generation, fitting and bias checks remain in the denominator."""
    if len(rows) != 256:
        raise ValueError("The frozen coverage gate requires all 256 cases")
    count = sum(
        "failure" not in row
        and row.get("bias_condition_met", False)
        and row.get("nominal_ball_contains_truth", False)
        for row in rows
    )
    interval = [
        float(beta.ppf(0.025, count, 257 - count)) if count else 0.0,
        float(beta.ppf(0.975, count + 1, 256 - count)) if count < 256 else 1.0,
    ]
    return {
        "valid_covered": count,
        "n": 256,
        "exact_95_two_sided_interval": interval,
        "coverage_gate": "SUPPORTED" if interval[0] >= 0.95 else "NOT_SUPPORTED",
    }


def selfcheck():
    valid = {"bias_condition_met": True, "nominal_ball_contains_truth": True}
    failed = valid | {"failure": "synthetic failed fit"}
    assert verdict([valid] * 251 + [failed] * 5)["coverage_gate"] == "SUPPORTED"
    assert verdict([valid] * 250 + [failed] * 6)["coverage_gate"] == "NOT_SUPPORTED"
    assert verdict([valid | {"bias_condition_met": False}] * 256)["valid_covered"] == 0
    return check_confidence()


def run_case(arguments):
    destination, case, amplitude, rotation, seed = arguments
    path = destination / f"case_{case}.npz"
    truth = destination / f"truth_{case}.npz"
    full = destination / f"case_{case}_full.json"
    started = time.perf_counter()
    row = {
        "case": case,
        "amplitude": amplitude,
        "rotation_rate": rotation,
        "seed": seed,
    }
    try:
        row["generation"] = generate(path, truth, amplitude, rotation, seed)
        result = evaluate(path, truth, model_order=5)
        full.write_text(json.dumps(result, indent=2) + "\n")
        row.update(result["evaluation"])
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        row["failure"] = f"{type(error).__name__}: {error}"
    row["seconds"] = time.perf_counter() - started
    row["artifact_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (path, truth, full, destination / f"case_{case}_selection_failed.npz")
        if p.exists()
    }
    (destination / f"case_{case}_result.json").write_text(
        json.dumps(row, indent=2) + "\n"
    )
    return row


def run(destination):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit protocol and code before independent data")
    if destination.exists():
        raise FileExistsError("Preserve all earlier coverage outputs")
    checked = {
        "gate": selfcheck(),
        "inverse": selfcheck_jet(),
        "generator": check_generator(),
    }
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    rng = np.random.default_rng(2026092950)
    population = list(
        zip(
            rng.choice([0.0, 0.25, 0.75], 256).tolist(),
            rng.choice([-0.5, 0.0, 0.5], 256).tolist(),
            strict=True,
        )
    )
    sources = (
        [Path(__file__)]
        + [
            Path(__file__).with_name(name)
            for name in (
                "intrinsic_weyl_rotating_pilot.py",
                "rotating_plane_wave.py",
                "intrinsic_weyl_joint.py",
                "intrinsic_weyl.py",
                "intrinsic_weyl_coverage.py",
                "intrinsic_weyl_joint_audit.py",
                "intrinsic_weyl_joint_pilot.py",
                "window_response.py",
                "p14_plane_wave.py",
            )
        ]
        + [
            ROOT / "docs/paper/paper_a/intrinsic_weyl" / name
            for name in (
                "jet_coverage_protocol.md",
                "rotating_generator.md",
            )
        ]
    )
    manifest = {
        "code_version": stamp,
        "selfcheck": checked,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "population_seed": 2026092950,
        "data_seed": 2026100100,
        "population": population,
        "selection_n": 4096,
        "measurement_n": 262144,
        "model_order": 5,
        "workers": 4,
        "source_sha256_lf": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }
    destination.mkdir(parents=True)
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    arguments = [
        (destination, i, a, w, 2026100100 + i) for i, (a, w) in enumerate(population)
    ]
    rows = []
    # Independent streams and case-specific files; no shared fitting state.
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row in pool.map(run_case, arguments):
            rows.append(row)
            print(
                json.dumps(
                    {
                        k: v
                        for k, v in row.items()
                        if k
                        not in ("artifact_sha256", "truth_coefficients_xy_positive")
                    }
                ),
                flush=True,
            )
    summary = (
        manifest
        | verdict(rows)
        | {
            "cases": rows,
            "scope": (
                "End-to-end coverage in the frozen rotating-wave mixture; "
                "not accuracy or a uniform regularity theorem"
            ),
        }
    )
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(verdict(rows)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/intrinsic_weyl_jet_coverage_v1"
    )
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    else:
        run(args.output)
