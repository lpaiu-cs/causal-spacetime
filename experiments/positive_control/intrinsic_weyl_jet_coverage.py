"""Frozen independent coverage check for the blind fifth-degree inverse."""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import numpy as np
import scipy
from intrinsic_weyl_coverage import selfcheck as check_confidence
from intrinsic_weyl_joint import selfcheck_jet
from intrinsic_weyl_joint_pilot import ROOT
from intrinsic_weyl_rotating_pilot import evaluate, generate
from rotating_plane_wave import selfcheck as check_generator
from scipy.stats import beta, binom


def verdict(rows, expected_n=256):
    """Failed generation, fitting and bias checks remain in the denominator."""
    if len(rows) != expected_n:
        raise ValueError(f"The frozen coverage gate requires all {expected_n} cases")
    count = sum(
        "failure" not in row
        and row.get("bias_condition_met", False)
        and row.get("nominal_ball_contains_truth", False)
        for row in rows
    )
    interval = [
        float(beta.ppf(0.025, count, expected_n + 1 - count)) if count else 0.0,
        float(beta.ppf(0.975, count + 1, expected_n - count))
        if count < expected_n
        else 1.0,
    ]
    return {
        "valid_covered": count,
        "n": expected_n,
        "exact_95_two_sided_interval": interval,
        "coverage_gate": "SUPPORTED" if interval[0] >= 0.95 else "NOT_SUPPORTED",
    }


def selfcheck():
    valid = {"bias_condition_met": True, "nominal_ball_contains_truth": True}
    failed = valid | {"failure": "synthetic failed fit"}
    assert verdict([valid] * 251 + [failed] * 5)["coverage_gate"] == "SUPPORTED"
    assert verdict([valid] * 250 + [failed] * 6)["coverage_gate"] == "NOT_SUPPORTED"
    assert verdict([valid | {"bias_condition_met": False}] * 256)["valid_covered"] == 0
    assert verdict([valid] * 497 + [failed] * 15, 512)["coverage_gate"] == "SUPPORTED"
    assert (
        verdict([valid] * 496 + [failed] * 16, 512)["coverage_gate"] == "NOT_SUPPORTED"
    )
    return check_confidence()


def run_case(arguments, nested=False):
    destination, case, amplitude, rotation, seed = arguments
    suffix = "_nested" if nested else ""
    path = destination / f"case_{case}{suffix}.npz"
    truth = destination / f"truth_{case}{suffix}.npz"
    full = destination / f"case_{case}_full.json"
    started = time.perf_counter()
    row = {
        "case": case,
        "amplitude": amplitude,
        "rotation_rate": rotation,
        "seed": seed,
    }
    try:
        if nested:
            from intrinsic_weyl_design_pilot import generate as generate_nested

            arms = {"nested": {}}
            row["generation"] = generate_nested(
                destination, case, amplitude, rotation, seed, arms, ("nested",)
            )
            if "failure" in arms["nested"]:
                raise ValueError(arms["nested"]["failure"])
            row["generation"]["landmarks"] = arms["nested"]["landmarks"]
        else:
            row["generation"] = generate(path, truth, amplitude, rotation, seed)
        result = evaluate(path, truth, model_order=5, square_correction=nested)
        full.write_text(json.dumps(result, indent=2) + "\n")
        row.update(result["evaluation"])
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        row["failure"] = f"{type(error).__name__}: {error}"
    row["seconds"] = time.perf_counter() - started
    row["artifact_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (
            path,
            truth,
            full,
            destination / f"case_{case}{suffix}_selection_failed.npz",
        )
        if p.exists()
    }
    (destination / f"case_{case}_result.json").write_text(
        json.dumps(row, indent=2) + "\n"
    )
    return row


def run(destination, *, nested=False):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit protocol and code before independent data")
    if destination.exists():
        raise FileExistsError("Preserve all earlier coverage outputs")
    n, population_seed, data_seed, protocol = (
        (512, 2026100950, 2026101000, "nested_coverage_protocol.md")
        if nested
        else (256, 2026092950, 2026100100, "jet_coverage_protocol.md")
    )
    checked = {
        "gate": selfcheck(),
        "inverse": selfcheck_jet(),
        "generator": check_generator(),
    }
    if nested:
        from intrinsic_weyl_design_pilot import selfcheck as check_design

        checked["nested_design"] = check_design(
            ("connected", "nested"), model_order=5, square_correction=True
        )
        checked["nested_512_gate"] = {
            "required_successes": 497,
            "power_at_098": float(binom.sf(496, 512, 0.98)),
        }
        assert checked["nested_512_gate"]["power_at_098"] >= 0.9
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    rng = np.random.default_rng(population_seed)
    population = list(
        zip(
            rng.choice([0.0, 0.25, 0.75], n).tolist(),
            rng.choice([-0.5, 0.0, 0.5], n).tolist(),
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
                protocol,
                "rotating_generator.md",
            )
        ]
    )
    if nested:
        sources.extend(
            [
                Path(__file__).with_name("intrinsic_weyl_design_pilot.py"),
                Path(__file__).with_name("intrinsic_weyl_square.py"),
                ROOT / "docs/paper/paper_a/intrinsic_weyl/square_correction.md",
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
        "population_seed": population_seed,
        "data_seed": data_seed,
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
    if nested:
        manifest.update(
            selection_design="nested",
            curvature_square_correction=True,
            worker_thread_environment={
                name: os.environ.get(name)
                for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
            },
            thread_settings_changed=False,
        )
    destination.mkdir(parents=True)
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    arguments = [
        (destination, i, a, w, data_seed + i) for i, (a, w) in enumerate(population)
    ]
    rows = []
    # Independent streams and case-specific files; no shared fitting state.
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row in pool.map(run_case, arguments, repeat(nested)):
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
        | verdict(rows, n)
        | {
            "cases": rows,
            "scope": (
                "End-to-end coverage in the frozen rotating-wave mixture; "
                "not accuracy or a uniform regularity theorem"
            ),
        }
    )
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(verdict(rows, n)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--nested", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    else:
        output = args.output or ROOT / (
            "outputs/intrinsic_weyl_nested_coverage_v1"
            if args.nested
            else "outputs/intrinsic_weyl_jet_coverage_v1"
        )
        run(output, nested=args.nested)
