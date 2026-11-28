"""Post-result audit of the fixed 24-assignment high-count continuation."""

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEVELS = (262144, 1048576, 4194304)
CASES = (
    (0, 0),
    (0, 0.5),
    (0.25, 0),
    (0.75, 0),
    (0.25, 0.5),
    (0.75, 0.5),
    (0.25, -0.5),
    (0.75, -0.5),
)
FORBIDDEN = (
    "rotating_plane_wave",
    "window_response",
    "p14_plane_wave",
    "intrinsic_weyl_joint_audit",
    "intrinsic_weyl_joint_pilot",
    "intrinsic_weyl_rotating_pilot",
    "torch",
    "cupy",
    "jax",
    "tensorflow",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def pure_replay(folder):
    """Fresh process; only the predetermined case-7 small order archive opens."""
    import intrinsic_weyl_joint as joint
    import numpy as np

    archive = folder / "case_7_n262144.npz"
    original_load = np.load
    opened = []

    def guarded_load(path, *args, **kwargs):
        assert Path(path).resolve() == archive.resolve(), str(path)
        opened.append(str(Path(path).resolve()))
        return original_load(path, *args, **kwargs)

    np.load = guarded_load
    wall = time.perf_counter()
    assert not any(name in sys.modules for name in FORBIDDEN)
    fit = joint.fit_order_archive(archive, model_order=5, square_correction=True)
    assert not any(name in sys.modules for name in FORBIDDEN)
    original = read(folder / "case_7_n262144_pure.json")["fit"]
    assert fit == original, "Pure replay differs from the saved fit"
    return {
        "case": 7,
        "measurement_n": 262144,
        "all_fit_fields_exact": True,
        "max_absolute_difference": 0.0,
        "opened_archives": opened,
        "oracle_and_gpu_modules_loaded": [],
        "fresh_events": 0,
        "inverse_module": joint.__file__,
        "wall_seconds": time.perf_counter() - wall,
    }


def audit(folder, pilot, output):
    import numpy as np
    import scipy
    from check_intrinsic_weyl_coverage import check_nested_evaluation
    from intrinsic_weyl import WEYL
    from rotating_plane_wave import RotatingWave
    from scipy.stats import chi2

    wall = time.perf_counter()
    assert not output.exists()
    status = read(folder / "coordinator_status.json")
    assert status["state"] == "COMPLETED" and status["assignments"] == 24
    manifest = read(folder / "manifest.json")
    continuation = read(folder / "continuation.json")
    summary = read(folder / "summary.json")
    assert all(summary[k] == v for k, v in manifest.items() if k != "scope")
    assert summary["continuation"] == continuation
    assert manifest["cases"] == [list(c) for c in CASES]
    assert manifest["measurement_levels"] == list(LEVELS)
    assert manifest["measurement_seed_base"] == 2026103000
    assert manifest["chunk"] == 4096 and manifest["fresh_selection_events"] == 0
    assert manifest["runtime"] == {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    assert {k: os.environ.get(k) for k in manifest["thread_environment"]} == manifest[
        "thread_environment"
    ]
    assert continuation["thread_environment"] == manifest["thread_environment"]
    assert continuation["gpu_used"] is False and continuation["seed"] == 2026103007
    assert continuation["completed_cases_reused"] == list(range(7))
    assert continuation["remaining_case"] == 7
    original = Path(continuation["original_directory"])
    evidence_path = original.with_name(original.name + ".interruption_evidence.json")
    assert sha(evidence_path) == continuation["evidence_sha256"]
    evidence = read(evidence_path)
    preserved = evidence["original_file_sha256"]
    assert preserved == continuation["preserved_sha256"] and len(preserved) == 132
    assert set(p.name for p in original.iterdir()) == set(preserved)
    for name, digest in preserved.items():
        assert sha(original / name) == sha(folder / name) == digest, name
    for name, digest in manifest["source_sha256"].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in manifest["pilot_input_sha256"].items():
        assert sha(pilot / name) == digest, name
    for name, digest in continuation["coordinator_source_sha256"].items():
        assert sha(ROOT / name) == digest, name
    for stamp in (manifest["code_version"], continuation["coordinator_code_version"]):
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=ROOT, check=True
        )
    preflight = read(folder / "continuation_replay_comparison.json")
    assert preflight["within_frozen_tolerance"] and preflight["all_full_fields_exact"]
    assert preflight["differences"] == {} and preflight["fresh_events"] == 0
    expected = {(i, n) for i in range(8) for n in LEVELS}
    assert summary["completed_assignments"] == len(summary["results"]) == 24
    assert {(r["case"], r["measurement_n"]) for r in summary["results"]} == expected
    assert len(list(folder.glob("case_*_result.json"))) == 24
    generation = [read(folder / f"case_{i}_generation.json") for i in range(8)]
    assert generation == summary["generation"]
    assert all(
        g["sampled_events"] == g["classified_events"] == LEVELS[-1] for g in generation
    )
    assert generation[7]["saved_prefixes_reconstructed_bit_exact"] == list(LEVELS[:2])
    assert summary["unique_archived_measurement_events"] == 33554432
    assert (
        summary["sample_operations_lower_bound"]
        == 33554432 + evidence["last_saved_status"]["events"]
    )

    epsilon = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        epsilon[p] = (-1) ** sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
    electric = WEYL[:, 0, 1:, 0, 1:]
    magnetic = 0.5 * np.einsum("ikl,aklj->aij", epsilon, WEYL[:, 1:, 1:, 0, 1:])
    physical_map = np.concatenate(
        (electric.reshape(10, 9), magnetic.reshape(10, 9)), axis=1
    ).T

    def eb(coefficients):
        return tuple(
            np.einsum("a,aij->ij", coefficients, basis)
            for basis in (electric, magnetic)
        )

    def canonical(fit):
        e, b = eb(fit["coefficients"])
        points = np.asarray(fit["coordinates_inferred_from_counts"])
        d = np.diag([np.sign(points[1, 1]), np.sign(points[2, 2]), 1])
        return d @ e @ d, np.linalg.det(d) * d @ b @ d

    def distance(left, right):
        e, b = left
        et, bt = right
        d = np.diag([1.0, 1.0, -1.0])
        return min(
            float(np.sqrt(np.sum((e - a) ** 2) + np.sum((b - c) ** 2)))
            for a, c in (right, (d @ et @ d, -d @ bt @ d))
        )

    rows = []
    allowed = {
        "selection_order",
        "selection_design",
        "selected_ids",
        "after",
        "before",
        "measurement_n",
    }
    for case, (amplitude, rotation) in enumerate(CASES):
        with np.load(pilot / f"case_{case}_nested.npz", allow_pickle=False) as old:
            selected = {
                k: old[k]
                for k in ("selection_order", "selection_design", "selected_ids")
            }
        with np.load(
            folder / f"case_{case}_n{LEVELS[-1]}.npz", allow_pickle=False
        ) as large:
            assert (
                set(large.files) == allowed
                and int(large["measurement_n"]) == LEVELS[-1]
            )
            assert all(np.array_equal(large[k], v) for k, v in selected.items())
            after, before, ids = (
                large["after"],
                large["before"],
                selected["selected_ids"],
            )
            assert after.shape == before.shape == (len(ids), LEVELS[-1] // 8)
            assert after.dtype == before.dtype == np.uint8
            assert np.all(after[0] == 255) and np.all(before[-1] == 255)
            assert (
                not np.any(before[0])
                and not np.any(after[-1])
                and not np.any(after & before)
            )
            for n in LEVELS[:-1]:
                with np.load(
                    folder / f"case_{case}_n{n}.npz", allow_pickle=False
                ) as small:
                    assert (
                        set(small.files) == allowed and int(small["measurement_n"]) == n
                    )
                    assert all(np.array_equal(small[k], v) for k, v in selected.items())
                    assert np.array_equal(small["after"], after[:, : n // 8])
                    assert np.array_equal(small["before"], before[:, : n // 8])
            order = selected["selection_order"][ids[:, None], ids]
            for a, b in np.argwhere(order):
                assert not np.any(after[b] & ~after[a]) and not np.any(
                    before[a] & ~before[b]
                )
            for a, b in np.argwhere(~order):
                assert not np.any(after[a] & before[b])
        with np.load(pilot / f"truth_{case}_nested.npz", allow_pickle=False) as truth:
            assert (float(truth["amplitude"]), float(truth["rotation_rate"])) == (
                amplitude,
                rotation,
            )
            normal = RotatingWave(amplitude, rotation).normal_coordinates(
                truth["points"]
            )
        x = normal[1, 1:].copy()
        x /= np.linalg.norm(x)
        y = normal[2, 1:] - normal[2, 1:] @ x * x
        y /= np.linalg.norm(y)
        q = np.column_stack([x, y, np.cross(x, y)])
        e0 = np.diag([-amplitude, amplitude, 0.0])
        b0 = np.array([[0.0, -amplitude, 0.0], [-amplitude, 0.0, 0.0], [0.0, 0.0, 0.0]])
        truth_eb = q.T @ e0 @ q, q.T @ b0 @ q
        for n in LEVELS:
            path = folder / f"case_{case}_n{n}.json"
            record = read(path)
            full = record["full"]
            fit, evaluation = full["fit"], full["evaluation"]
            pure = read(folder / f"case_{case}_n{n}_pure.json")
            row = read(folder / f"case_{case}_n{n}_result.json")
            assert row == next(
                r
                for r in summary["results"]
                if (r["case"], r["measurement_n"]) == (case, n)
            )
            assert "failure" not in row
            assert pure["fit"] == fit and pure["oracle_modules_before_return"] == []
            assert row["full_record_sha256"] == sha(path)
            assert (
                row["archive_sha256"]
                == record["archive_sha256"]
                == pure["archive_sha256"]
                == fit["order_archive_sha256"]
                == sha(folder / f"case_{case}_n{n}.npz")
            )
            assert (
                row["archive_bytes"]
                == (folder / f"case_{case}_n{n}.npz").stat().st_size
            )
            assert (
                row["seed"] == 2026103000 + case
                and row["measurement_n"] == fit["measurement_n"] == n
            )
            assert row["amplitude"] == amplitude and row["rotation_rate"] == rotation
            assert (
                record["thread_environment"]
                == row["thread_environment"]
                == manifest["thread_environment"]
            )
            assert all(row[k] == v for k, v in evaluation.items())
            saved_truth = eb(evaluation["truth_coefficients_xy_positive"])
            assert (
                max(
                    np.max(abs(u - v))
                    for u, v in zip(truth_eb, saved_truth, strict=True)
                )
                < 1e-12
            )
            error, bias = (
                distance(canonical(fit), truth_eb),
                distance(canonical(full["continuum_fit"]), truth_eb),
            )
            norm = float(np.linalg.norm(physical_map @ np.asarray(fit["coefficients"])))
            variance = float(
                np.linalg.eigvalsh(
                    physical_map
                    @ np.asarray(fit["nominal_component_covariance"])
                    @ physical_map.T
                )[-1]
            )
            radius = (math.sqrt(chi2.ppf(0.95, 10) * variance) + 0.02 * norm) / 0.98
            for key, value in (
                ("full_tidal_error_mod_reflection", error),
                ("continuum_bias_mod_reflection", bias),
                ("confidence_radius", radius),
                ("largest_tidal_variance", variance),
            ):
                assert abs(evaluation[key] - value) < 1e-12, (case, n, key)
            valid = check_nested_evaluation(full)
            assert valid == (bias <= 0.04 * amplitude + 1e-10 and error <= radius)
            assert row["variance_price_V"] == n * row["largest_tidal_variance"]
            rows.append(
                {
                    "case": case,
                    "A": amplitude,
                    "w": rotation,
                    "N": n,
                    "error": error,
                    "radius": radius,
                    "bias": bias,
                    "relative_error": error / (2 * amplitude) if amplitude else None,
                    "relative_radius": radius / (2 * amplitude) if amplitude else None,
                    "relative_bias": bias / (2 * amplitude) if amplitude else None,
                    "variance_times_N": n * variance,
                    "valid_and_covered": bool(valid),
                    "triangle_upper_relative_error_plus_radius": (error + radius)
                    / (2 * amplitude)
                    if amplitude
                    else None,
                    "full_file_sha256": sha(path),
                }
            )
        print(
            json.dumps({"case_audited": case, "results_checked": len(rows)}), flush=True
        )
    replay = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--folder",
            str(folder),
            "--replay",
        ],
        check=True,
        capture_output=True,
        text=True,
        creationflags=subprocess.CREATE_NO_WINDOW,
    )
    replay_result = json.loads(replay.stdout)
    result = {
        "completed_assignments_audited": len(rows),
        "source_files_checked": len(manifest["source_sha256"]),
        "old_pilot_files_checked": len(manifest["pilot_input_sha256"]),
        "preserved_original_and_copied_files_checked": len(preserved),
        "all_24_input_schemas_selection_and_nested_bits_checked": True,
        "all_8_landmark_event_transitivity_checked": True,
        "summary_sha256": sha(folder / "summary.json"),
        "checker_sha256": sha(Path(__file__)),
        "post_completion_pure_replay": replay_result,
        "rows": rows,
        "gpu_used": False,
        "fresh_events": 0,
        "audit_wall_seconds": time.perf_counter() - wall,
        "scope": (
            "Saved-fit metrics independently recomputed via E/B tensors, covariance "
            "and reflection parity; established truth normal-coordinate map reused. "
            "One separate pure inverse replay. Does not regenerate all order "
            "streams or establish powered accuracy."
        ),
    }
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "output": str(output),
                "sha256": sha(output),
                "rows": len(rows),
                "pure_replay_exact": True,
                "wall_seconds": result["audit_wall_seconds"],
            }
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=Path, required=True)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    assert os.environ.get("CUDA_VISIBLE_DEVICES") == "-1"
    if args.replay:
        print(json.dumps(pure_replay(args.folder)))
    else:
        if args.pilot is None or args.output is None:
            parser.error("pilot and output are required for the full audit")
        audit(args.folder, args.pilot, args.output)
