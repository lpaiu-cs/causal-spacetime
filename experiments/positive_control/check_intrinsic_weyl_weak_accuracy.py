"""Independent post-completion audit of the eight frozen precision assignments."""

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
CASES = (
    (0, 0, 0, 4000, 25108480),
    (0, 0.5, 1, 4000, 25108480),
    (0.25, 0, 2, 5500, 49311744),
    (0.25, 0.5, 4, 5700, 51101696),
    (0.25, -0.5, 6, 2900, 26001408),
    (0.25, -0.5, None, 8000, 71720960),
    (0.25, 0, None, 8000, 71720960),
    (0.25, 0.5, None, 8000, 71720960),
)
FORBIDDEN = (
    "rotating_plane_wave",
    "rotating_plane_wave_gpu",
    "window_response",
    "p14_plane_wave",
    "intrinsic_weyl_joint_audit",
    "intrinsic_weyl_joint_pilot",
    "intrinsic_weyl_rotating_pilot",
    "intrinsic_weyl_weak_accuracy",
)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def quantile10(probability):
    """Invert the even-df chi-square CDF without scipy.stats or runner code."""
    lo, hi = 0.0, 128.0
    for _ in range(100):
        mid = (lo + hi) / 2
        survival = math.exp(-mid / 2) * sum(
            (mid / 2) ** j / math.factorial(j) for j in range(5)
        )
        if survival > 1 - probability:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def conclusion(rows):
    if len(rows) != 8 or {row["case"] for row in rows} != set(range(8)):
        raise ValueError("Every distinct frozen assignment is required")
    return (
        "SUPPORTED"
        if all(not row.get("failure") and row.get("passed", False) for row in rows)
        else "NOT_SUPPORTED"
    )


def selfcheck():
    from scipy.stats import chi2

    for p in (0.95, 0.99375, 0.9875):
        assert abs(quantile10(p) - chi2.ppf(p, 10)) < 1e-12
    factor = (math.sqrt(quantile10(0.99375)) + math.sqrt(quantile10(0.9875))) ** 2
    assert abs(factor - 94.18944697204289) < 1e-12
    for a, _, _, cap, n in CASES:
        margin = 0.205 * 2 * a if a else 0.98 * 0.125
        assert 4096 * math.ceil(cap * factor / (margin - 2e-10) ** 2 / 4096) == n
    rows = [dict(case=i, passed=True) for i in range(8)]
    assert conclusion(rows) == "SUPPORTED"
    for broken in (rows[:-1], rows[:-1] + [rows[0]]):
        try:
            conclusion(broken)
        except ValueError:
            pass
        else:
            raise AssertionError("Missing or duplicated assignment was accepted")
    rows[5]["failure"] = "Failure cannot be erased by a favorable flag"
    assert conclusion(rows) == "NOT_SUPPORTED"
    return {
        "independent_quantiles_and_power_counts": True,
        "missing_duplicate_and_failed_cases_rejected": True,
    }


def replay(folder, case):
    """Fresh process: only the chosen order archive opens until the pure return."""
    import intrinsic_weyl_joint as joint
    import numpy as np
    from intrinsic_weyl_gpu_counts import clock_second_moment

    archive = folder / f"case_{case}.npz"
    original_load, opened = np.load, []

    def guard(path, *args, **kwargs):
        assert Path(path).resolve() == archive.resolve(), str(path)
        opened.append(str(archive.resolve()))
        return original_load(path, *args, **kwargs)

    joint._clock_second_moment = clock_second_moment
    np.load = guard
    assert not any(name in sys.modules for name in FORBIDDEN)
    wall = time.perf_counter()
    fit = joint.fit_order_archive(archive, model_order=5, square_correction=True)
    assert not any(name in sys.modules for name in FORBIDDEN)
    assert fit == read(folder / f"case_{case}_full_pure.json")["fit"]
    return {
        "case": case,
        "measurement_n": fit["measurement_n"],
        "all_fit_fields_exact": True,
        "opened_archives": opened,
        "oracle_modules_loaded": [],
        "new_scientific_events": 0,
        "wall_seconds": time.perf_counter() - wall,
    }


def audit(folder, pilot, output):
    import numpy as np
    import scipy
    from check_intrinsic_weyl_coverage import check_nested_evaluation
    from intrinsic_weyl import WEYL
    from intrinsic_weyl_joint import selection
    from rotating_plane_wave import RotatingWave

    wall = time.perf_counter()
    assert not output.exists()
    checked = selfcheck()
    manifest, summary, status = (
        read(folder / f"{n}.json") for n in ("manifest", "summary", "status")
    )
    assert status["state"] == "COMPLETED" and status["assignments"] == 8
    assert summary["completed_assignments"] == len(summary["results"]) == 8
    assert all(summary[k] == v for k, v in manifest.items())
    assert summary["audit_status"] == "PENDING_INDEPENDENT_POST_COMPLETION_AUDIT"
    assert len(list(folder.glob("case_*_result.json"))) == 8
    assert {r["case"] for r in summary["results"]} == set(range(8))
    assert manifest["alpha_family"] == 0.05 and manifest["beta_joint_power"] == 0.1
    assert manifest["relative_bias_allowance"] == 0.02
    assert manifest["absolute_bias_slack"] == 1e-10
    assert manifest["weak_relative_accuracy_target"] == 0.25
    assert manifest["flat_absolute_accuracy_target"] == 0.125
    assert (
        manifest["sampling_chunk"] == 4096 and manifest["classification_batch"] == 32768
    )
    assert manifest["maximum_measurement_events"] == 391794688
    assert manifest["maximum_new_selection_events"] == 12288
    assert manifest["supplied_new_tips"] == 6 and manifest["gpu_used"]
    for key, value in dict(
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__
    ).items():
        assert manifest["runtime"][key] == value
    assert {k: os.environ.get(k) for k in manifest["thread_environment"]} == manifest[
        "thread_environment"
    ]
    for name, digest in manifest["source_sha256"].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in manifest["pilot_input_sha256"].items():
        assert sha(pilot / name) == digest, name
    for name, digest in manifest["planning_evidence_sha256"].items():
        assert sha(Path(name)) == digest, name
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", manifest["code_version"], "HEAD"],
        cwd=ROOT,
        check=True,
    )
    assert manifest["code_version"] == "0ff98023c90ebbdfbf74b5577e319b95eb0f58c9"
    preflight = read(folder / "preflight_comparison.json")
    assert preflight["all_full_fields_exact"] and preflight["within_frozen_tolerance"]
    assert preflight["differences"] == {} and preflight["scientific_fresh_events"] == 0
    reference = next(
        Path(name).parent
        for name in manifest["planning_evidence_sha256"]
        if "high_count_v1_continued" in name
    )
    assert read(folder / "preflight_worker.json")["full"] == read(
        reference / "case_4_n262144.json"
    )["full"]

    epsilon = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        epsilon[p] = (-1) ** sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
    electric = WEYL[:, 0, 1:, 0, 1:]
    magnetic = 0.5 * np.einsum("ikl,aklj->aij", epsilon, WEYL[:, 1:, 1:, 0, 1:])
    physical = np.concatenate(
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

    rows, hashes = [], 0
    allowed = {
        "selection_order",
        "selection_design",
        "selected_ids",
        "after",
        "before",
        "measurement_n",
    }
    for i, (amplitude, rotation, retained, cap, n) in enumerate(CASES):
        assignment = dict(
            case=i,
            amplitude=amplitude,
            rotation_rate=rotation,
            retained_pilot_case=retained,
            variance_cap=cap,
            measurement_n=n,
            seed=2026112000 + i,
        )
        assert assignment == manifest["cases"][i]
        row = read(folder / f"case_{i}_result.json")
        assert row == next(r for r in summary["results"] if r["case"] == i)
        assert all(row[k] == v for k, v in assignment.items())
        generation = read(folder / f"case_{i}_generation.json")
        assert generation == row["generation"]
        paths = set(folder.glob(f"case_{i}_*")) | {
            folder / f"case_{i}.npz",
            folder / f"truth_{i}.npz",
        }
        expected_files = {
            p.name for p in paths if p.is_file() and not p.name.endswith("_result.json")
        }
        assert set(row["artifact_sha256"]) == expected_files
        for name, digest in row["artifact_sha256"].items():
            assert sha(folder / name) == digest, name
            hashes += 1
        archive = folder / f"case_{i}.npz"
        if not archive.exists():
            assert row.get("failure")
            rows.append(dict(case=i, failure=row["failure"], passed=False))
            continue
        assert generation["sampled_events"] == generation["classified_events"] == n
        assert generation["selection_events_drawn"] == (4096 if retained is None else 0)
        assert row["archive_bytes"] == archive.stat().st_size
        with np.load(archive, allow_pickle=False) as saved:
            assert set(saved.files) == allowed and int(saved["measurement_n"]) == n
            assert str(saved["selection_design"]) == "nested"
            order, ids, after, before = (
                saved[k] for k in ("selection_order", "selected_ids", "after", "before")
            )
        assert order.shape == (4098, 4098) and order.dtype == np.bool_
        expected_ids, mask = selection(order, design="nested")
        assert np.array_equal(ids, expected_ids)
        assert after.shape == before.shape == (len(ids), n // 8)
        assert after.dtype == before.dtype == np.uint8
        assert np.all(after[0] == 255) and np.all(before[-1] == 255)
        assert (
            not np.any(before[0])
            and not np.any(after[-1])
            and not np.any(after & before)
        )
        assert generation["landmarks_including_tips"] == len(ids)
        assert generation["directional_predicates_classified"] == 2 * len(ids) * n
        restricted = order[ids[:, None], ids]
        for a, b in np.argwhere(restricted):
            assert not np.any(after[b] & ~after[a]) and not np.any(
                before[a] & ~before[b]
            )
        for a, b in np.argwhere(~restricted):
            assert not np.any(after[a] & before[b])
        with np.load(folder / f"truth_{i}.npz", allow_pickle=False) as truth:
            assert set(truth.files) == {"points", "amplitude", "rotation_rate"}
            assert (float(truth["amplitude"]), float(truth["rotation_rate"])) == (
                amplitude,
                rotation,
            )
            points = truth["points"]
        wave = RotatingWave(amplitude, rotation)
        if retained is not None:
            with np.load(
                pilot / f"case_{retained}_nested.npz", allow_pickle=False
            ) as old:
                assert np.array_equal(order, old["selection_order"])
                assert np.array_equal(ids, old["selected_ids"])
            assert sha(folder / f"truth_{i}.npz") == sha(
                pilot / f"truth_{retained}_nested.npz"
            )
        else:
            rng = np.random.default_rng([assignment["seed"], 0])
            selected = np.vstack(
                ([0, 0.5, 0, 0], wave.sample(4096, rng), [1, -0.5, 0, 0])
            )
            permutation = np.r_[
                0, rng.permutation(np.arange(1, len(selected) - 1)), len(selected) - 1
            ]
            selected = selected[permutation]
            assert np.array_equal(wave.order(selected), order)
            assert np.array_equal(selected[ids], points)
        events = wave.sample(4096, np.random.default_rng([assignment["seed"], 1]))
        future, past = wave.relations(points[:, None], events[None])
        assert np.array_equal(np.packbits(future, axis=1), after[:, :512])
        assert np.array_equal(np.packbits(past, axis=1), before[:, :512])
        del after, before
        if row.get("failure"):
            rows.append(dict(case=i, failure=row["failure"], passed=False))
            continue
        record = read(folder / f"case_{i}_full.json")
        full = record["full"]
        fit, evaluation = full["fit"], full["evaluation"]
        pure = read(folder / f"case_{i}_full_pure.json")
        boundary = read(folder / f"case_{i}_full_blind_boundary.json")
        assert pure["fit"] == fit and pure["oracle_modules_before_return"] == []
        assert boundary["oracle_modules_loaded"] == []
        assert boundary["opened_archives"] == [str(archive.resolve())]
        assert (
            row["archive_sha256"]
            == record["archive_sha256"]
            == pure["archive_sha256"]
            == boundary["archive_sha256"]
            == fit["order_archive_sha256"]
            == row["artifact_sha256"][archive.name]
        )
        assert fit["selected_ids"] == ids.tolist() and fit["measurement_n"] == n
        assert all(row[k] == v for k, v in evaluation.items())
        assert (
            record["thread_environment"]
            == row["thread_environment"]
            == manifest["thread_environment"]
        )
        for returned in (fit, full["continuum_fit"]):
            assert returned["optimizer_success"] and returned["boundary_excess"] <= 0
            assert returned["coordinate_rank"] == returned["coordinate_parameters"]
            assert returned["projected_curvature_rank"] == 10
            assert returned["nuisance_rank"] == returned["nuisance_parameters"]
            assert returned["measured_intervals"] == int(np.triu(mask, 1).sum())
            assert returned["landmarks"] == len(ids) - 2
            coordinates = np.asarray(returned["coordinates_inferred_from_counts"])
            assert min(abs(coordinates[1, 1]), abs(coordinates[2, 2])) > 1e-12
            covariance = np.asarray(returned["nominal_component_covariance"])
            assert (
                np.isfinite(covariance).all()
                and np.isfinite(returned["coefficients"]).all()
            )
            assert np.linalg.eigvalsh(covariance)[0] > 0
        normal = wave.normal_coordinates(points)
        x = normal[1, 1:].copy()
        x /= np.linalg.norm(x)
        y = normal[2, 1:] - normal[2, 1:] @ x * x
        y /= np.linalg.norm(y)
        q = np.column_stack((x, y, np.cross(x, y)))
        e0 = np.diag([-amplitude, amplitude, 0.0])
        b0 = np.array([[0, -amplitude, 0], [-amplitude, 0, 0], [0, 0, 0]])
        truth_eb = q.T @ e0 @ q, q.T @ b0 @ q
        assert (
            max(
                np.max(abs(a - b))
                for a, b in zip(
                    truth_eb,
                    eb(evaluation["truth_coefficients_xy_positive"]),
                    strict=True,
                )
            )
            < 1e-12
        )
        error = distance(canonical(fit), truth_eb)
        bias = distance(canonical(full["continuum_fit"]), truth_eb)
        estimate_norm = float(
            np.linalg.norm(physical @ np.asarray(fit["coefficients"]))
        )
        variance = float(
            np.linalg.eigvalsh(
                physical @ np.asarray(fit["nominal_component_covariance"]) @ physical.T
            )[-1]
        )
        statistical = math.sqrt(quantile10(0.99375) * variance)
        radius = (statistical + 0.02 * estimate_norm + 1e-10) / 0.98
        old_radius = (
            math.sqrt(quantile10(0.95) * variance) + 0.02 * estimate_norm
        ) / 0.98
        for key, value in (
            ("full_tidal_error_mod_reflection", error),
            ("continuum_bias_mod_reflection", bias),
            ("largest_tidal_variance", variance),
            ("confidence_radius", old_radius),
        ):
            assert abs(evaluation[key] - value) < 1e-12, (i, key)
        assert check_nested_evaluation(full) == (
            bias <= 0.04 * amplitude + 1e-10 and error <= old_radius
        )
        volume = wave.volumes(points)
        coarse = wave.volumes(points, 32)
        refinement = float(np.max(abs(coarse[mask] / volume[mask] - 1)))
        assert abs(refinement - evaluation["quadrature_relative_change"]) < 1e-13
        threshold = 0.5 * amplitude if amplitude else 0.125
        conditions = dict(
            bias=bias <= 0.04 * amplitude + 1e-10,
            quadrature=0 <= refinement < 1e-10,
            contains_truth=error <= radius,
            whole_ball_within_accuracy_target=error + radius <= threshold,
            variance_cap=n * variance <= cap,
        )
        precision = row["precision"]
        assert precision["conditions"] == conditions and precision["passed"] == all(
            conditions.values()
        )
        for key, value in dict(
            nominal_case_coverage=0.99375,
            family_coverage=0.95,
            target_joint_power=0.9,
            radius=radius,
            statistical_radius=statistical,
            error=error,
            truth_norm=2 * amplitude,
            accuracy_threshold=threshold,
            error_plus_radius=error + radius,
            variance_times_N=n * variance,
            variance_cap=cap,
        ).items():
            assert math.isclose(precision[key], value, rel_tol=1e-12, abs_tol=1e-12), (
                i,
                key,
            )
        relative = (error + radius) / (2 * amplitude) if amplitude else None
        if relative is None:
            assert precision["relative_error_plus_radius"] is None
        else:
            assert abs(precision["relative_error_plus_radius"] - relative) < 1e-12
        rows.append(
            dict(
                case=i,
                amplitude=amplitude,
                rotation_rate=rotation,
                measurement_n=n,
                error=error,
                radius=radius,
                bias=bias,
                error_plus_radius=error + radius,
                relative_error_plus_radius=relative,
                variance_times_N=n * variance,
                conditions=conditions,
                passed=all(conditions.values()),
            )
        )
        print(json.dumps({"case_audited": i, "passed": rows[-1]["passed"]}), flush=True)
    verdict = conclusion(rows)
    assert verdict == summary["precision_gate"] == status["precision_gate"]
    for total, field in (
        ("sampled_measurement_events", "sampled_events"),
        ("classified_measurement_events", "classified_events"),
        ("new_selection_events", "selection_events_drawn"),
    ):
        assert summary[total] == sum(
            row["generation"][field] for row in summary["results"]
        )
    returned = [r["case"] for r in summary["results"] if not r.get("failure")]
    replay_case = 5 if 5 in returned else (min(returned) if returned else None)
    replay_receipt = None
    if replay_case is not None:
        print(
            json.dumps({"phase": "post_completion_pure_replay", "case": replay_case}),
            flush=True,
        )
        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--folder",
                str(folder),
                "--replay",
                str(replay_case),
            ],
            check=True,
            capture_output=True,
            text=True,
            creationflags=subprocess.CREATE_NO_WINDOW,
        )
        replay_receipt = json.loads(result.stdout)
    result = dict(
        audit="VERIFIED",
        precision_gate=verdict,
        completed_assignments=8,
        summary_sha256=sha(folder / "summary.json"),
        manifest_sha256=sha(folder / "manifest.json"),
        checker_sha256=sha(Path(__file__)),
        checker_code_version=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        reused_checker_sha256=sha(
            Path(__file__).with_name("check_intrinsic_weyl_coverage.py")
        ),
        preflight_artifact_sha256={
            p.name: sha(p) for p in sorted(folder.glob("preflight_*"))
        },
        selfcheck=checked,
        source_files_checked=len(manifest["source_sha256"]),
        pilot_files_checked=len(manifest["pilot_input_sha256"]),
        planning_evidence_files_checked=len(manifest["planning_evidence_sha256"]),
        case_artifact_hashes_checked=hashes,
        rows=rows,
        post_completion_pure_replay=replay_receipt,
        wall_seconds=time.perf_counter() - wall,
        scope=(
            "All saved tensor metrics and raw-bit transitivity audited; fresh "
            "selections and first 4096 measurement events per returned archive "
            "reproduced on CPU. One separate full pure inverse replay. Existing "
            "truth normal map and volume routines reused; not a full independent "
            "regeneration or proof of Gaussian coverage."
        ),
    )
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "output": str(output),
                "sha256": sha(output),
                "verdict": verdict,
                "wall_seconds": result["wall_seconds"],
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--replay", type=int)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    else:
        assert os.environ.get("CUDA_VISIBLE_DEVICES") == "0"
        if args.replay is not None:
            print(json.dumps(replay(args.folder, args.replay)))
        else:
            if any(x is None for x in (args.folder, args.pilot, args.output)):
                parser.error("folder, pilot and output are required")
            audit(args.folder, args.pilot, args.output)
