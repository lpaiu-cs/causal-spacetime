"""Read-only verification of the completed frozen coverage campaign."""

import hashlib
import json
import subprocess
from pathlib import Path

import intrinsic_weyl_joint as joint
import numpy as np
from intrinsic_weyl_joint import confidence_ball, fit_order_archive, tidal_distance
from scipy.stats import beta


def coverage_interval(covered, expected):
    """Independent exact interval; incomplete and failed cases cannot vanish."""
    if len(covered) != expected:
        raise ValueError(f"Expected all {expected} assigned cases")
    count = sum(covered)
    return [
        float(beta.ppf(0.025, count, expected + 1 - count)) if count else 0.0,
        float(beta.ppf(0.975, count + 1, expected - count))
        if count < expected
        else 1.0,
    ]


def check_nested_evaluation(full):
    """Recompute both errors and their flags from the two saved full fits."""
    finite = confidence_ball(full["fit"])
    continuum = confidence_ball(full["continuum_fit"])
    assert finite == full["inference"]
    for fit in (full["fit"], full["continuum_fit"]):
        assert fit["model_order"] == 5 and fit["curvature_square_correction"]
    evaluation = full["evaluation"]
    truth = evaluation["truth_coefficients_xy_positive"]
    norm = float(np.linalg.norm(joint.TIDAL_ROOT @ np.asarray(truth)))
    error = tidal_distance(finite["coefficients_xy_positive"], truth)
    bias = tidal_distance(continuum["coefficients_xy_positive"], truth)
    assert abs(norm - evaluation["truth_norm"]) < 1e-13
    assert abs(error - evaluation["full_tidal_error_mod_reflection"]) < 1e-13
    assert abs(bias - evaluation["continuum_bias_mod_reflection"]) < 1e-13
    for field, value in (
        ("relative_full_tidal_error", error),
        ("relative_continuum_bias", bias),
    ):
        if norm:
            assert abs(evaluation[field] - value / norm) < 1e-13
        else:
            assert evaluation[field] is None
    # Preserve the frozen evaluator's numerical slack, including curved cases.
    bias_ok = bias <= 0.02 * evaluation["truth_norm"] + 1e-10
    covered = error <= finite["radius"]
    assert evaluation["bias_condition_met"] == bias_ok
    assert evaluation["nominal_ball_contains_truth"] == covered
    assert evaluation["confidence_radius"] == finite["radius"]
    assert evaluation["largest_tidal_variance"] == finite["largest_tidal_variance"]
    assert 0 <= evaluation["quadrature_relative_change"] < 1e-10
    return bias_ok and covered


def selfcheck():
    assert coverage_interval([True] * 497 + [False] * 15, 512)[0] >= 0.95
    assert coverage_interval([True] * 496 + [False] * 16, 512)[0] < 0.95
    assert coverage_interval([False] * 512, 512)[0] == 0
    assert coverage_interval([True] * 512, 512)[1] == 1
    try:
        coverage_interval([True] * 511, 512)
    except ValueError:
        pass
    else:
        raise AssertionError("Incomplete campaigns must be rejected")
    fit = {
        "optimizer_success": True,
        "boundary_excess": 0,
        "coordinate_rank": 3,
        "coordinate_parameters": 3,
        "projected_curvature_rank": 10,
        "nuisance_rank": 27,
        "nuisance_parameters": 27,
        "coordinates_inferred_from_counts": [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]],
        "coefficients": [0.0] * 10,
        "nominal_component_covariance": np.eye(10).tolist(),
        "model_order": 5,
        "curvature_square_correction": True,
    }
    inference = confidence_ball(fit)
    evaluation = {
        "truth_coefficients_xy_positive": [0.0] * 10,
        "truth_norm": 0.0,
        "full_tidal_error_mod_reflection": 0.0,
        "relative_full_tidal_error": None,
        "continuum_bias_mod_reflection": 0.0,
        "relative_continuum_bias": None,
        "bias_condition_met": True,
        "nominal_ball_contains_truth": True,
        "confidence_radius": inference["radius"],
        "largest_tidal_variance": inference["largest_tidal_variance"],
        "quadrature_relative_change": 0.0,
    }
    full = dict(fit=fit, continuum_fit=fit, inference=inference, evaluation=evaluation)
    assert check_nested_evaluation(full)
    for field in ("bias_condition_met", "nominal_ball_contains_truth"):
        corrupted = full | {"evaluation": evaluation | {field: False}}
        try:
            check_nested_evaluation(corrupted)
        except AssertionError:
            pass
        else:
            raise AssertionError(f"A corrupted {field} flag was accepted")
    return {"gate_boundaries": "passed", "incomplete_and_corrupt_records": "rejected"}


def check(folder):
    # Validate the actual imported inverse, including when this read-only
    # checker is launched from an isolated worktree against frozen sources.
    root = Path(joint.__file__).resolve().parents[2]
    summary = json.loads((folder / "summary.json").read_text())
    manifest = json.loads((folder / "manifest.json").read_text())
    nested = manifest.get("selection_design") == "nested"
    n = 512 if nested else 256
    suffix = "_nested" if nested else ""
    data_seed = 2026101000 if nested else 2026092300
    if nested:
        assert manifest["population_seed"] == 2026100950
        assert manifest["data_seed"] == data_seed
        assert manifest["selection_n"] == 4096
        assert manifest["measurement_n"] == 262144
        assert manifest["model_order"] == 5
        assert manifest["curvature_square_correction"] is True
        assert manifest["workers"] == 4
        assert manifest["thread_settings_changed"] is False
        assert manifest["worker_thread_environment"] == {
            "OPENBLAS_NUM_THREADS": None,
            "OMP_NUM_THREADS": None,
        }
        rng = np.random.default_rng(2026100950)
        amplitudes = rng.choice([0.0, 0.25, 0.75], n)
        rotations = rng.choice([-0.5, 0.0, 0.5], n)
        assert (
            manifest["population"] == np.column_stack((amplitudes, rotations)).tolist()
        )
    else:
        if manifest.get("model_order", 4) != 4:
            raise ValueError("This checker supports the original and nested campaigns")
        amplitudes = np.random.default_rng(2026092250).choice([0.0, 0.25, 0.75], n)
    stamp = manifest["code_version"]
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=root, check=True
    )
    for name, digest in manifest["source_sha256_lf"].items():
        historical = subprocess.check_output(
            ["git", "show", f"{stamp}:{name}"], cwd=root
        )
        assert hashlib.sha256(historical.replace(b"\r\n", b"\n")).hexdigest() == digest
        assert (
            hashlib.sha256(
                (root / name).read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            == digest
        )
    for key, value in manifest.items():
        assert summary[key] == value
    if "continuation" in summary:
        receipt = json.loads((folder / "continuation.json").read_text())
        assert summary["continuation"] == receipt
        assert receipt["source_code_version"] == stamp
        assert receipt["original_status"] == "INCOMPLETE_AFTER_INTERRUPTION"
        assert receipt["completed_cases_reused"] == list(range(226))
        assert receipt["saved_inputs_reused"] == list(range(226, 230))
        assert receipt["new_cases"] == list(range(230, 512))
        original_directory = Path(receipt["original_directory"])
        assert not (original_directory / "summary.json").exists()
        preserved = receipt["preserved_sha256"]
        assert len(preserved) == 913
        assert set(preserved) == {p.name for p in original_directory.iterdir()}
        for name, digest in preserved.items():
            assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == digest
            assert (
                hashlib.sha256((original_directory / name).read_bytes()).hexdigest()
                == digest
            )
        coordinator_root = Path(__file__).resolve().parents[2]
        coordinator_stamp = receipt["coordinator_code_version"]
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", coordinator_stamp, "HEAD"],
            cwd=coordinator_root,
            check=True,
        )
        historical = subprocess.check_output(
            ["git", "show", f"{coordinator_stamp}:{receipt['coordinator_source']}"],
            cwd=coordinator_root,
        )
        assert (
            hashlib.sha256(historical.replace(b"\r\n", b"\n")).hexdigest()
            == (receipt["coordinator_sha256_lf"])
        )
    rows = summary["cases"]
    assert len(rows) == len(list(folder.glob("case_*_result.json"))) == n
    if nested:
        assert summary["n"] == n
    checked_hashes = 0
    covered = []
    replay_candidates = []
    for i, row in enumerate(rows):
        assert row == json.loads((folder / f"case_{i}_result.json").read_text())
        assert (row["case"], row["seed"], row["amplitude"]) == (
            i,
            data_seed + i,
            amplitudes[i],
        )
        if nested:
            assert row["rotation_rate"] == rotations[i]
        for name, digest in row["artifact_sha256"].items():
            assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == digest
            checked_hashes += 1
        archive = folder / f"case_{i}{suffix}.npz"
        if not archive.exists():
            assert nested and row.get("failure")
            covered.append(False)
            continue
        assert archive.name in row["artifact_sha256"]
        with np.load(archive, allow_pickle=False) as saved:
            assert set(saved.files) == {
                "selection_order",
                "selection_design",
                "selected_ids",
                "after",
                "before",
                "measurement_n",
            }
            assert int(saved["measurement_n"]) == 262144
            assert str(saved["selection_design"]) == (
                "nested" if nested else "connected"
            )
            assert saved["selection_order"].shape == (4098, 4098)
            assert saved["selection_order"].dtype == np.bool_
            ids, after, before = saved["selected_ids"], saved["after"], saved["before"]
            assert after.shape == before.shape == (len(ids), 32768)
            assert after.dtype == before.dtype == np.uint8
            assert np.all(after[0] == 255) and np.all(before[-1] == 255)
            assert not np.any(before[0]) and not np.any(after[-1])
        if "failure" in row:
            assert row["failure"]
            assert not row.get("valid_and_covered", False)
            covered.append(False)
            continue
        full = json.loads((folder / f"case_{i}_full.json").read_text())
        assert f"case_{i}_full.json" in row["artifact_sha256"]
        assert f"truth_{i}{suffix}.npz" in row["artifact_sha256"]
        assert full["fit"]["selected_ids"] == ids.tolist()
        assert full["fit"]["measurement_n"] == 262144
        assert (
            full["fit"]["order_archive_sha256"] == row["artifact_sha256"][archive.name]
        )
        inference = confidence_ball(full["fit"])
        assert inference == full["inference"]
        evaluation = full["evaluation"]
        if nested:
            assert evaluation["truth_norm"] == 2 * row["amplitude"]
        error = tidal_distance(
            inference["coefficients_xy_positive"],
            evaluation["truth_coefficients_xy_positive"],
        )
        assert abs(error - row["full_tidal_error_mod_reflection"]) < 1e-13
        assert all(row[key] == value for key, value in evaluation.items())
        valid = error <= inference["radius"] and row["bias_condition_met"]
        if nested:
            assert check_nested_evaluation(full) == valid
        else:
            assert row["valid_and_covered"] == valid
        covered.append(valid)
        replay_candidates.append((not valid, i))
    successes = sum(covered)
    interval = coverage_interval(covered, n)
    assert summary["valid_covered"] == successes
    assert summary["exact_95_two_sided_interval"] == interval
    assert summary["coverage_gate"] == (
        "SUPPORTED" if interval[0] >= 0.95 else "NOT_SUPPORTED"
    )
    # Prefer an unfavorable returned fit, otherwise the last returned fit.
    # No generator or oracle evaluator is imported; truth archives stay closed.
    if not replay_candidates:
        raise ValueError("No returned fit is available for the required pure replay")
    replay_case = max(replay_candidates)[1] if nested else 255
    replay = fit_order_archive(
        folder / f"case_{replay_case}{suffix}.npz",
        model_order=5 if nested else 4,
        square_correction=nested,
    )
    original = json.loads((folder / f"case_{replay_case}_full.json").read_text())["fit"]
    difference = max(
        float(np.max(abs(np.array(replay[key]) - np.array(original[key]))))
        for key in (
            "coefficients",
            "nominal_component_covariance",
            "coordinates_inferred_from_counts",
        )
    )
    assert difference < 1e-11
    if nested:
        assert replay == original
    return {
        "code_version": stamp,
        "cases_checked": len(rows),
        "artifact_hashes_checked": checked_hashes,
        "source_hashes_checked": len(manifest["source_sha256_lf"]),
        "replay_case": replay_case,
        "replay_max_absolute_difference": difference,
        "valid_covered": successes,
        "exact_95_two_sided_interval": interval,
        "coverage_gate": summary["coverage_gate"],
        "summary_sha256": hashlib.sha256(
            (folder / "summary.json").read_bytes()
        ).hexdigest(),
        "checker_sha256_lf": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest(),
        "scope": (
            "Available estimator-input schemas, receipts and saved confidence "
            "calculations; nested bias flags recomputed from saved continuum fits; "
            "one pure-estimator replay. Does not regenerate orders or oracle truth."
        ),
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path, nargs="?")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    elif args.folder is None:
        parser.error("A completed campaign folder is required")
    else:
        print(json.dumps(check(args.folder), indent=2))
