"""Read-only verification of the completed frozen coverage campaign."""

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from intrinsic_weyl_joint import confidence_ball, fit_order_archive, tidal_distance
from scipy.stats import beta


def check(folder):
    root = Path(__file__).resolve().parents[2]
    summary = json.loads((folder / "summary.json").read_text())
    manifest = json.loads((folder / "manifest.json").read_text())
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
    amplitudes = np.random.default_rng(2026092250).choice([0.0, 0.25, 0.75], 256)
    rows = summary["cases"]
    assert len(rows) == len(list(folder.glob("case_*_result.json"))) == 256
    checked_hashes = 0
    for i, row in enumerate(rows):
        assert row == json.loads((folder / f"case_{i}_result.json").read_text())
        assert (row["case"], row["seed"], row["amplitude"]) == (
            i,
            2026092300 + i,
            amplitudes[i],
        )
        for name, digest in row["artifact_sha256"].items():
            assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == digest
            checked_hashes += 1
        with np.load(folder / f"case_{i}.npz", allow_pickle=False) as saved:
            assert set(saved.files) == {
                "selection_order",
                "selection_design",
                "selected_ids",
                "after",
                "before",
                "measurement_n",
            }
            assert int(saved["measurement_n"]) == 262144
            assert str(saved["selection_design"]) == "connected"
            assert saved["selection_order"].shape == (4098, 4098)
            assert saved["selection_order"].dtype == np.bool_
            ids, after, before = saved["selected_ids"], saved["after"], saved["before"]
            assert after.shape == before.shape == (len(ids), 32768)
            assert after.dtype == before.dtype == np.uint8
            assert np.all(after[0] == 255) and np.all(before[-1] == 255)
            assert not np.any(before[0]) and not np.any(after[-1])
        if "failure" in row:
            assert not row["valid_and_covered"]
            continue
        full = json.loads((folder / f"case_{i}_full.json").read_text())
        assert full["fit"]["selected_ids"] == ids.tolist()
        inference = confidence_ball(full["fit"])
        assert inference == full["inference"]
        evaluation = full["evaluation"]
        error = tidal_distance(
            inference["coefficients_xy_positive"],
            evaluation["truth_coefficients_xy_positive"],
        )
        assert abs(error - row["full_tidal_error_mod_reflection"]) < 1e-13
        assert all(row[key] == value for key, value in evaluation.items())
        assert row["valid_and_covered"] == (
            error <= inference["radius"] and row["bias_condition_met"]
        )
    successes = sum(row["valid_and_covered"] for row in rows)
    interval = [
        float(beta.ppf(0.025, successes, 257 - successes)),
        float(beta.ppf(0.975, successes + 1, 256 - successes)),
    ]
    assert summary["valid_covered"] == successes
    assert summary["exact_95_two_sided_interval"] == interval
    assert summary["coverage_gate"] == (
        "SUPPORTED" if interval[0] >= 0.95 else "NOT_SUPPORTED"
    )
    # The final case is replayed through the pure estimator module. This checker
    # imports no generator or oracle evaluator; no truth archive is opened.
    replay = fit_order_archive(folder / "case_255.npz")
    original = json.loads((folder / "case_255_full.json").read_text())["fit"]
    difference = max(
        float(np.max(abs(np.array(replay[key]) - np.array(original[key]))))
        for key in (
            "coefficients",
            "nominal_component_covariance",
            "coordinates_inferred_from_counts",
        )
    )
    assert difference < 1e-11
    return {
        "code_version": stamp,
        "cases_checked": len(rows),
        "artifact_hashes_checked": checked_hashes,
        "source_hashes_checked": len(manifest["source_sha256_lf"]),
        "replay_case": 255,
        "replay_max_absolute_difference": difference,
        "valid_covered": successes,
        "exact_95_two_sided_interval": interval,
        "coverage_gate": summary["coverage_gate"],
        "summary_sha256": hashlib.sha256(
            (folder / "summary.json").read_bytes()
        ).hexdigest(),
        "scope": (
            "All archive schemas, receipts and saved confidence calculations; "
            "one independent pure-estimator replay. Does not regenerate all orders."
        ),
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.folder), indent=2))
