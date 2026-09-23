"""Finish the fixed high-count ladder without replacing its saved prefixes."""

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from intrinsic_weyl_high_count import (
    CASES,
    CHUNK,
    LEVELS,
    ROOT,
    SEED,
    THREAD_VARIABLES,
    differences,
    peak_working_set,
    sha,
    write_json,
    write_prefix,
)

EVIDENCE_SHA = "5c68aef9d1f47c0aeb8e1fee88991932e768cf4047af333921e29a3b98dfee91"
PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/high_count_continuation.md"


def check_prefix(path, selected, after, before, n):
    with np.load(path, allow_pickle=False) as saved:
        assert set(saved.files) == set(selected) | {"after", "before", "measurement_n"}
        assert int(saved["measurement_n"]) == n
        assert all(np.array_equal(saved[k], v) for k, v in selected.items())
        assert np.array_equal(saved["after"], after[:, : n // 8])
        assert np.array_equal(saved["before"], before[:, : n // 8])


def selfcheck():
    import tempfile

    selected = {
        "selection_order": np.triu(np.ones((3, 3), bool), 1),
        "selection_design": np.array("nested"),
        "selected_ids": np.arange(3),
    }
    after = np.arange(12, dtype=np.uint8).reshape(3, 4)
    before = 255 - after
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "prefix.npz"
        write_prefix(path, selected, after, before, 24)
        check_prefix(path, selected, after, before, 24)
        changed = after.copy()
        changed[1, 1] ^= 1
        try:
            check_prefix(path, selected, changed, before, 24)
        except AssertionError:
            pass
        else:
            raise AssertionError("A changed saved prefix was accepted")
    return {"unchanged_prefix_accepted": True, "changed_bit_rejected": True}


def run(source, pilot, destination, evidence_path):
    assert os.name == "nt" and os.environ.get("CUDA_VISIBLE_DEVICES") == "-1"
    assert sha(evidence_path) == EVIDENCE_SHA
    evidence = json.loads(evidence_path.read_text())
    preserved = evidence["original_file_sha256"]
    assert len(preserved) == 132
    assert {p.name for p in source.iterdir()} == set(preserved)
    assert evidence["classification"] == "INCOMPLETE_AFTER_INTERRUPTION"
    assert evidence["completed_assignments"] == 21 and evidence["remaining_case"] == 7
    assert not (source / "summary.json").exists()
    manifest = json.loads((source / "manifest.json").read_text())
    assert manifest["cases"] == [list(c) for c in CASES]
    assert manifest["measurement_levels"] == list(LEVELS)
    assert manifest["measurement_seed_base"] == SEED and manifest["chunk"] == CHUNK
    assert manifest["runtime"] == {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    assert {k: os.environ.get(k) for k in THREAD_VARIABLES} == manifest[
        "thread_environment"
    ]
    assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", manifest["code_version"], stamp],
        cwd=ROOT,
        check=True,
    )
    coordinator = {
        str(p.relative_to(ROOT)).replace("\\", "/"): sha(p)
        for p in (Path(__file__), ROOT / PROTOCOL)
    }

    def verify(copied=False):
        for name, digest in preserved.items():
            assert sha(source / name) == digest, name
            if copied:
                assert sha(destination / name) == digest, name
        for name, digest in manifest["source_sha256"].items():
            assert sha(ROOT / name) == digest, name
        for name, digest in manifest["pilot_input_sha256"].items():
            assert sha(pilot / name) == digest, name
        for name, digest in coordinator.items():
            assert sha(ROOT / name) == digest, name

    verify()
    destination.mkdir(parents=True, exist_ok=False)
    for name in preserved:
        shutil.copyfile(source / name, destination / name)
    verify(copied=True)
    started = time.time()

    def status(**fields):
        (destination / "coordinator_status.json").write_text(
            json.dumps(
                {"pid": os.getpid(), "started_unix": started, **fields}, indent=2
            )
            + "\n"
        )

    def invoke(archive, truth, output):
        with output.with_suffix(".log").open("xb") as log:
            child = subprocess.Popen(
                [
                    sys.executable,
                    str(
                        ROOT
                        / "experiments/positive_control/intrinsic_weyl_high_count.py"
                    ),
                    "--worker",
                    str(archive),
                    "--truth",
                    str(truth),
                    "--output",
                    str(output),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
            status(
                state="RUNNING", phase="fit", output=output.name, child_pid=child.pid
            )
            code = child.wait()
        if code:
            return {
                "failure": f"Worker exit {code}; see {output.with_suffix('.log').name}"
            }
        return json.loads(output.read_text())

    receipt = {
        "coordinator_code_version": stamp,
        "coordinator_source_sha256": coordinator,
        "original_directory": str(source.resolve()),
        "evidence_sha256": EVIDENCE_SHA,
        "preserved_sha256": preserved,
        "completed_cases_reused": list(range(7)),
        "remaining_case": 7,
        "seed": SEED + 7,
        "thread_environment": manifest["thread_environment"],
        "gpu_used": False,
        "selfcheck": selfcheck(),
    }
    write_json(destination / "continuation.json", receipt)
    status(state="RUNNING", phase="preflight")
    try:
        replay = invoke(
            destination / "case_6_n262144.npz",
            pilot / "truth_6_nested.npz",
            destination / "continuation_replay.json",
        )
        original = json.loads((source / "case_6_n262144.json").read_text())
        mismatch = differences(replay.get("full"), original["full"])
        write_json(
            destination / "continuation_replay_comparison.json",
            {
                "case": 6,
                "measurement_n": 262144,
                "differences": mismatch,
                "within_frozen_tolerance": not mismatch,
                "all_full_fields_exact": replay.get("full") == original["full"],
                "fresh_events": 0,
            },
        )
        if mismatch:
            raise RuntimeError("Saved-data replay failed before regeneration")

        from rotating_plane_wave import RotatingWave

        with np.load(pilot / "case_7_nested.npz", allow_pickle=False) as saved:
            selected = {
                k: saved[k]
                for k in ("selection_order", "selection_design", "selected_ids")
            }
        amplitude, rotation = CASES[7]
        with np.load(pilot / "truth_7_nested.npz", allow_pickle=False) as truth:
            assert (float(truth["amplitude"]), float(truth["rotation_rate"])) == (
                amplitude,
                rotation,
            )
            points = truth["points"]
        after = np.empty((len(points), LEVELS[-1] // 8), dtype=np.uint8)
        before = np.empty_like(after)
        wave, rng = (
            RotatingWave(amplitude, rotation),
            np.random.default_rng([SEED + 7, 1]),
        )
        generation = {
            "case": 7,
            "sampled_events": 0,
            "classified_events": 0,
            "prefix_costs": [],
            "saved_prefixes_reconstructed_bit_exact": [],
        }
        wall, cpu = time.perf_counter(), time.process_time()
        status(state="RUNNING", phase="generate", case=7, events=0)
        try:
            for start in range(0, LEVELS[-1], CHUNK):
                events = wave.sample(CHUNK, rng)
                generation["sampled_events"] += CHUNK
                future, past = wave.relations(points[:, None], events[None])
                end = start + CHUNK
                after[:, start // 8 : end // 8] = np.packbits(future, axis=1)
                before[:, start // 8 : end // 8] = np.packbits(past, axis=1)
                generation["classified_events"] = end
                if end in LEVELS:
                    path = destination / f"case_7_n{end}.npz"
                    if path.exists():
                        check_prefix(path, selected, after, before, end)
                        generation["saved_prefixes_reconstructed_bit_exact"].append(end)
                    else:
                        write_prefix(path, selected, after, before, end)
                    generation["prefix_costs"].append(
                        {
                            "measurement_n": end,
                            "wall_seconds": time.perf_counter() - wall,
                            "process_cpu_seconds": time.process_time() - cpu,
                        }
                    )
                if end % (64 * CHUNK) == 0:
                    status(state="RUNNING", phase="generate", case=7, events=end)
        except (ValueError, ArithmeticError, MemoryError) as error:
            generation["failure"] = f"{type(error).__name__}: {error}"
        generation.update(
            wall_seconds=time.perf_counter() - wall,
            process_cpu_seconds=time.process_time() - cpu,
            parent_lifetime_peak_working_set_bytes=peak_working_set(),
            directional_predicates_classified=2
            * len(points)
            * generation["classified_events"],
            directional_predicates_attempted=2
            * len(points)
            * generation["sampled_events"],
        )
        write_json(destination / "case_7_generation.json", generation)
        del after, before
        rows = [
            json.loads((destination / f"case_{i}_n{n}_result.json").read_text())
            for i in range(7)
            for n in LEVELS
        ]
        for n in LEVELS:
            archive = destination / f"case_7_n{n}.npz"
            row = {
                "case": 7,
                "amplitude": amplitude,
                "rotation_rate": rotation,
                "seed": SEED + 7,
                "measurement_n": n,
            }
            if not archive.exists():
                row["failure"] = generation.get(
                    "failure", "Measurement prefix not completed"
                )
            else:
                output = archive.with_suffix(".json")
                result = invoke(archive, pilot / "truth_7_nested.npz", output)
                row.update({k: v for k, v in result.items() if k != "full"})
                if "full" in result:
                    row.update(result["full"]["evaluation"])
                    row["variance_price_V"] = n * row["largest_tidal_variance"]
                row.update(
                    archive_sha256=sha(archive), archive_bytes=archive.stat().st_size
                )
                if output.exists():
                    row["full_record_sha256"] = sha(output)
            write_json(destination / f"case_7_n{n}_result.json", row)
            rows.append(row)
            print(json.dumps(row), flush=True)
        verify(copied=True)
        generated = [
            json.loads((destination / f"case_{i}_generation.json").read_text())
            for i in range(8)
        ]
        largest_archived = max(
            [LEVELS[1]] + [p["measurement_n"] for p in generation["prefix_costs"]]
        )
        write_json(
            destination / "summary.json",
            {
                **manifest,
                "results": rows,
                "generation": generated,
                "completed_assignments": len(rows),
                "continuation": receipt,
                "unique_archived_measurement_events": 7 * LEVELS[-1] + largest_archived,
                "sample_operations_lower_bound": sum(
                    g["sampled_events"] for g in generated
                )
                + evidence["last_saved_status"]["events"],
                "interrupted_case_7_cost": (
                    "At least 3670016 draws/classifications were observed before "
                    "interruption; its complete CPU/wall cost is unavailable. "
                    "Replayed work is not free."
                ),
                "scope": (
                    "All fixed assignments processed across an interruption; "
                    "descriptive diagnostic, not powered accuracy or coverage"
                ),
            },
        )
        status(state="COMPLETED", assignments=len(rows))
    except BaseException as error:
        status(state="STOPPED", error=f"{type(error).__name__}: {error}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    else:
        if any(
            v is None
            for v in (args.source, args.pilot, args.destination, args.evidence)
        ):
            parser.error("source, pilot, destination and evidence are required")
        run(args.source, args.pilot, args.destination, args.evidence)
