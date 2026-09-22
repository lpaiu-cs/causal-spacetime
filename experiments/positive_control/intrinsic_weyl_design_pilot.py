"""Fixed-event paired landmark comparisons with order-only fits."""

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import scipy
from intrinsic_weyl_joint import TIDAL_GRAM, fit_joint, select_landmarks, selection
from intrinsic_weyl_joint_pilot import ROOT
from intrinsic_weyl_rotating_pilot import (
    CASES,
    CHUNK,
    MEASUREMENT_N,
    SELECTION_N,
    evaluate,
)
from rotating_plane_wave import RotatingWave

DESIGNS = ("connected", "wide")


def selfcheck(designs=DESIGNS, *, model_order=4, square_correction=False):
    clouds = []
    for t, number, low, high in (
        (-0.3, 64, 0.12, 0.19),
        (0, 96, 0.10, 0.29),
        (0.3, 64, 0.12, 0.19),
    ):
        k = np.arange(number)
        z = 1 - 2 * (k + 0.5) / number
        phi = k * np.pi * (3 - np.sqrt(5))
        direction = np.column_stack(
            (np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z)
        )
        clouds.append(
            np.column_stack(
                (
                    np.full(number, t),
                    np.linspace(low, high, number)[:, None] * direction,
                )
            )
        )
    points = np.vstack(([-0.5, 0, 0, 0], *clouds, [0.5, 0, 0, 0]))
    delta = points[:, None] - points[None]
    tau2 = delta[..., 0] ** 2 - np.sum(delta[..., 1:] ** 2, axis=-1)
    counts = 1e7 * np.maximum(tau2, 0) ** 2
    chosen = [select_landmarks(counts, design=design) for design in designs]
    old, new = chosen
    assert np.array_equal(old[:3], new[:3])
    if designs[1] == "nested":
        assert np.array_equal(old[:-1], new[: len(old) - 1]) and new[-1] == old[-1]
        assert len(set(new)) == len(new) and 66 < len(new) <= 194
    else:
        for band in (-0.3, 0.3):
            assert np.array_equal(
                old[points[old, 0] == band], new[points[new, 0] == band]
            )
        assert len(set(new)) == len(new) == 66
    radii = [
        float(np.median(np.linalg.norm(points[ids[points[ids, 0] == 0], 1:], axis=1)))
        for ids in chosen
    ]
    if designs[1] == "wide":
        assert radii[1] > radii[0]
    fitted = fit_joint(
        counts[new[:, None], new],
        minimum_count=2,
        starts=2,
        model_order=model_order,
        square_correction=square_correction,
    )
    assert fitted["projected_curvature_rank"] == 10
    assert np.linalg.norm(fitted["coefficients"]) < 1e-7
    if model_order == 5:
        assert fitted["nuisance_rank"] == fitted["nuisance_parameters"]
        assert np.linalg.norm(fitted["gradient_coefficients"]) < 1e-7
    return {
        "same_frame_vectors": True,
        "baseline_is_subset": set(old) <= set(new),
        "selected_landmarks_including_tips": [len(old), len(new)],
        "central_median_radii": radii,
        "exact_flat_projected_rank": 10,
        "exact_flat_coefficient_norm": float(np.linalg.norm(fitted["coefficients"])),
        "exact_flat_nuisance_rank": fitted.get(
            "nuisance_rank", fitted["coordinate_rank"]
        ),
    }


def generate(folder, case, amplitude, rotation_rate, seed, arms, designs=DESIGNS):
    wave = RotatingWave(amplitude, rotation_rate)
    rng = np.random.default_rng([seed, 0])
    points = np.vstack(([0, 0.5, 0, 0], wave.sample(SELECTION_N, rng), [1, -0.5, 0, 0]))
    permutation = np.r_[
        0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
    ]
    points = points[permutation]
    order = wave.order(points)
    selections = {}
    for design in designs:
        try:
            selections[design] = selection(order, design)[0]
        except ValueError as error:
            arms[design]["failure"] = f"selection: {error}"
            np.savez_compressed(
                folder / f"case_{case}_{design}_selection_failed.npz",
                selection_order=order,
                selection_design=design,
            )
    if not selections:
        return {"selection_n": SELECTION_N, "measurement_n_generated": 0}
    union = np.unique(np.concatenate(list(selections.values())))
    after = np.empty((len(union), MEASUREMENT_N // 8), dtype=np.uint8)
    before = np.empty_like(after)
    rng = np.random.default_rng([seed, 1])
    for start in range(0, MEASUREMENT_N, CHUNK):
        events = wave.sample(CHUNK, rng)
        future, past = wave.relations(points[union, None], events[None])
        after[:, start // 8 : (start + CHUNK) // 8] = np.packbits(future, axis=1)
        before[:, start // 8 : (start + CHUNK) // 8] = np.packbits(past, axis=1)
    for design, ids in selections.items():
        indices = np.searchsorted(union, ids)
        np.savez_compressed(
            folder / f"case_{case}_{design}.npz",
            selection_order=order,
            selection_design=design,
            selected_ids=ids,
            after=after[indices],
            before=before[indices],
            measurement_n=MEASUREMENT_N,
        )
        np.savez_compressed(
            folder / f"truth_{case}_{design}.npz",
            points=points[ids],
            amplitude=amplitude,
            rotation_rate=rotation_rate,
        )
        arms[design]["landmarks"] = len(ids) - 2
    return {
        "selection_n": SELECTION_N,
        "measurement_n_generated": MEASUREMENT_N,
        "union_landmarks_including_tips": len(union),
    }


def run(
    folder,
    *,
    design="wide",
    model_order=4,
    square_correction=False,
    data_seed=2026092700,
    protocol="design_pilot_protocol.md",
):
    designs = ("connected", design)
    if design not in ("wide", "nested"):
        raise ValueError("Unknown paired landmark comparison")
    if square_correction and model_order != 5:
        raise ValueError("Curvature-square correction requires model order five")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit the paired protocol and code before data")
    if folder.exists():
        raise FileExistsError("Preserve existing paired pilots")
    checked = selfcheck(
        designs, model_order=model_order, square_correction=square_correction
    )
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    sources = [
        Path(__file__),
        *[
            Path(__file__).with_name(name)
            for name in (
                "intrinsic_weyl_joint.py",
                "intrinsic_weyl_rotating_pilot.py",
                "rotating_plane_wave.py",
                "intrinsic_weyl_joint_audit.py",
                "intrinsic_weyl_joint_pilot.py",
                "intrinsic_weyl.py",
                "window_response.py",
                "p14_plane_wave.py",
            )
        ],
        ROOT / "docs/paper/paper_a/intrinsic_weyl" / protocol,
    ]
    if square_correction:
        sources.append(Path(__file__).with_name("intrinsic_weyl_square.py"))
    manifest = {
        "code_version": stamp,
        "selfcheck": checked,
        "designs": designs,
        "model_order": model_order,
        "curvature_square_correction": square_correction,
        "data_seed": data_seed,
        "cases_fixed_before_data": CASES,
        "selection_n": SELECTION_N,
        "measurement_n": MEASUREMENT_N,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }
    folder.mkdir(parents=True)
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    for case, (a, w) in enumerate(CASES):
        started = time.perf_counter()
        arms = {name: {} for name in designs}
        row = {
            "case": case,
            "amplitude": a,
            "rotation_rate": w,
            "seed": data_seed + case,
            "arms": arms,
        }
        try:
            row["generation"] = generate(folder, case, a, w, row["seed"], arms, designs)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            row["generation_failure"] = f"{type(error).__name__}: {error}"
            for arm in arms.values():
                arm.setdefault("failure", row["generation_failure"])
        for design, arm in arms.items():
            if "failure" in arm:
                continue
            try:
                result = evaluate(
                    folder / f"case_{case}_{design}.npz",
                    folder / f"truth_{case}_{design}.npz",
                    model_order=model_order,
                    square_correction=square_correction,
                )
                (folder / f"case_{case}_{design}_full.json").write_text(
                    json.dumps(result, indent=2) + "\n"
                )
                arm.update(result["evaluation"])
                covariance = np.array(result["fit"]["nominal_component_covariance"])
                for name, slice_ in (
                    ("electric", slice(0, 5)),
                    ("magnetic", slice(5, 10)),
                ):
                    arm[f"{name}_tidal_variance"] = float(
                        np.trace(
                            TIDAL_GRAM[slice_, slice_] @ covariance[slice_, slice_]
                        )
                    )
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                arm["failure"] = f"{type(error).__name__}: {error}"
        if all("failure" not in arm for arm in arms.values()):
            frame_difference = float(
                np.max(
                    abs(
                        np.array(arms[designs[1]]["truth_coefficients_xy_positive"])
                        - np.array(arms["connected"]["truth_coefficients_xy_positive"])
                    )
                )
            )
            if frame_difference >= 1e-10:
                raise ArithmeticError("Paired designs changed the true internal frame")
            row["shared_frame_truth_max_difference"] = frame_difference
            row[f"{designs[1]}_to_connected_ratios"] = {
                k: arms[designs[1]][k] / arms["connected"][k]
                for k in (
                    "magnetic_tidal_variance",
                    "electric_tidal_variance",
                    "largest_tidal_variance",
                    "full_tidal_error_mod_reflection",
                )
            }
        row["seconds"] = time.perf_counter() - started
        row["artifact_sha256"] = {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in folder.glob(f"*_{case}_*")
            if p.suffix in (".npz", ".json")
        }
        rows.append(row)
        (folder / f"case_{case}_result.json").write_text(
            json.dumps(row, indent=2) + "\n"
        )
        print(
            json.dumps({k: v for k, v in row.items() if k != "artifact_sha256"}),
            flush=True,
        )
    summary = manifest | {
        "cases": rows,
        "scope": "Prospective paired design pilot; no accuracy or coverage verdict",
    }
    (folder / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Written {folder / 'summary.json'}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--design", choices=("wide", "nested"), default="wide")
    parser.add_argument("--model-order", type=int, choices=(4, 5), default=4)
    parser.add_argument("--square-correction", action="store_true")
    parser.add_argument("--data-seed", type=int, default=2026092700)
    parser.add_argument("--protocol", default="design_pilot_protocol.md")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/intrinsic_weyl_design_pilot_v1"
    )
    args = parser.parse_args()
    if args.check:
        print(
            json.dumps(
                selfcheck(
                    ("connected", args.design),
                    model_order=args.model_order,
                    square_correction=args.square_correction,
                )
            )
        )
    else:
        run(
            args.output,
            design=args.design,
            model_order=args.model_order,
            square_correction=args.square_correction,
            data_seed=args.data_seed,
            protocol=args.protocol,
        )
