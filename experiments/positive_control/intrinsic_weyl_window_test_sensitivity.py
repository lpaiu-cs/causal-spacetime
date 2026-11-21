"""Post hoc: how much the five-window test depends on its side windows and links.

Not part of the frozen test (``window_test_protocol.md``, commit ``1ca7af2``).
It was added after an independent review and leaves the record's verdicts as
they are. For every replicate of the recorded run it reads the saved pure
fits, the saved correction and the evaluation, and splits the primary reading

    a = <c, sum_s W_s R_s c_s> / <c, c>

into the contributions of the five windows: ``W_s`` are the stencil weights,
``R_s`` the order-only links, ``c_s`` the windows' fits and ``c`` the truth in
the centre window's inferred frame, exactly as the run evaluated it. It then
recomputes the reading for counterfactual corrections: the centre term alone
(side windows deleted), all links replaced by the identity (parity ignored),
and the side links replaced by Haar-random elements of O(3).

Run it from a clean worktree; the output records the commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "outputs" / "intrinsic_weyl_window_test_v1"
RECORDS = ROOT / "docs" / "paper" / "paper_a" / "intrinsic_weyl"
SUMMARY = RECORDS / "window_test_summary.json"
OUTPUT = RECORDS / "window_test_sensitivity.json"
LEVEL = 33554432
REPLICATES = 64
SHIFTS = (-1.0, -0.5, 0.0, 0.5, 1.0)
CENTRE = 0.0
TOLERANCE = 0.05
COVERAGE = 0.95
RANDOM_SEED = 2026112150
RANDOM_DRAWS = 64
AGREEMENT = 1e-9


def git(*args):
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()


def lf_sha(path):
    """Source digest under the repository's LF storage-boundary rule."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def t_interval(values):
    """Mean, sample sd and the run's two-sided Student-t interval."""
    from scipy.stats import t as student

    values = np.asarray(values, dtype=float)
    mean, sd = float(values.mean()), float(values.std(ddof=1))
    quantile = float(student.ppf(0.5 + COVERAGE / 2, len(values) - 1))
    half = quantile * sd / math.sqrt(len(values))
    return {"mean": mean, "sd": sd, "interval": [mean - half, mean + half]}


def accuracy_verdict(interval):
    """G1 of the protocol: strict containment in [1 - t, 1 + t]."""
    low, high = interval
    if low > 1 - TOLERANCE and high < 1 + TOLERANCE:
        return "ACCURATE"
    if high < 1 - TOLERANCE or low > 1 + TOLERANCE:
        return "INACCURATE"
    return "UNRESOLVED"


def judged(values):
    interval = t_interval(values)
    return {**interval, "g1": accuracy_verdict(interval["interval"])}


def haar_o3(rng):
    """A Haar-random element of O(3), parity included."""
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    return q * np.sign(np.diag(r))


def key(shift):
    return f"{shift:+.1f}"


def load(run, replicate):
    prefix = f"rep{replicate:02d}"
    paths = {
        "correction": run / f"{prefix}_n{LEVEL}_combined_pure.json",
        "evaluation": run / f"{prefix}_n{LEVEL}_combined.json",
    }
    for window in range(len(SHIFTS)):
        paths[f"fit{window}"] = run / f"{prefix}_win{window}_n{LEVEL}_pure.json"
    data = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in paths.items()
    }
    digests = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths.values()
    }
    return data, digests


def replicate_row(data, recorded, frames):
    """Per-window split, counterfactuals and checks for one replicate."""
    from intrinsic_weyl_joint import TIDAL_GRAM
    from intrinsic_weyl_window_study import weyl_rotation

    correction = data["correction"]["correction"]
    evaluation = data["evaluation"]["evaluation"]
    truth = np.asarray(evaluation["truth_coefficients_centre_frame"], dtype=float)
    gram = TIDAL_GRAM @ truth
    norm = float(truth @ gram)
    primary = correction["corrections"]["primary"]
    weights = {
        float(shift): np.asarray(matrix, dtype=float)
        for shift, matrix in zip(primary["shifts"], primary["weights"], strict=True)
    }
    rotations = correction["links"]["rotations"]
    links = {s: weyl_rotation(np.asarray(rotations[key(s)])) for s in SHIFTS}
    fits = {
        s: np.asarray(data[f"fit{window}"]["fit"]["coefficients"], dtype=float)
        for window, s in enumerate(SHIFTS)
    }

    def along(vector):
        return float(gram @ vector) / norm

    contribution = {s: along(weights[s] @ links[s] @ fits[s]) for s in SHIFTS}
    reading = {s: along(links[s] @ fits[s]) for s in SHIFTS}
    weight = {s: along(weights[s] @ truth) for s in SHIFTS}
    total = sum(contribution.values())
    single = evaluation["readings"]["single"]["value"]
    checks = {
        "primary_vs_evaluation": abs(
            total - evaluation["readings"]["primary"]["value"]
        ),
        "primary_vs_summary": abs(total - recorded["primary"]["value"]),
        "single_vs_centre_reading": abs(reading[CENTRE] - single),
    }
    if max(checks.values()) > AGREEMENT:
        raise ArithmeticError(f"Decomposition does not reproduce the run: {checks}")
    sides = [s for s in SHIFTS if s != CENTRE]
    random_links = [
        contribution[CENTRE]
        + sum(along(weights[s] @ weyl_rotation(frame[s]) @ fits[s]) for s in sides)
        for frame in frames
    ]
    return {
        "primary": total,
        "single": single,
        "centre_term_only": contribution[CENTRE],
        "identity_links": sum(along(weights[s] @ fits[s]) for s in SHIFTS),
        "contribution": {key(s): contribution[s] for s in SHIFTS},
        "reading_in_linked_frame": {key(s): reading[s] for s in SHIFTS},
        "weight_along_truth": {key(s): weight[s] for s in SHIFTS},
        "largest_check_difference": max(checks.values()),
    }, random_links


def analyse(run, summary_path):
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    recorded = {row["replicate"]: row["readings"] for row in summary["results"]}
    if sorted(recorded) != list(range(REPLICATES)):
        raise ValueError("The summary does not hold the 64 recorded replicates")
    rng = np.random.default_rng(RANDOM_SEED)
    sides = [s for s in SHIFTS if s != CENTRE]
    rows, digests = [], {}
    random_readings = np.empty((RANDOM_DRAWS, REPLICATES))
    for replicate in range(REPLICATES):
        data, files = load(run, replicate)
        digests.update(files)
        frames = [{s: haar_o3(rng) for s in sides} for _ in range(RANDOM_DRAWS)]
        row, random_links = replicate_row(data, recorded[replicate], frames)
        rows.append({"replicate": replicate, **row})
        random_readings[:, replicate] = random_links

    def column(name, shift=None):
        if shift is None:
            return [row[name] for row in rows]
        return [row[name][key(shift)] for row in rows]

    random_draws = [judged(values) for values in random_readings]
    lift = np.subtract(column("primary"), column("single"))
    centre_part = np.subtract(column("centre_term_only"), column("single"))
    side_part = lift - centre_part
    listing = "\n".join(f"{name} {digests[name]}" for name in sorted(digests))
    return {
        "post_hoc": True,
        "description": (
            "Split of the recorded primary reading into the contributions of "
            "the five windows, and the reading under counterfactual corrections. "
            "Added after an independent review; not part of the frozen test, "
            "and the record's verdicts stand as frozen."
        ),
        "frozen_test": {
            "protocol": "docs/paper/paper_a/intrinsic_weyl/window_test_protocol.md",
            "code_version": summary["code_version"],
            "primary_claim": summary["primary_claim"],
        },
        "code_version": git("rev-parse", "--short=7", "HEAD"),
        "source_sha256_lf": lf_sha(Path(__file__)),
        "inputs": {
            "run": run.relative_to(ROOT).as_posix(),
            "summary": summary_path.relative_to(ROOT).as_posix(),
            "files": len(digests),
            "sha256_of_file_digests": hashlib.sha256(listing.encode()).hexdigest(),
        },
        "agreement_with_run": max(row["largest_check_difference"] for row in rows),
        "tolerance": TOLERANCE,
        "coverage_two_sided": COVERAGE,
        "variants": {
            "recorded_primary": judged(column("primary")),
            "single_window": judged(column("single")),
            "centre_term_only": judged(column("centre_term_only")),
            "identity_links": judged(column("identity_links")),
        },
        "random_side_links": {
            "seed": RANDOM_SEED,
            "draws": RANDOM_DRAWS,
            "law": "Haar on O(3), independent per replicate and side window",
            "g1_accurate": sum(d["g1"] == "ACCURATE" for d in random_draws),
            "mean_range": [
                min(d["mean"] for d in random_draws),
                max(d["mean"] for d in random_draws),
            ],
            "lowest_lower_bound": min(d["interval"][0] for d in random_draws),
            "highest_upper_bound": max(d["interval"][1] for d in random_draws),
        },
        "windows": {
            key(s): {
                "contribution": t_interval(column("contribution", s)),
                "reading_in_linked_frame": t_interval(
                    column("reading_in_linked_frame", s)
                ),
                "weight_along_truth": t_interval(column("weight_along_truth", s)),
            }
            for s in SHIFTS
        },
        "lift": {
            "total": t_interval(lift),
            "centre_term_minus_single": t_interval(centre_part),
            "side_windows": t_interval(side_part),
        },
        "replicates": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", type=Path, default=RUN)
    parser.add_argument("--summary", type=Path, default=SUMMARY)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--allow-dirty", action="store_true", help="debug runs only; recorded false"
    )
    args = parser.parse_args()
    clean = git("status", "--porcelain", "--untracked-files=all") == ""
    if not clean and not args.allow_dirty:
        raise SystemExit("Refusing to run from a dirty worktree")
    result = analyse(args.run.resolve(), args.summary.resolve())
    result = {"worktree_clean_at_start": clean, **result}
    args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    for name, value in result["variants"].items():
        low, high = value["interval"]
        print(
            f"{name:18s} {value['mean']:.4f} [{low:.4f}, {high:.4f}] "
            f"sd {value['sd']:.4f} {value['g1']}"
        )
    block = result["random_side_links"]
    print(
        f"random side links: {block['g1_accurate']}/{block['draws']} ACCURATE, "
        f"means {block['mean_range'][0]:.4f}-{block['mean_range'][1]:.4f}"
    )
    for shift, value in result["windows"].items():
        print(
            f"window {shift}: contribution {value['contribution']['mean']:+.4f}, "
            f"reading {value['reading_in_linked_frame']['mean']:+.4f}, "
            f"weight {value['weight_along_truth']['mean']:+.4f}"
        )
    print(f"agreement with the run: {result['agreement_with_run']:.1e}")


if __name__ == "__main__":
    main()
