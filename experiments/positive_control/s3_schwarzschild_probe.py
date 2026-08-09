"""S3: Schwarzschild paired exploration probe (S-track, exploratory).

Measures the paired same-point-set shift Delta = f_M - f_0 between the
Schwarzschild causal predicate (the S1 solver on its frozen domain)
and the flat predicate on IDENTICAL events. The design rests on the
measure identity: in Schwarzschild coordinates sqrt(-g) = r^2
sin(theta), independent of M and equal to the flat spherical
Minkowski volume element, so one sprinkle serves both geometries --
the exact analogue of det g = -1 in the plane-wave chain (S2 design
memo). What the identity controls is the sampling MEASURE; individual
diamond volumes differ under the two predicates, and that difference
is part of the signal, not a confound.

Exploration only -- no preregistered gate. The output (sign, size and
spread of Delta) sizes a later operationally-anchored confirmation in
the P3-C pattern, which needs no diamond-volume oracle.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import s1_schwarzschild_cost as s1
from p14_probe_p2 import student_t_crit
from seed_windows import (
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_point_seeds_fresh,
)

# ---------------------------------------------------------------------
# Exploration constants (NOT frozen -- this is a probe, not a stage)
# ---------------------------------------------------------------------

N_EVENTS = 300               # events per reading (S1 bench scale)
N_READINGS = 300             # paired readings
SEED = 40_000_201            # fresh vs every ledger, see SPENT_SCALARS
TOL = s1.DEFAULT_TOL         # 1e-8, the priced S1 operating point
TOL_ESCALATED = 1e-10        # last rung of s1.TOL_LADDER

#: Every scalar seed the program has spent, across all chains: the
#: P14 probe/campaign families (777-781 including the (781, k)
#: SeedSequence roots, 2026xxxx probe seeds) and the 40M block
#: (C1/C2 execution, S1 bench).
SPENT_SCALARS = frozenset({
    777, 778, 779, 780, 781,
    20_260_808, 20_260_811, 20_260_812, 20_260_813, 20_260_814,
    20_260_821, 20_260_822, 20_260_823, 20_260_824,
    20_260_851, 20_260_852,
    40_000_061, 40_000_071, 40_000_072, 40_000_101,
})

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_s3_probe_results.json")


# ---------------------------------------------------------------------
# Flat predicate: Minkowski in spherical coordinates, closed form
# ---------------------------------------------------------------------


def flat_chord(p: np.ndarray, q: np.ndarray) -> float:
    """Euclidean chord between the spatial points of two events
    (t, r, theta, phi) -- the exact flat light-cone radius."""

    cosang = (math.sin(p[2]) * math.sin(q[2]) * math.cos(p[3] - q[3])
              + math.cos(p[2]) * math.cos(q[2]))
    cosang = max(-1.0, min(1.0, cosang))
    return math.sqrt(max(0.0,
                         p[1] * p[1] + q[1] * q[1]
                         - 2.0 * p[1] * q[1] * cosang))


def flat_relation(p: np.ndarray, q: np.ndarray) -> bool:
    """`p prec q` under Minkowski, events time-ordered by caller."""

    return (q[0] - p[0]) > flat_chord(p, q)


# ---------------------------------------------------------------------
# One paired reading: sprinkle once, census under both geometries
# ---------------------------------------------------------------------


def reading(rng: np.random.Generator,
            n_events: int = N_EVENTS) -> tuple[float, float, int, int]:
    """(f_M, f_0, ambiguous, escalated) on one shared sprinkle.

    The Schwarzschild side inherits the undecided-is-never-a-silent-
    False contract: an undecided pair at TOL escalates to
    TOL_ESCALATED; a pair still undecided there is COUNTED as
    ambiguous (and surfaced in the artifact) before being left
    unrelated, bounding the bias by ambiguous/pairs."""

    pts = s1.sample_events(n_events, rng)
    pts = pts[np.argsort(pts[:, 0], kind="stable")]
    pairs = n_events * (n_events - 1) // 2
    rel_m = rel_0 = ambiguous = escalated = 0
    for i in range(n_events - 1):
        p = pts[i]
        for j in range(i + 1, n_events):
            q = pts[j]
            if flat_relation(p, q):
                rel_0 += 1
            rel = s1.causal_relation(p, q, tol=TOL)
            if rel is None:
                escalated += 1
                rel = s1.causal_relation(p, q, tol=TOL_ESCALATED)
                if rel is None:
                    ambiguous += 1
                    rel = False
            if rel:
                rel_m += 1
    return rel_m / pairs, rel_0 / pairs, ambiguous, escalated


# ---------------------------------------------------------------------
# The probe
# ---------------------------------------------------------------------


def _git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"],
                           capture_output=True, text=True).stdout.strip()
    return {"rev": rev, "dirty": bool(dirty)}


def main() -> None:
    smoke = "--smoke" in sys.argv
    n_events = 40 if smoke else N_EVENTS
    n_readings = 3 if smoke else N_READINGS

    assert_point_seeds_fresh(
        {"s3_exploration": SEED}, SPENT_SCALARS,
        P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,), "S3")

    rng = np.random.default_rng(SEED)
    fm = np.empty(n_readings)
    f0 = np.empty(n_readings)
    ambiguous = escalated = 0
    start = time.perf_counter()
    for k in range(n_readings):
        fm[k], f0[k], amb, esc = reading(rng, n_events)
        ambiguous += amb
        escalated += esc
        if (k + 1) % 10 == 0 or k == 0 or (k + 1) == n_readings:
            dt = time.perf_counter() - start
            eta = dt / (k + 1) * (n_readings - k - 1)
            print(f"  s3 {k + 1}/{n_readings}  "
                  f"delta so far {np.mean(fm[:k + 1] - f0[:k + 1]):+.6f}  "
                  f"elapsed {dt / 60:.1f} min, eta {eta / 60:.1f} min",
                  flush=True)
    runtime = time.perf_counter() - start

    delta = fm - f0
    n = n_readings
    mean = float(np.mean(delta))
    sd = float(np.std(delta, ddof=1))
    sem = sd / math.sqrt(n)
    tcrit = student_t_crit(n - 1)
    ci = [mean - tcrit * sem, mean + tcrit * sem]

    result = {
        "probe": "S3 Schwarzschild paired exploration",
        "design": "one sprinkle, dual census (measure identity, S2 memo)",
        "domain": {"m": s1.M, "r_shell": [s1.R_MIN, s1.R_MAX],
                   "cap_half_angle": s1.CAP_HALF_ANGLE,
                   "t_extent": s1.T_EXTENT},
        "params": {"n_events": n_events, "n_readings": n_readings,
                   "seed": SEED, "tol": TOL,
                   "tol_escalated": TOL_ESCALATED, "smoke": smoke},
        "delta": {"mean": mean, "sd": sd, "sem": sem,
                  "ci95_student_t": ci,
                  "per_reading": [float(d) for d in delta]},
        "f_flat": {"mean": float(np.mean(f0)), "sd": float(np.std(f0, ddof=1))},
        "f_schwarzschild": {"mean": float(np.mean(fm)),
                            "sd": float(np.std(fm, ddof=1))},
        "ambiguity": {"ambiguous": ambiguous, "escalated": escalated,
                      "pairs_per_reading": n_events * (n_events - 1) // 2},
        "runtime_seconds": runtime,
        "code": _git_state(),
    }

    print(f"\nS3: delta mean {mean:+.6f}  sd {sd:.6f}  "
          f"CI95 [{ci[0]:+.6f}, {ci[1]:+.6f}]  "
          f"f0 {np.mean(f0):.6f} -> fM {np.mean(fm):.6f}  "
          f"ambiguous {ambiguous}, escalated {escalated}  "
          f"({runtime / 3600:.2f} h)")

    if smoke:
        print("smoke run -- artifact NOT written")
        return
    _ARTIFACT.write_text(json.dumps(result, indent=2) + "\n",
                         encoding="utf-8")
    print(f"wrote {_ARTIFACT}")


if __name__ == "__main__":
    main()
