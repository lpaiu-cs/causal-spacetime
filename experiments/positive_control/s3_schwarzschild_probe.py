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

Scope discipline: Delta is a statement about the FROZEN Schwarzschild
coordinates and S1 domain (exterior shell, polar cap, coordinate-time
box). The Shapiro-delay sign argument behind the expected Delta < 0
is coordinate-tied, so every downstream sentence keeps the frozen
coordinates and domain explicit.

Sampling law: the event count per reading is N ~ Poisson(E_N), so the
point process is the Poisson process of the common volume measure
(S2/PR review R1 -- a fixed N would be a binomial process conditioned
on N, and its variance would not transfer to Poisson-based sizing).

Ambiguity: an undecided pair after escalation is COUNTED, and the
census reports interval bounds -- f_M in [related/pairs,
(related + ambiguous)/pairs] -- so Delta is reported as
[Delta_lower, Delta_upper]. With zero ambiguity the bounds coincide.

Provenance: the git state is captured at ENTRY; official runs refuse
a dirty tree, re-capture at EXIT, and refuse to write on any
mismatch. `--pilot` permits a dirty tree and writes to the pilot
artifact path with both states recorded.

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
from probe_seed_ledger import S3_SEED, assert_probe_seed_fresh

# ---------------------------------------------------------------------
# Exploration constants (NOT frozen -- this is a probe, not a stage)
# ---------------------------------------------------------------------

E_N = 300                    # Poisson mean events per reading
N_READINGS = 300             # paired readings
SEED = S3_SEED               # 40_000_201, ledger-asserted at entry
TOL = s1.DEFAULT_TOL         # 1e-8, the priced S1 operating point
TOL_ESCALATED = 1e-10        # last rung of s1.TOL_LADDER

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_s3_probe_results.json")
_PILOT_ARTIFACT = _ARTIFACT.with_name("p14_s3_pilot_results.json")


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
            n_events: int) -> tuple[float, float, float, int, int]:
    """(f_M_lower, f_M_upper, f_0, ambiguous, escalated) on one shared
    sprinkle of `n_events` events (the caller draws the Poisson count).

    A pair undecided at TOL escalates to TOL_ESCALATED; a pair still
    undecided is counted as ambiguous and enters the census as the
    interval [unrelated, related] -- never a silent False."""

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
            elif rel:
                rel_m += 1
    return (rel_m / pairs, (rel_m + ambiguous) / pairs, rel_0 / pairs,
            ambiguous, escalated)


# ---------------------------------------------------------------------
# The probe
# ---------------------------------------------------------------------


def _git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"],
                           capture_output=True, text=True).stdout.strip()
    return {"rev": rev, "dirty": bool(dirty)}


def _summary(values: np.ndarray) -> dict:
    n = len(values)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1))
    sem = sd / math.sqrt(n)
    tcrit = student_t_crit(n - 1)
    return {"mean": mean, "sd": sd, "sem": sem,
            "ci95_student_t": [mean - tcrit * sem, mean + tcrit * sem],
            "per_reading": [float(v) for v in values]}


def main() -> None:
    smoke = "--smoke" in sys.argv
    pilot = "--pilot" in sys.argv
    n_readings = 3 if smoke else N_READINGS
    e_n = 40 if smoke else E_N

    assert_probe_seed_fresh("s3_exploration")
    state_start = _git_state()
    if state_start["dirty"] and not (pilot or smoke):
        raise SystemExit(
            "S3: refusing an official run from a dirty tree "
            f"(rev {state_start['rev']}); commit first or pass --pilot.")

    rng = np.random.default_rng(SEED)
    fm_lo = np.empty(n_readings)
    fm_hi = np.empty(n_readings)
    f0 = np.empty(n_readings)
    counts = np.empty(n_readings, dtype=int)
    ambiguous = escalated = 0
    start = time.perf_counter()
    for k in range(n_readings):
        n = int(rng.poisson(e_n))
        if n < 2:
            raise SystemExit(f"S3: degenerate sprinkle n={n}")
        counts[k] = n
        fm_lo[k], fm_hi[k], f0[k], amb, esc = reading(rng, n)
        ambiguous += amb
        escalated += esc
        if (k + 1) % 10 == 0 or k == 0 or (k + 1) == n_readings:
            dt = time.perf_counter() - start
            eta = dt / (k + 1) * (n_readings - k - 1)
            print(f"  s3 {k + 1}/{n_readings}  "
                  f"delta so far {np.mean(fm_lo[:k + 1] - f0[:k + 1]):+.6f}  "
                  f"elapsed {dt / 60:.1f} min, eta {eta / 60:.1f} min",
                  flush=True)
    runtime = time.perf_counter() - start

    state_end = _git_state()
    if not (pilot or smoke):
        if state_end != state_start:
            raise SystemExit(
                "S3: refusing to write -- the tree changed during the "
                f"run ({state_start} -> {state_end}).")

    delta_lower = fm_lo - f0
    delta_upper = fm_hi - f0
    result = {
        "probe": "S3 Schwarzschild paired exploration"
                 + (" (PILOT)" if pilot else ""),
        "design": ("one sprinkle, dual census (measure identity, S2 "
                   "memo); N ~ Poisson(E_N); Delta statements are tied "
                   "to the frozen Schwarzschild coordinates and domain"),
        "domain": {"m": s1.M, "r_shell": [s1.R_MIN, s1.R_MAX],
                   "cap_half_angle": s1.CAP_HALF_ANGLE,
                   "t_extent": s1.T_EXTENT},
        "params": {"e_n": e_n, "n_readings": n_readings, "seed": SEED,
                   "tol": TOL, "tol_escalated": TOL_ESCALATED,
                   "smoke": smoke, "pilot": pilot},
        "event_counts": {"mean": float(np.mean(counts)),
                         "min": int(counts.min()),
                         "max": int(counts.max())},
        "delta_lower": _summary(delta_lower),
        "delta_upper": _summary(delta_upper),
        "bounds_note": ("delta_lower counts undecided pairs as "
                        "unrelated, delta_upper as related; identical "
                        "when ambiguity is 0"),
        "f_flat": {"mean": float(np.mean(f0)),
                   "sd": float(np.std(f0, ddof=1))},
        "f_schwarzschild_lower": {"mean": float(np.mean(fm_lo)),
                                  "sd": float(np.std(fm_lo, ddof=1))},
        "ambiguity": {"ambiguous": ambiguous, "escalated": escalated},
        "runtime_seconds": runtime,
        "code": {"start": state_start, "end": state_end},
    }

    lo, hi = result["delta_lower"], result["delta_upper"]
    print(f"\nS3: delta_lower mean {lo['mean']:+.6f} "
          f"CI95 [{lo['ci95_student_t'][0]:+.6f}, "
          f"{lo['ci95_student_t'][1]:+.6f}]  "
          f"delta_upper mean {hi['mean']:+.6f}  "
          f"ambiguous {ambiguous}, escalated {escalated}  "
          f"({runtime / 3600:.2f} h)")

    if smoke:
        print("smoke run -- artifact NOT written")
        return
    target = _PILOT_ARTIFACT if pilot else _ARTIFACT
    target.write_text(json.dumps(result, indent=2) + "\n",
                      encoding="utf-8")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
