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

Provenance: the git state is captured at ENTRY; the run refuses a
dirty tree, re-captures at EXIT, and refuses to write on any
mismatch. The historical fixed-N pilot is preserved separately with a
provenance sidecar; its seed 40_000_201 is ledger-spent (observed).
The official artifact was observed on stream 40_000_231 (allocated
fresh, moved to OBSERVED in the results commit), so any rerun of the
official path is a deterministic REPLAY of that stream -- seed
freshness and replay are distinct ledger operations (PR review R2).
`--smoke` draws only the dedicated observed smoke stream 40_000_221
(the first official allocation, demoted after its smoke output was
observed during review-R2 validation), so validation can never spend
an official seed (PR review R3).

Claim boundary (PR review R2): the paired estimand Delta = f_M - f_0
on one point set is C1-class counterfactual sensitivity in the P14
taxonomy. It sizes a C1-style paired confirmation with an
operationally anchored margin (the epsilon_Delta pattern), which
needs no diamond-volume oracle. It is NOT a C2-class single-poset
discrimination; that would need independently sprinkled
Schwarzschild/flat arms and an unpaired protocol, designed
separately. The artifact stores the per-reading raw f arrays so such
a design (and unpaired power models) can be sized from this
exploration.

Exploration only -- no preregistered gate.
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
from probe_seed_ledger import (
    FRESH_PROBE_SCALARS,
    S3_SEED,
    S3_SMOKE_SEED,
    replay_scalar,
)

# ---------------------------------------------------------------------
# Exploration constants (NOT frozen -- this is a probe, not a stage)
# ---------------------------------------------------------------------

E_N = 300                    # Poisson mean events per reading
N_READINGS = 300             # paired readings
SEED = S3_SEED               # 40_000_231, OBSERVED (official artifact
                             # exists) -- reruns are deterministic
                             # replays of that stream, never fresh
SMOKE_SEED = S3_SMOKE_SEED   # 40_000_221, the observed smoke stream --
                             # smoke may never draw the official seed
TOL = s1.DEFAULT_TOL         # 1e-8, the priced S1 operating point
TOL_ESCALATED = 1e-10        # last rung of s1.TOL_LADDER

_REPO = Path(__file__).resolve().parents[2]
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s3_probe_results.json"


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
    """The state of THIS script's repository, pinned via cwd=_REPO --
    never the process working directory, where a different (or no)
    repository could fabricate a clean lineage -- and with check=True
    so a failed git call refuses instead of recording rev='' clean
    (PR review)."""

    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         capture_output=True, text=True,
                         check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO,
                           capture_output=True, text=True,
                           check=True).stdout.strip()
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
    n_readings = 3 if smoke else N_READINGS
    e_n = 40 if smoke else E_N

    if smoke:
        seed = SMOKE_SEED
    else:
        seed = replay_scalar("s3_exploration")
        if SEED != seed:
            raise SystemExit("S3: SEED drifted from the observed "
                             "stream; a rerun must be a replay of it.")
    state_start = _git_state()
    if state_start["dirty"] and not smoke:
        raise SystemExit(
            "S3: refusing an official run from a dirty tree "
            f"(rev {state_start['rev']}); commit first.")

    rng = np.random.default_rng(seed)
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
    if not smoke and state_end != state_start:
        raise SystemExit(
            "S3: refusing to write -- the tree changed during the "
            f"run ({state_start} -> {state_end}).")

    delta_lower = fm_lo - f0
    delta_upper = fm_hi - f0
    lo_sum = _summary(delta_lower)
    hi_sum = _summary(delta_upper)
    run_kind = ("fresh_observation"
                if "s3_exploration" in FRESH_PROBE_SCALARS
                else "replay")
    result = {
        "probe": "S3 Schwarzschild paired exploration",
        "run_kind": run_kind,
        "replay_of": (None if run_kind == "fresh_observation" else
                      "observed stream 40000231; the fresh observation "
                      "is the artifact committed at fe3f353 -- a replay "
                      "reproduces its readings byte-identically and "
                      "must not silently replace it (contract test)"),
        "design": ("one sprinkle, dual census (measure identity, S2 "
                   "memo); N ~ Poisson(E_N); Delta statements are tied "
                   "to the frozen Schwarzschild coordinates and domain"),
        "claim_class": ("C1 counterfactual sensitivity (paired, same "
                        "point set); NOT a C2 single-poset "
                        "discrimination -- raw per-reading f arrays "
                        "are stored so an unpaired C2-style design "
                        "can be sized separately"),
        "domain": {"m": s1.M, "r_shell": [s1.R_MIN, s1.R_MAX],
                   "cap_half_angle": s1.CAP_HALF_ANGLE,
                   "t_extent": s1.T_EXTENT},
        "params": {"e_n": e_n, "n_readings": n_readings, "seed": seed,
                   "tol": TOL, "tol_escalated": TOL_ESCALATED,
                   "smoke": smoke},
        "event_counts": {"mean": float(np.mean(counts)),
                         "min": int(counts.min()),
                         "max": int(counts.max()),
                         "per_reading": [int(c) for c in counts]},
        "delta_lower": lo_sum,
        "delta_upper": hi_sum,
        # Each endpoint carries 97.5% one-sided coverage, so taking
        # the outer endpoints is a conservative joint 95% interval by
        # the union bound.
        "identified_ci95": [lo_sum["ci95_student_t"][0],
                            hi_sum["ci95_student_t"][1]],
        "bounds_note": ("delta_lower counts undecided pairs as "
                        "unrelated, delta_upper as related; identical "
                        "when ambiguity is 0; identified_ci95 spans "
                        "the lower CI's bottom to the upper CI's top"),
        "f_flat": {"mean": float(np.mean(f0)),
                   "sd": float(np.std(f0, ddof=1)),
                   "per_reading": [float(v) for v in f0]},
        "f_schwarzschild_lower": {"mean": float(np.mean(fm_lo)),
                                  "sd": float(np.std(fm_lo, ddof=1)),
                                  "per_reading": [float(v) for v in fm_lo]},
        "f_schwarzschild_upper": {"mean": float(np.mean(fm_hi)),
                                  "sd": float(np.std(fm_hi, ddof=1)),
                                  "per_reading": [float(v) for v in fm_hi]},
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
    _ARTIFACT.write_text(json.dumps(result, indent=2) + "\n",
                         encoding="utf-8")
    print(f"wrote {_ARTIFACT}")


if __name__ == "__main__":
    main()
