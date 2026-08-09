"""W1: Weyl/Ricci channel-decomposition exploration probe (W-track).

Same Brinkmann family, same slab, same measure -- det g = -1 holds for
EVERY profile H(u, x, y), so one sprinkle serves every arm. Three
curved channels are censused against the same-point flat census at the
P3-C frozen operating point (aniso-a1.0). The legs integrate
x'' = +w^2 x (defocus) and y'' = -w^2 y (focus), and with the
transverse geodesic equation a'' = (1/2) dH/da the matching profiles
are (PR review R1 -- NOT w^2/2):

  vacuum  H = +w^2 (x^2 - y^2)   x defocus, y focus
  ricci   H = -w^2 (x^2 + y^2)   x focus,   y focus
  mixed   H = -w^2 y^2           x flat,    y focus

Channel normalization (Brinkmann: R_uiuj = -1/2 d_i d_j H,
R_uu = -1/2 Lap H; cf. Harte & Drivas, PRD 85, 124039):

  vacuum  tidal eigenvalues (-w^2, +w^2): R_uu = 0, pure Weyl
  ricci   tidal eigenvalues (+w^2, +w^2): R_uu = +2 w^2 > 0 (NEC,
          null dust), zero Weyl
  mixed   tidal eigenvalues (0, +w^2): R_uu = +w^2 and Weyl
          eigenvalues (+-)w^2/2 -- EXACTLY half of each pure arm,
          which is what licenses the midpoint contrast below

Both pure arms share the same tidal eigenvalue magnitude w^2; every
between-channel comparison ("ricci reads stronger") is a statement AT
THAT NORMALIZATION, not an invariant ranking. The decomposition is
pinned by tests/test_w1_channel_probe.py via `channel_curvature`.

The vacuum and flat censuses use the FROZEN predicate unchanged. The
new channel predicate is needed only for ricci and mixed; it is
composed from the frozen chain's own per-leg primitives
(cost_defocusing / cost_focusing / _leg_error_bound / _series), keeps
the escalate-rather-than-guess contract including the R9.1 rule
(Decimal separations formed from STORED coordinates), and must
reproduce the frozen verdicts bit-identically on the vacuum signature
-- sealed by tests. The frozen sources are digest-pinned, which is WHY
this lives in a new module and imports their private per-leg names
instead of editing them.

Provenance: git state captured at entry, official runs refuse a dirty
tree and refuse to write if the exit state differs. A nonzero
ambiguity count anywhere refuses artifact promotion (the probe's f
would otherwise silently treat undecided as False).

Seed semantics (PR review R2): W1's stream 40_000_211 is OBSERVED --
its results exist. Any rerun of this script is therefore a
deterministic REPLAY of that observed stream (the committed artifact
is one such replay, regenerated for corrected entry/exit provenance
with byte-identical readings), never a new observation. The ledger
entry point is `replay_scalar`, not a freshness assert.

Exploration only -- no preregistered gate.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np
import p14_plane_wave as pw
import p14_probe_p3c as p3c
from p14_plane_wave import Relation, Slab, arms, sprinkle
from p14_probe_p1 import relation_census
from p14_probe_p2 import student_t_crit
from probe_seed_ledger import FRESH_PROBE_SCALARS, W1_SEED, replay_scalar

# ---------------------------------------------------------------------
# Exploration constants (NOT frozen -- this is a probe, not a stage)
# ---------------------------------------------------------------------

N_READINGS = 300
SEED = W1_SEED               # 40_000_211, OBSERVED -- reruns replay it

#: Channel signatures: (x-leg, y-leg), each "defocus"|"focus"|"flat".
VACUUM = ("defocus", "focus")
RICCI = ("focus", "focus")
MIXED = ("flat", "focus")

#: H-coefficient per leg kind: H = c_x x^2 + c_y y^2 with
#: a'' = (1/2) dH/da = c_a a.
_H_COEFF = {"defocus": +1.0, "focus": -1.0, "flat": 0.0}

_REPO = Path(__file__).resolve().parents[2]
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_w1_channel_results.json"


def channel_curvature(channel: tuple[str, str], w: float) -> dict:
    """The exact curvature decomposition implied by a channel
    signature, in the Brinkmann conventions of the module docstring
    (R_uiuj = -1/2 d_i d_j H): tidal eigenvalues (-c_x, -c_y) w^2,
    R_uu = -(c_x + c_y) w^2, Weyl eigenvalues -+ (c_x - c_y)/2 w^2."""

    cx = _H_COEFF[channel[0]] * w * w
    cy = _H_COEFF[channel[1]] * w * w
    return {"h_coefficients": (cx, cy),
            "tidal_eigenvalues": (-cx, -cy),
            "ricci_uu": -(cx + cy),
            "weyl_eigenvalues": ((cy - cx) / 2.0, (cx - cy) / 2.0)}


# ---------------------------------------------------------------------
# The channel predicate: frozen per-leg primitives, new composition
# ---------------------------------------------------------------------


def _leg_cost(kind: str, s: float, ap: float, aq: float,
              w: float) -> float:
    if kind == "flat":
        return pw._leg_cost_flat(s, ap, aq)
    if kind == "defocus":
        return pw.cost_defocusing(s, ap, aq, w)
    return pw.cost_focusing(s, ap, aq, w)


def _leg_err(kind: str, s: float, ap: float, aq: float,
             w: float) -> float:
    if kind == "flat":
        return pw._OPS_PER_LEG * pw._EPS * abs(pw._leg_cost_flat(s, ap, aq))
    return pw._leg_error_bound(s, ap, aq, w, focusing=(kind == "focus"))


def _exact_leg(kind: str, ds: Decimal, ap: float, aq: float,
               w: float) -> Decimal:
    """One leg of the exact escalation cost, mirroring the per-leg
    structure of the frozen `_exact_cost` (alternating series = the
    focusing leg)."""

    dap, daq, dw = pw._dec(ap), pw._dec(aq), pw._dec(w)
    if kind == "flat" or dw == 0:
        return (daq - dap) ** 2 / (2 * ds)
    ws = dw * ds
    alt = kind == "focus"
    return (dw / (2 * pw._series(ws, odd=True, alternating=alt))) * (
        (dap * dap + daq * daq) * pw._series(ws, odd=False, alternating=alt)
        - 2 * dap * daq)


def channel_relation(channel: tuple[str, str], p: np.ndarray,
                     q: np.ndarray, w: float) -> Relation:
    """The frozen `causal_relation`, generalized to a per-leg channel
    signature. Identical float path, error model, and escalation
    discipline; on the vacuum signature it must agree bit-identically
    (test-sealed)."""

    xk, yk = channel
    du = float(q[0]) - float(p[0])
    dv = float(q[1]) - float(p[1])

    if du < 0.0:
        return Relation(False, -math.inf, -math.inf, False)
    if du == 0.0:
        # R6.3: at zero u-separation ds^2 = dx^2 + dy^2 for EVERY
        # profile (H multiplies du^2), so the frozen rule is
        # channel-independent.
        if (float(q[2]) == float(p[2]) and float(q[3]) == float(p[3])
                and dv < 0.0):
            return Relation(True, -dv, -dv, False)
        return Relation(False, -math.inf, -math.inf, False)

    if w > 0.0 and "focus" in channel and du >= pw.conjugate_du(w):
        raise ValueError(
            f"pair straddles a conjugate point: du={du} reaches "
            f"{pw.conjugate_du(w)} at w={w}; the probe's slab must "
            "forbid this by construction.")

    xp, xq = float(p[2]), float(q[2])
    yp, yq = float(p[3]), float(q[3])
    cost = _leg_cost(xk, du, xp, xq, w) + _leg_cost(yk, du, yp, yq, w)
    margin = -(dv + cost)
    error = (_leg_err(xk, du, xp, xq, w) + _leg_err(yk, du, yp, yq, w)
             + pw._OPS_PER_LEG * pw._EPS * abs(dv))

    verdict = pw._verdict(margin, error, escalated=False)
    if not verdict.ambiguous:
        return verdict

    # Escalate rather than guess, with the R9.1 rule: the separations
    # are rebuilt HERE as Decimal differences of the stored
    # coordinates.
    getcontext().prec = pw.ESCALATED_PRECISION
    du_exact = pw._dec(q[0]) - pw._dec(p[0])
    dv_exact = pw._dec(q[1]) - pw._dec(p[1])
    exact_margin = -(dv_exact + _exact_leg(xk, du_exact, xp, xq, w)
                     + _exact_leg(yk, du_exact, yp, yq, w))
    exact_error = Decimal(10) ** (-(pw.ESCALATED_PRECISION - 10))
    scale = abs(exact_margin) + abs(dv_exact) + pw._dec(abs(cost)) + 1
    return pw._verdict(float(exact_margin), float(exact_error * scale),
                       escalated=True)


def channel_census(channel: tuple[str, str], points: np.ndarray,
                   w: float) -> tuple[float, int, int]:
    """(f, ambiguous, escalated) for one channel on one sprinkle --
    an ambiguous pair is counted and left unrelated; `main` refuses
    artifact promotion on any nonzero count."""

    order = np.argsort(points[:, 0], kind="stable")
    pts = points[order]
    n = len(pts)
    related = ambiguous = escalated = 0
    for i in range(n):
        for j in range(i + 1, n):
            rel = channel_relation(channel, pts[i], pts[j], w)
            if rel.escalated:
                escalated += 1
            if rel.related is None:
                ambiguous += 1
            elif rel.related:
                related += 1
    return related / (n * (n - 1) / 2), ambiguous, escalated


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


def _summary(delta: np.ndarray) -> dict:
    n = len(delta)
    mean = float(np.mean(delta))
    sd = float(np.std(delta, ddof=1))
    sem = sd / math.sqrt(n)
    tcrit = student_t_crit(n - 1)
    return {"mean": mean, "sd": sd, "sem": sem,
            "ci95_student_t": [mean - tcrit * sem, mean + tcrit * sem],
            "per_reading": [float(d) for d in delta]}


def main() -> None:
    smoke = "--smoke" in sys.argv
    n_readings = 5 if smoke else N_READINGS

    if SEED != replay_scalar("w1_exploration"):
        raise SystemExit("W1: SEED drifted from the observed stream; "
                         "a rerun must be a replay of it.")
    state_start = _git_state()
    if state_start["dirty"] and not smoke:
        raise SystemExit(
            "W1: refusing an official run from a dirty tree "
            f"(rev {state_start['rev']}); commit first.")

    label, w, du, dv, dx, dy = p3c.POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    rho = p3c.E_N / slab.coordinate_volume
    curved, flat = arms(slab, w)
    assert slab.du < pw.conjugate_du(w)

    rng = np.random.default_rng(SEED)
    f = {"flat": np.empty(n_readings), "vacuum": np.empty(n_readings),
         "ricci": np.empty(n_readings), "mixed": np.empty(n_readings)}
    amb = {k: 0 for k in ("flat", "vacuum", "ricci", "mixed")}
    esc = {k: 0 for k in amb}
    start = time.perf_counter()
    for k in range(n_readings):
        pts = sprinkle(flat, rho, rng)
        n = len(pts)
        pairs = n * (n - 1) / 2
        cen0 = relation_census(flat, pts)
        cenv = relation_census(curved, pts)
        f["flat"][k] = cen0.related.sum() / pairs
        f["vacuum"][k] = cenv.related.sum() / pairs
        amb["flat"] += cen0.ambiguous
        esc["flat"] += cen0.escalated
        amb["vacuum"] += cenv.ambiguous
        esc["vacuum"] += cenv.escalated
        for name, channel in (("ricci", RICCI), ("mixed", MIXED)):
            f[name][k], a, e = channel_census(channel, pts, w)
            amb[name] += a
            esc[name] += e
        if (k + 1) % 25 == 0 or k == 0 or (k + 1) == n_readings:
            t = time.perf_counter() - start
            print(f"  w1 {k + 1}/{n_readings}  "
                  f"dv {np.mean(f['vacuum'][:k + 1] - f['flat'][:k + 1]):+.6f}  "
                  f"dr {np.mean(f['ricci'][:k + 1] - f['flat'][:k + 1]):+.6f}  "
                  f"dm {np.mean(f['mixed'][:k + 1] - f['flat'][:k + 1]):+.6f}  "
                  f"({t / 60:.1f} min)", flush=True)
    runtime = time.perf_counter() - start

    state_end = _git_state()
    if not smoke and state_end != state_start:
        raise SystemExit(
            "W1: refusing to write -- the tree changed during the run "
            f"({state_start} -> {state_end}).")
    if sum(amb.values()) > 0 and not smoke:
        raise SystemExit(
            f"W1: refusing artifact promotion -- nonzero ambiguity "
            f"{amb}; the census would silently treat undecided as "
            "False. Rerun with interval bounds before promoting.")

    deltas = {name: f[name] - f["flat"]
              for name in ("vacuum", "ricci", "mixed")}
    contrasts = {"vacuum_minus_ricci": deltas["vacuum"] - deltas["ricci"],
                 "vacuum_minus_mixed": deltas["vacuum"] - deltas["mixed"],
                 "ricci_minus_mixed": deltas["ricci"] - deltas["mixed"],
                 "midpoint_residual": (deltas["mixed"]
                                       - 0.5 * (deltas["vacuum"]
                                                + deltas["ricci"]))}

    run_kind = ("fresh_observation"
                if "w1_exploration" in FRESH_PROBE_SCALARS
                else "replay")
    result = {
        "probe": "W1 Weyl/Ricci channel decomposition",
        "run_kind": run_kind,
        "replay_of": (None if run_kind == "fresh_observation" else
                      "observed stream 40000211; readings are "
                      "byte-identical to the first full observation "
                      "(the superseded pre-provenance run) -- the "
                      "deterministic-replay label the fresh/replay "
                      "ledger split requires at the artifact boundary"),
        "design": ("one sprinkle per reading (det g = -1 for every H), "
                   "four censuses on the same points; vacuum and flat "
                   "by the frozen predicate, ricci and mixed by the "
                   "channel predicate composed from frozen legs"),
        "operating_point": {"label": label, "w": w,
                            "slab": [du, dv, dx, dy], "e_n": p3c.E_N},
        "channels": {name: {"signature": list(sig),
                            **{k: list(v) if isinstance(v, tuple) else v
                               for k, v in
                               channel_curvature(sig, w).items()}}
                     for name, sig in (("vacuum", VACUUM),
                                       ("ricci", RICCI),
                                       ("mixed", MIXED))},
        "normalization_note": ("both pure arms share tidal eigenvalue "
                               "magnitude w^2; between-channel "
                               "comparisons hold at that "
                               "normalization, not invariantly"),
        "params": {"n_readings": n_readings, "seed": SEED,
                   "smoke": smoke},
        "f_means": {k: float(np.mean(v)) for k, v in f.items()},
        "delta": {k: _summary(v) for k, v in deltas.items()},
        "contrasts": {k: _summary(v) for k, v in contrasts.items()},
        "ambiguity": {"ambiguous": amb, "escalated": esc},
        "runtime_seconds": runtime,
        "code": {"start": state_start, "end": state_end},
    }

    print("\nW1:")
    for name in ("vacuum", "ricci", "mixed"):
        s = result["delta"][name]
        print(f"  {name:7s} delta {s['mean']:+.6f}  sd {s['sd']:.6f}  "
              f"CI95 [{s['ci95_student_t'][0]:+.6f}, "
              f"{s['ci95_student_t'][1]:+.6f}]")
    r = result["contrasts"]["midpoint_residual"]
    print(f"  midpoint residual {r['mean']:+.7f}  "
          f"CI95 [{r['ci95_student_t'][0]:+.7f}, "
          f"{r['ci95_student_t'][1]:+.7f}]")
    print(f"  ambiguity {amb}  escalated {esc}  "
          f"({runtime / 60:.1f} min)")

    if smoke:
        print("smoke run -- artifact NOT written")
        return
    _ARTIFACT.write_text(json.dumps(result, indent=2) + "\n",
                         encoding="utf-8")
    print(f"wrote {_ARTIFACT}")


if __name__ == "__main__":
    main()
