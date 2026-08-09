"""W1: Weyl/Ricci channel-decomposition exploration probe (W-track).

Same Brinkmann family, same slab, same measure -- det g = -1 holds for
EVERY profile H(u, x, y), so one sprinkle serves every arm. Three
curved channels are censused against the same-point flat census at the
P3-C frozen operating point (aniso-a1.0):

  vacuum  H = +(w^2/2)(x^2 - y^2)  Weyl-only   (type N, vacuum)
  ricci   H = -(w^2/2)(x^2 + y^2)  Ricci-only  (type O, null dust --
                                    both transverse directions focus,
                                    the attractive/NEC-consistent sign)
  mixed   H = -(w^2/2) y^2         equal Weyl and Ricci parts
                                    (c_x - c_y = c_x + c_y = w^2)

per-leg channel signature: vacuum (defocus, focus), ricci (focus,
focus), mixed (flat, focus).

The vacuum and flat censuses use the FROZEN predicate unchanged. The
new channel predicate is needed only for ricci and mixed; it is
composed from the frozen chain's own per-leg primitives
(cost_defocusing / cost_focusing / _leg_error_bound / _series), keeps
the escalate-rather-than-guess contract including the R9.1 rule
(Decimal separations formed from STORED coordinates), and must
reproduce the frozen verdicts bit-identically on the vacuum signature
-- sealed by tests/test_w1_channel_probe.py. The frozen sources are
digest-pinned, which is WHY this lives in a new module and imports
their private per-leg names instead of editing them.

Exploration only -- no preregistered gate. Output sizes the channel
story: does the census respond to Weyl-sourced and Ricci-sourced cone
bending identically, additively, or with structure (the P3-E 3a sign
reversal is the standing hint that channels are not interchangeable)?
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
from seed_windows import (
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_point_seeds_fresh,
)

# ---------------------------------------------------------------------
# Exploration constants (NOT frozen -- this is a probe, not a stage)
# ---------------------------------------------------------------------

N_READINGS = 300
SEED = 40_000_211            # fresh vs every ledger, see SPENT_SCALARS

#: Every scalar seed the program has spent (the S3 probe's ledger plus
#: S3's own 40_000_201).
SPENT_SCALARS = frozenset({
    777, 778, 779, 780, 781,
    20_260_808, 20_260_811, 20_260_812, 20_260_813, 20_260_814,
    20_260_821, 20_260_822, 20_260_823, 20_260_824,
    20_260_851, 20_260_852,
    40_000_061, 40_000_071, 40_000_072, 40_000_101, 40_000_201,
})

#: Channel signatures: (x-leg, y-leg), each "defocus"|"focus"|"flat".
VACUUM = ("defocus", "focus")
RICCI = ("focus", "focus")
MIXED = ("flat", "focus")

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_w1_channel_results.json")


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
    the undecided-is-never-a-silent-False contract as in the frozen
    census: an ambiguous pair is counted and left unrelated, bounding
    the bias by ambiguous/pairs."""

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
    rev = subprocess.run(["git", "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"],
                           capture_output=True, text=True).stdout.strip()
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

    assert_point_seeds_fresh(
        {"w1_exploration": SEED}, SPENT_SCALARS,
        P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,), "W1")

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

    deltas = {name: f[name] - f["flat"]
              for name in ("vacuum", "ricci", "mixed")}
    contrasts = {"vacuum_minus_ricci": deltas["vacuum"] - deltas["ricci"],
                 "vacuum_minus_mixed": deltas["vacuum"] - deltas["mixed"],
                 "ricci_minus_mixed": deltas["ricci"] - deltas["mixed"]}

    result = {
        "probe": "W1 Weyl/Ricci channel decomposition",
        "design": ("one sprinkle per reading (det g = -1 for every H), "
                   "four censuses on the same points; vacuum and flat "
                   "by the frozen predicate, ricci and mixed by the "
                   "channel predicate composed from frozen legs"),
        "operating_point": {"label": label, "w": w,
                            "slab": [du, dv, dx, dy], "e_n": p3c.E_N},
        "channels": {"vacuum": list(VACUUM), "ricci": list(RICCI),
                     "mixed": list(MIXED)},
        "params": {"n_readings": n_readings, "seed": SEED,
                   "smoke": smoke},
        "f_means": {k: float(np.mean(v)) for k, v in f.items()},
        "delta": {k: _summary(v) for k, v in deltas.items()},
        "contrasts": {k: _summary(v) for k, v in contrasts.items()},
        "ambiguity": {"ambiguous": amb, "escalated": esc},
        "runtime_seconds": runtime,
        "code": _git_state(),
    }

    print("\nW1:")
    for name in ("vacuum", "ricci", "mixed"):
        s = result["delta"][name]
        print(f"  {name:7s} delta {s['mean']:+.6f}  sd {s['sd']:.6f}  "
              f"CI95 [{s['ci95_student_t'][0]:+.6f}, "
              f"{s['ci95_student_t'][1]:+.6f}]")
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
