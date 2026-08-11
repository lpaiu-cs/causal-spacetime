"""O4: auditing the S1 predicate's volume response against the O3
certified oracle -- campaign runner.

EXECUTION IS NOT YET APPROVED. The design ruling approves preparing
this freeze; the campaign itself needs a separate approval after the
freeze PR review converges and an exact-checkout preflight passes.

What this stage is, in the ruling's own words: a three-way consistency
audit of `sampler + S1 flight-time volume response + O3 oracle`. It is
NOT a Poisson-count audit -- the estimand is the deterministic
functional

    V_S1 = integral over the certified box of L_S1 dmu_3,
    L_S1(X) = [dt - T1_S1(X) - T2_S1(X)]_+,

with T_S1 the S1 solver's DETERMINISTIC point output at tol 1e-8. The
solver's `err` never enters the inference: it is a Gauss-Legendre
stopping heuristic, not a certified enclosure, so an interval built
from it would not be a conservative bound (review P1-5). Anything the
solver cannot answer -- exception, non-finite, negative, or an
undecided `causal_relation` -- aborts the run with no result.

Nothing here upgrades Paper A. A concordant result says the
instrument stack agrees with the certified continuum volume to within
the frozen margin at ONE frozen configuration; a discordant result
indicts the stack (sampler, predicate, or certification), not causal
set theory and not the S4/S5 verdicts.

Gates (familywise 0.04 + 0.01 = 0.05, so the composed G1+G2 sentence
carries simultaneous coverage >= 95%):

  G1 primary       equivalence on V_S1 inside the certified box
  G2 co-primary    leakage bound outside the box
  G3 prerequisite  `causal_relation` wrapper contract (NOT a
                   scientific verdict, never composed into the
                   3.25% sentence)

Run (after separate execution approval, from a clean exact checkout):

    python experiments/oracle/o4_volume_audit.py --preflight
    python experiments/oracle/o4_volume_audit.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import empirical_bernstein as eb  # noqa: E402
import exact_binomial as xb  # noqa: E402
import gmpy2  # noqa: E402
import o4_sizing as sz  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402
from probe_seed_ledger import (  # noqa: E402
    O4_SMOKE_SEED,
    assert_fresh_scalar,
)

_MANIFEST = _REPO / "docs" / "prereg" / "p14_o4_freeze_manifest.json"
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o4_results.json"

#: The frozen configuration. Sizes come from `o4_sizing`, which
#: derives them from the sampler's own boundaries in one path -- they
#: are re-asserted here so a desynchronised edit fails loudly.
FROZEN = {
    "tau": sz.TAU,
    "delta_g1_per_side": sz.DELTA_G1,
    "alpha_g2_per_end": sz.ALPHA_G2,
    "alpha_g3": sz.ALPHA_G3,
    "leak_budget_frac": sz.LEAK_BUDGET,
    "n_g1": sz.N_G1,
    "n_g2": sz.g2_points(),
    "g3_clusters": sz.G3_CLUSTERS,
    "g3_stress_offset": 1e-6,
    "tol": s1.DEFAULT_TOL,
    "max_calls": 80_000_000,
    "max_wall_s": 86_400.0,
}

#: Cap headroom, stated separately because the two ratios differ
#: (review item 4): calls 80M / 54.9M = 1.46x, wall 24h / 11.7h = 2.05x.
CHUNK = 4096

#: Smoke runs execute the frozen ANCHORS, so a large one would
#: pre-observe the very quantity the freeze protects -- the O3 lesson
#: ("running the frozen configuration first collapses the freeze
#: boundary"). The cap keeps the smoke estimate uninformative: the
#: sampling standard deviation of V_S1 is ~415/sqrt(n), so at n = 2000
#: it is ~9.3, more than five times the equivalence band 1.70, while
#: at the campaign's 2.62e7 it is 0.081. Smoke validates plumbing,
#: never magnitude.
SMOKE_MAX = 2_000


# ---------------------------------------------------------------- io

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _environment() -> dict:
    return {
        "python": ".".join(platform.python_version_tuple()[:2]),
        "gmpy2": gmpy2.version(),
        "mpfr": gmpy2.mpfr_version(),
        "gmp": gmpy2.mp_version(),
        "numpy": np.__version__,
    }


def verify_digests(stage: str) -> dict:
    manifest = json.loads(_MANIFEST.read_text(encoding="utf-8"))
    bad = [rel for rel, want in manifest["files"].items()
           if _sha256(_REPO / rel) != want]
    if bad:
        raise SystemExit(
            f"freeze verification ({stage}): digest mismatch on "
            f"{bad} -- refusing to run on a drifted protocol surface")
    return manifest


def verify_environment(stage: str, manifest: dict) -> None:
    env, locked = _environment(), manifest["environment"]
    drift = {k: (locked[k], env.get(k)) for k in locked
             if env.get(k) != locked[k]}
    if drift:
        raise SystemExit(
            f"freeze verification ({stage}): environment drift "
            f"{drift} -- the solver and the rounding instrument are "
            f"part of the frozen apparatus")


def verify_freeze(stage: str) -> dict:
    manifest = verify_digests(stage)
    verify_environment(stage, manifest)
    return manifest


def _git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         check=True, capture_output=True,
                         text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO,
                           check=True, capture_output=True,
                           text=True).stdout.strip()
    return {"rev": rev, "dirty": bool(dirty)}


def preflight() -> dict:
    manifest = verify_freeze("preflight")
    state = _git_state()
    if state["dirty"]:
        raise SystemExit(
            "preflight: working tree is dirty -- the campaign runs "
            "from a clean exact checkout only")
    if _ARTIFACT.exists():
        raise SystemExit(
            f"preflight: {_ARTIFACT.name} already exists -- the audit "
            f"is write-once; a rerun may not overwrite or relabel the "
            f"first observation")
    for name in ("g1_audit", "g2_leakage", "g3_wrapper"):
        assert_fresh_scalar(name)
    return {"manifest": manifest, "git": state}


def _publish_write_once(path: Path, payload: str) -> None:
    """Atomic no-clobber publication (the O3 pattern): fsynced temp in
    the same directory, then `os.link`, which fails atomically if the
    destination exists."""

    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(payload)
            f.flush()
            os.fsync(f.fileno())
        try:
            os.link(tmp, path)
        except FileExistsError:
            raise SystemExit(
                f"publish: {path.name} already exists -- the audit is "
                f"write-once and the first observation stands"
            ) from None
    finally:
        tmp.unlink(missing_ok=True)


# ----------------------------------------------------------- sampling

def _draw(rng: np.random.Generator, k: int, r_lo: float, r_hi: float,
          psi: float) -> tuple[np.ndarray, np.ndarray]:
    """`k` points distributed as the r^2 sin(theta) measure on the
    polar-cap shell: r uniform in r^3, cos(theta) uniform. phi is not
    drawn -- the configuration is axisymmetric, so L depends only on
    (r, theta) and the phi integral is the 2 pi already inside B."""

    u = rng.random(k)
    r = np.cbrt(u * (r_hi ** 3 - r_lo ** 3) + r_lo ** 3)
    v = rng.random(k)
    theta = np.arccos(v * (1.0 - math.cos(psi)) + math.cos(psi))
    return r, theta


def _flight(r1: float, r2: float, dpsi: float, tol: float) -> float:
    """One S1 flight time, fail-closed. The deterministic point output
    IS the estimand's input; `err` is deliberately discarded."""

    try:
        t, _ = s1.flight_time(r1, r2, dpsi, s1.M, tol)
    except Exception as exc:                    # noqa: BLE001
        raise SystemExit(
            f"fail-closed: flight_time({r1!r}, {r2!r}, {dpsi!r}) "
            f"raised {exc!r} -- no result is published") from None
    if not math.isfinite(t) or t < 0.0:
        raise SystemExit(
            f"fail-closed: flight_time({r1!r}, {r2!r}, {dpsi!r}) "
            f"returned {t!r} -- no result is published")
    return t


def _ell(r: float, theta: float, tol: float) -> tuple[float, float,
                                                      float]:
    """(L_S1, T1, T2) at one spatial point."""

    t1 = _flight(sz.R_IN, r, theta, tol)
    t2 = _flight(r, sz.R_OUT, theta, tol)
    ell = sz.DT - t1 - t2
    if ell > sz.DT:
        raise SystemExit(
            f"fail-closed: L_S1={ell!r} exceeds dt={sz.DT!r}, so the "
            f"empirical Bernstein precondition Z in [0, 1] fails")
    return (ell if ell > 0.0 else 0.0), t1, t2


# -------------------------------------------------------------- gates

class _Budget:
    """Frozen caps. Raising them during or after the run is forbidden;
    if they bind, the run reports the fact and the verdict degrades."""

    def __init__(self, cfg: dict) -> None:
        self.calls = 0
        self.max_calls = cfg["max_calls"]
        self.t0 = time.perf_counter()
        self.max_wall = cfg["max_wall_s"]
        self.reason = "complete"

    def spend(self, k: int) -> bool:
        self.calls += k
        if self.calls >= self.max_calls:
            self.reason = "max-calls"
            return False
        if time.perf_counter() - self.t0 >= self.max_wall:
            self.reason = "max-wall"
            return False
        return True

    @property
    def wall_s(self) -> float:
        return time.perf_counter() - self.t0


def run_g1(cfg: dict, seed: int, budget: _Budget,
           progress=None) -> dict:
    """Primary: equivalence on V_S1 inside the certified box."""

    rng = np.random.default_rng(seed)
    acc = eb.Accumulator()
    stress: list[tuple[float, float, float, float]] = []
    n_target = cfg["n_g1"]
    over_ceiling = 0
    last = 0.0
    while acc.n < n_target:
        k = min(CHUNK, n_target - acc.n)
        rs, ths = _draw(rng, k, sz.R_LO, sz.R_HI, sz.PSI_MAX)
        for r, th in zip(rs, ths, strict=True):
            ell, t1, t2 = _ell(float(r), float(th), cfg["tol"])
            if ell > sz.L_MAX_UB:
                # a diagnostic outcome, independent of the verdict:
                # S1 claims to beat the certified optical-triangle
                # minimum. Recorded, never used to clip the estimand.
                over_ceiling += 1
            acc.add(ell / sz.DT)
            if ell > 0.0 and len(stress) < cfg["g3_clusters"]:
                stress.append((float(r), float(th), t1, t2))
        if not budget.spend(2 * k):
            break
        if progress is not None and budget.wall_s - last >= 300.0:
            last = budget.wall_s
            progress("g1", acc.n, budget)
    iv = eb.interval(acc, cfg["delta_g1_per_side"])
    c_lo, c_hi = iv.rescaled(sz.SCALE)
    ident = (c_lo - sz.V_HI, c_hi - sz.V_LO)
    band = cfg["tau"] * sz.V_REF
    if ident[0] >= -band and ident[1] <= band:
        status = "concordant"
    elif ident[1] < -band or ident[0] > band:
        status = "discordant"
    else:
        status = "inconclusive"
    return {
        "status": status, "n": acc.n, "n_target": n_target,
        "mean_z": acc.mean, "var_z": acc.var,
        "half_width_z": iv.half_width,
        "v_s1_lo": c_lo, "v_s1_hi": c_hi,
        "identified_discrepancy": list(ident),
        "band_abs": band,
        "over_ceiling_points": over_ceiling,
        "power_certified": acc.n >= n_target,
        "_stress": stress,
    }


def run_g2(cfg: dict, seed: int, budget: _Budget) -> dict:
    """Co-primary: leaked volume outside the certified box.

    Both ends of the decision carry one shared error budget -- a
    Clopper-Pearson upper bound at alpha/2 on the leaking fraction and
    an empirical Bernstein lower bound at alpha/2 on E[L_out] (review
    P1-3: the v0.1 draft charged only the upper end)."""

    rng = np.random.default_rng(seed)
    acc = eb.Accumulator()
    leaks = 0
    n_target = cfg["n_g2"]
    while acc.n < n_target:
        k = min(CHUNK, n_target - acc.n)
        rs, ths = _draw(rng, 2 * k, sz.PATCH_R[0], sz.PATCH_R[1],
                        sz.PATCH_CAP)
        used = 0
        for r, th in zip(rs, ths, strict=True):
            if acc.n >= n_target:
                break
            if (sz.R_LO <= r <= sz.R_HI) and th <= sz.PSI_MAX:
                continue                       # inside the box
            ell, _, _ = _ell(float(r), float(th), cfg["tol"])
            acc.add(ell / sz.DT)
            used += 1
            if ell > 0.0:
                leaks += 1
        if not budget.spend(2 * used):
            break
    alpha = cfg["alpha_g2_per_end"]
    p_up = xb.cp_upper(leaks, acc.n, alpha)
    leak_ub = sz.DT * sz.B_OUT * p_up
    leak_lb = sz.B_OUT * sz.DT * max(
        0.0, eb.lower_bound(acc, alpha))
    budget_abs = cfg["leak_budget_frac"] * sz.V_REF
    if leak_ub <= budget_abs:
        status = "concordant"
    elif leak_lb > budget_abs:
        status = "discordant"
    else:
        status = "inconclusive"
    return {"status": status, "n": acc.n, "n_target": n_target,
            "leaking_points": leaks, "cp_upper_rate": p_up,
            "leak_upper_abs": leak_ub, "leak_lower_abs": leak_lb,
            "budget_abs": budget_abs}


def run_g3(cfg: dict, stress: list, budget: _Budget) -> dict:
    """Instrumentation prerequisite: does `causal_relation` accept
    exactly the time window that L_S1 implies?

    The independent unit is the spatial point (CLUSTER): the two
    stress times at one point share that point's flight times, so
    they are not independent Bernoulli trials (review P2). A cluster
    counts as mismatching if either stress time disagrees.

    The reported rate is over the frozen boundary-stress
    distribution -- G1's sampling measure conditioned on L_S1 > 0,
    probed at the window midpoint and just outside its upper edge --
    not over a generic 4D distribution."""

    mismatched = 0
    n = 0
    off = cfg["g3_stress_offset"]
    for r, th, t1, t2 in stress:
        lo_edge, hi_edge = t1, sz.DT - t2
        p = np.array([0.0, sz.R_IN, 0.0, 0.0])
        q = np.array([sz.DT, sz.R_OUT, 0.0, 0.0])
        bad = False
        for t_x, want in ((0.5 * (lo_edge + hi_edge), True),
                          (hi_edge + off, False)):
            x = np.array([t_x, r, th, 0.0])
            a = s1.causal_relation(p, x, s1.M, cfg["tol"])
            b = s1.causal_relation(x, q, s1.M, cfg["tol"])
            if a is None or b is None:
                raise SystemExit(
                    "fail-closed: causal_relation returned undecided "
                    "at a G3 stress point -- no result is published")
            if bool(a and b) is not want:
                bad = True
        n += 1
        mismatched += int(bad)
        if not budget.spend(4):
            break
    rate_ub = (xb.cp_upper(mismatched, n, cfg["alpha_g3"])
               if n else 1.0)
    return {"valid": mismatched == 0 and n == cfg["g3_clusters"],
            "clusters": n, "clusters_target": cfg["g3_clusters"],
            "mismatching_clusters": mismatched,
            "cp_upper_rate": rate_ub,
            "distribution": ("G1 sampling measure conditioned on "
                             "L_S1 > 0, probed at the window midpoint "
                             "and at the upper edge + 1e-6")}


def compose_verdict(g1: dict, g2: dict, g3: dict) -> str:
    """The frozen composition (review item 3). G3 is an
    instrumentation contract: its failure invalidates the stage and is
    NEVER reported as scientific discordance."""

    if not g3["valid"]:
        return "INVALID"
    if g1["status"] == "discordant" or g2["status"] == "discordant":
        return "DISCORDANT"
    if not g1["power_certified"]:
        return "INCONCLUSIVE"
    if g1["status"] == "concordant" and g2["status"] == "concordant":
        return "CONCORDANT"
    return "INCONCLUSIVE"


# --------------------------------------------------------------- main

def assemble(cfg: dict, seeds: dict, progress=None) -> dict:
    budget = _Budget(cfg)
    g1 = run_g1(cfg, seeds["g1"], budget, progress)
    stress = g1.pop("_stress")
    g2 = run_g2(cfg, seeds["g2"], budget)
    g3 = run_g3(cfg, stress, budget)
    return {
        "verdict": compose_verdict(g1, g2, g3),
        "g1": g1, "g2": g2, "g3": g3,
        "termination_reason": budget.reason,
        "calls": budget.calls, "wall_s": budget.wall_s,
        "total_sentence_abs": (cfg["tau"] + cfg["leak_budget_frac"])
        * sz.V_REF,
    }


def main() -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "--preflight", action="store_true",
        help="run every pre-execution check and observe nothing")
    parser.add_argument(
        "--smoke", type=int, default=0, metavar="N",
        help="tiny validation run on the dedicated smoke stream; "
             "never touches the campaign seeds or the artifact")
    args = parser.parse_args()

    if args.smoke:
        if args.smoke > SMOKE_MAX:
            raise SystemExit(
                f"--smoke {args.smoke} exceeds {SMOKE_MAX}: a larger "
                f"smoke run on the FROZEN anchors would pre-observe "
                f"V_S1 at a precision comparable to the equivalence "
                f"band and collapse the freeze boundary")
        cfg = FROZEN | {"n_g1": args.smoke,
                        "n_g2": max(2, args.smoke // 4),
                        "g3_clusters": max(1, args.smoke // 100),
                        "max_calls": 10 ** 9, "max_wall_s": 3600.0}
        res = assemble(cfg, {"g1": O4_SMOKE_SEED,
                             "g2": O4_SMOKE_SEED + 1,
                             "g3": O4_SMOKE_SEED + 2})
        print(json.dumps({k: v for k, v in res.items()
                          if k != "curve"}, indent=1, default=str))
        return

    checks = preflight()
    if args.preflight:
        print("preflight PASS: freeze digests, environment lock, "
              f"clean tree at {checks['git']['rev'][:9]}, no result "
              "artifact, campaign seeds fresh. Nothing was observed.")
        return

    start = _git_state()
    seeds = {"g1": assert_fresh_scalar("g1_audit"),
             "g2": assert_fresh_scalar("g2_leakage"),
             "g3": assert_fresh_scalar("g3_wrapper")}

    def show(stage: str, done: int, budget: _Budget) -> None:
        print(f"  {stage} {done:,} pts  calls={budget.calls:,}  "
              f"t={budget.wall_s:.0f}s", flush=True)

    t0 = time.perf_counter()
    res = assemble(FROZEN, seeds, progress=show)
    verify_freeze("exit")
    end = _git_state()
    if end["rev"] != start["rev"] or end["dirty"]:
        raise SystemExit(
            "the tree changed underneath the campaign -- refusing to "
            "write a result with broken lineage")

    art = {
        "stage": ("O4: S1 predicate volume-response audit against the "
                  "O3 certified oracle, frozen anchors (12, 18, 8.5)M"),
        "frozen_config": FROZEN,
        "seeds": seeds,
        "geometry": sz.summary()["geometry"],
        "oracle": {"v_lo": sz.V_LO, "v_hi": sz.V_HI,
                   "v_ref": sz.V_REF},
        "environment": _environment(),
        "host": {"machine": platform.machine(),
                 "system": platform.system()},
        "code": {"start": start, "end": end},
        "result": res,
        "total_wall_s": time.perf_counter() - t0,
    }
    _publish_write_once(_ARTIFACT,
                        json.dumps(art, ensure_ascii=False,
                                   indent=1) + "\n")
    print(f"verdict: {res['verdict']}  "
          f"(G1 {res['g1']['status']}, G2 {res['g2']['status']}, "
          f"G3 {'valid' if res['g3']['valid'] else 'INVALID'})")
    print(f"artifact: {_ARTIFACT}")


if __name__ == "__main__":
    main()
