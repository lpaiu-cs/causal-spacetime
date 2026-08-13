"""O4 G3 abort: a REPLAY-ONLY diagnostic. It observes nothing.

The O4 campaign drawn from the freeze `1eb9461` completed G1's points
and G2's call, entered G3, and stopped on the frozen fail-closed path
with no verdict (`docs/prereg/p14_o4_incident.json`). The frozen runner
recorded no failure coordinates, so this module reproduces the stress
points and reads out what the predicate actually returned at each one.

WHAT THIS IS NOT. A replay is a reproduction, never an observation. No
value produced here may be promoted into a scientific result, and
nothing here changes a frozen rule: the G3 redesign happens in a
separate prereg re-opening, not as a consequence of these numbers. The
module therefore computes NO G1 statistic, NO G2 leakage, and NO
verdict, and its artifact carries a schema disjoint from the results
and gate schemas (`run_kind: replay`).

Scope, fixed by the PI's boundary list:

  * `replay_scalar("o4_aborted_g1")` only. The fresh pool is empty and
    stays empty; `40,000,301` stays unallocated; the remote reservation
    ref is read by nothing here and is neither claimed nor deleted.
  * The G1 prefix is reproduced only as far as the ORIGINAL G3 stress
    points required -- the leading `g3_clusters` accepted points -- and
    not one draw further. G2 is not re-run at all.
  * The first undecided is recorded in full: cluster index, (r, theta),
    probe kind, leg, `dt`, `t_min`, `err`, `|dt - t_min|`, and the
    decision margin `|dt - t_min| - err`.
  * Per-cause counts over the whole fixed stress set are reported when
    the run covers it. They are DIAGNOSTIC FREQUENCIES over one frozen
    stress set, not estimates of anything.

Interleaving is exact, not an approximation. The frozen runner builds
the stress list inside G1 and probes it later in G3, but each cluster's
four probe calls are a deterministic function of that cluster's own
`(r, theta, T1, T2)` alone, so probing a cluster as soon as it is
accepted produces byte-identical inputs to `causal_relation`. Doing it
this way is what lets a bounded `--clusters` reach the first undecided
without paying for the whole prefix.

Unlike G3, this walk does NOT stop at the first undecided. That is the
point: the frozen runner could only ever see the first one.

Run (no artifact is written unless `--out` is given):

    python experiments/oracle/o4_replay_diagnostic.py --clusters 2000
    python experiments/oracle/o4_replay_diagnostic.py \
        --out docs/prereg/p14_o4_replay_diagnostic.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_sizing as sz  # noqa: E402
import o4_volume_audit as o4  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402
from probe_seed_ledger import (  # noqa: E402
    FRESH_PROBE_SCALARS,
    replay_scalar,
    spent_scalars,
)

#: The snapshot of the freeze manifest AS EXECUTED. The current
#: manifest has moved on (the abort record retired the campaign
#: scalars, which re-pinned the ledger's digest), so a replay must
#: check itself against the executed bytes, not against HEAD's.
_EXECUTED = (_REPO / "docs" / "prereg"
             / "p14_o4_executed_freeze_manifest.json")
_INCIDENT = _REPO / "docs" / "prereg" / "p14_o4_incident.json"

#: The one file allowed to differ from the executed snapshot, and the
#: reason it may. `probe_seed_ledger` is a name-to-integer table; the
#: only thing it contributes to this computation is the seed, and that
#: value is pinned below and checked, so its drift cannot reach a
#: number. Every other file in the snapshot must match exactly.
_LEDGER_REL = "experiments/positive_control/probe_seed_ledger.py"

#: The replayed stream, pinned here so the ledger is verified rather
#: than trusted.
_STREAM = "o4_aborted_g1"
_STREAM_SEED = 40_000_281

#: Withdrawn unspent by the O4 R1 review and left that way by the
#: abort. Nothing in this module may draw from it; the check is here so
#: a future edit that reaches for "a spare seed" fails loudly.
_WITHDRAWN_UNSPENT = 40_000_301

RUN_KIND = "replay"
REPLAY_OF = "PR #70 incident"

MIDPOINT, OUTSIDE = "midpoint", "outside"
LEG_PX, LEG_XQ = "p->x", "x->q"

CAUSE_MIDPOINT = "midpoint-in-error-band"
CAUSE_OUTSIDE = "outside-in-error-band"
CAUSE_MISMATCH = "boolean-mismatch"
CAUSES = (CAUSE_MIDPOINT, CAUSE_OUTSIDE, CAUSE_MISMATCH)


# ------------------------------------------------------- replay surface

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_replay_surface() -> dict:
    """The replayed code must BE the executed code.

    Every file of the executed snapshot is compared byte for byte,
    with exactly one declared exception whose contribution is checked
    directly rather than waved through."""

    if not _EXECUTED.exists():
        raise SystemExit(
            f"replay: {_EXECUTED.name} is missing -- without the "
            f"executed snapshot there is nothing to reproduce against")
    snapshot = json.loads(_EXECUTED.read_text(encoding="utf-8"))
    drift = [rel for rel, want in snapshot["files"].items()
             if rel != _LEDGER_REL and _sha256(_REPO / rel) != want]
    if drift:
        raise SystemExit(
            f"replay: {drift} differ from the executed freeze -- these "
            f"files determine the replayed values, so this would be a "
            f"different computation, not a reproduction")

    env, locked = o4._environment(), snapshot["environment"]
    env_drift = {k: (locked[k], env.get(k)) for k in locked
                 if env.get(k) != locked[k]}
    if env_drift:
        raise SystemExit(
            f"replay: environment drift {env_drift} -- the solver and "
            f"the rounding instrument are part of the apparatus whose "
            f"output is being reproduced")

    seed = replay_scalar(_STREAM)
    if seed != _STREAM_SEED:
        raise SystemExit(
            f"replay: the ledger maps {_STREAM!r} to {seed}, not "
            f"{_STREAM_SEED} -- the one thing the ledger contributes "
            f"to this replay does not match the incident record")
    if FRESH_PROBE_SCALARS:
        raise SystemExit(
            f"replay: the fresh pool is not empty "
            f"({sorted(FRESH_PROBE_SCALARS)}) -- a diagnostic must not "
            f"run while an unobserved allocation is live, and it must "
            f"never draw from one")
    if _WITHDRAWN_UNSPENT in spent_scalars():
        raise SystemExit(
            f"replay: {_WITHDRAWN_UNSPENT} has been spent -- it was "
            f"withdrawn unspent and this diagnostic must leave it so")
    return {
        "executed_snapshot": _EXECUTED.relative_to(_REPO).as_posix(),
        "executed_snapshot_sha256": _sha256(_EXECUTED),
        "files_verified": sorted(
            r for r in snapshot["files"] if r != _LEDGER_REL),
        "declared_exception": {
            "path": _LEDGER_REL,
            "why_inert": ("a name-to-integer table; its only "
                          "contribution here is the replayed seed, "
                          "which is pinned in this module and checked"),
        },
        "environment": env,
    }


# -------------------------------------------------- the G1 prefix again

def _dpsi(p: np.ndarray, q: np.ndarray) -> float:
    """The separation angle exactly as `causal_relation` forms it.

    Not `theta`: the predicate recovers the angle through
    `acos(cos theta)`, which need not round back to the `theta` that
    `_ell` handed to `flight_time`. The difference is at most an ulp
    and is recorded rather than assumed away."""

    cosang = (math.sin(p[2]) * math.sin(q[2]) * math.cos(p[3] - q[3])
              + math.cos(p[2]) * math.cos(q[2]))
    return math.acos(max(-1.0, min(1.0, cosang)))


def stress_clusters(seed: int, clusters: int, tol: float,
                    progress=None):
    """Yield `(index, r, theta, T1, T2)` for the leading accepted
    points of the replayed G1 stream, in the frozen order.

    This mirrors `run_g1`'s draw loop exactly -- same generator, same
    chunking against the frozen `n_g1`, same acceptance test -- while
    accumulating no statistic of any kind. `_ell` is the frozen
    function, called for its side-effect-free value only."""

    rng = np.random.default_rng(seed)
    n_target = o4.FROZEN["n_g1"]
    seen = found = 0
    t0 = time.perf_counter()
    last = 0.0
    while found < clusters and seen < n_target:
        k = min(o4.CHUNK, n_target - seen)
        rs, ths = o4._draw(rng, k, sz.R_LO, sz.R_HI, sz.PSI_MAX)
        for r, th in zip(rs, ths, strict=True):
            seen += 1
            try:
                ell, t1, t2 = o4._ell(float(r), float(th), tol)
            except SystemExit as exc:
                raise SystemExit(
                    f"REPLAY DIVERGENCE at prefix point {seen}: the "
                    f"frozen fail-closed path fired where the campaign "
                    f"passed ({exc}) -- the stream is not reproducing"
                ) from None
            if ell > 0.0:
                yield found, float(r), float(th), t1, t2
                found += 1
                if found >= clusters:
                    break
        if progress is not None:
            now = time.perf_counter() - t0
            if now - last >= 60.0:
                last = now
                progress(seen, found, now)


# ------------------------------------------------------- the G3 probes

def _leg(r_from: float, r_to: float, dt: float, dpsi: float,
         tol: float, leg: str, returned) -> dict:
    """One `causal_relation` call, opened up.

    `t_min` and `err` are re-derived from the same `flight_time` call
    the predicate makes, at the same `dpsi`, so the record is what the
    predicate saw and not a paraphrase. The re-derivation is then
    checked against the value the predicate actually returned; a
    disagreement would mean this module has misread the predicate, and
    the walk counts every one.

    `dt < 0` is the predicate's one short circuit: it answers `False`
    without consulting the solver at all, so no margin exists there."""

    got = ("undecided" if returned is None
           else ("true" if returned else "false"))
    if dt < 0.0:
        return {"leg": leg, "dt": dt, "dpsi": dpsi,
                "short_circuited_negative_dt": True,
                "undecided": False, "returned": got,
                "rederivation_agrees": got == "false"}
    t_min, err = s1.flight_time(r_from, r_to, dpsi, s1.M, tol)
    gap = abs(dt - t_min)
    undecided = gap <= err
    expect = "undecided" if undecided else (
        "true" if dt > t_min else "false")
    return {
        "leg": leg, "dt": dt, "dpsi": dpsi, "t_min": t_min,
        "err": err, "abs_dt_minus_t_min": gap,
        "decision_margin": gap - err,
        "short_circuited_negative_dt": False,
        "undecided": undecided,
        "returned": got,
        "rederivation_agrees": got == expect,
    }


def probe_cluster(index: int, r: float, theta: float, t1: float,
                  t2: float, tol: float) -> dict:
    """The frozen G3 probes at one cluster, with the tri-state opened.

    Geometry, probe order and offsets are `run_g3`'s; what differs is
    that an undecided return is recorded instead of aborting."""

    off = o4.FROZEN["g3_stress_offset"]
    lo_edge, hi_edge = t1, sz.DT - t2
    p = np.array([0.0, sz.R_IN, 0.0, 0.0])
    q = np.array([sz.DT, sz.R_OUT, 0.0, 0.0])
    probes = []
    for kind, t_x, want in ((MIDPOINT, 0.5 * (lo_edge + hi_edge), True),
                            (OUTSIDE, hi_edge + off, False)):
        x = np.array([t_x, r, theta, 0.0])
        a = s1.causal_relation(p, x, s1.M, tol)
        b = s1.causal_relation(x, q, s1.M, tol)
        legs = [
            _leg(sz.R_IN, r, t_x, _dpsi(p, x), tol, LEG_PX, a),
            _leg(r, sz.R_OUT, sz.DT - t_x, _dpsi(x, q), tol,
                 LEG_XQ, b),
        ]
        undecided = a is None or b is None
        probes.append({
            "probe": kind, "t_x": t_x, "want": want,
            "undecided": undecided,
            "boolean_mismatch": (not undecided
                                 and bool(a and b) is not want),
            "legs": legs,
        })
    return {
        "cluster_index": index,
        "r": r, "theta": theta, "T1": t1, "T2": t2,
        "window": {"lo": lo_edge, "hi": hi_edge,
                   "L": hi_edge - lo_edge},
        "probes": probes,
    }


def _cause_of(probe: dict) -> str | None:
    """The one cause a probe carries, if any. Undecided wins: the
    frozen runner never reaches the boolean test when either leg is
    undecided, so the two are mutually exclusive by construction."""

    if probe["undecided"]:
        return (CAUSE_MIDPOINT if probe["probe"] == MIDPOINT
                else CAUSE_OUTSIDE)
    return CAUSE_MISMATCH if probe["boolean_mismatch"] else None


def cluster_causes(record: dict) -> list[str]:
    """The causes present at one cluster, in the frozen probe order,
    each listed once."""

    found: list[str] = []
    for probe in record["probes"]:
        cause = _cause_of(probe)
        if cause is not None and cause not in found:
            found.append(cause)
    return found


def site(record: dict, cause: str) -> dict:
    """A full failure site: cluster index, coordinates, probe, and --
    for an undecided -- the leg with every quantity the decision turned
    on. A boolean mismatch is a property of the probe rather than of
    one leg, so both legs are carried."""

    probe = next(p for p in record["probes"]
                 if _cause_of(p) == cause)
    where = {"probe": probe["probe"], "t_x": probe["t_x"],
             "want": probe["want"]}
    if cause == CAUSE_MISMATCH:
        where["legs"] = probe["legs"]
    else:
        where |= next(lg for lg in probe["legs"] if lg["undecided"])
    return {
        "cause": cause,
        "cluster_index": record["cluster_index"],
        "r": record["r"], "theta": record["theta"],
        "T1": record["T1"], "T2": record["T2"],
        "window": record["window"],
        "site": where,
    }


# ---------------------------------------------------------------- walk

def walk(clusters: int, tol: float, progress=None) -> dict:
    """Probe the leading `clusters` stress points and count causes.

    The walk continues past every undecided -- G3 could not."""

    seed = replay_scalar(_STREAM)
    per_cause = dict.fromkeys(CAUSES, 0)
    per_cause_legs = dict.fromkeys(CAUSES, 0)
    firsts: dict[str, dict] = {}
    first_undecided = None
    disagreements = clean = n = 0
    for index, r, th, t1, t2 in stress_clusters(
            seed, clusters, tol, progress):
        rec = probe_cluster(index, r, th, t1, t2, tol)
        n += 1
        for probe in rec["probes"]:
            disagreements += sum(1 for lg in probe["legs"]
                                 if not lg["rederivation_agrees"])
            cause = _cause_of(probe)
            if cause is None:
                continue
            per_cause_legs[cause] += (
                sum(1 for lg in probe["legs"] if lg["undecided"])
                if cause != CAUSE_MISMATCH else 1)
        causes = cluster_causes(rec)
        if not causes:
            clean += 1
            continue
        for cause in causes:
            per_cause[cause] += 1
            firsts.setdefault(cause, site(rec, cause))
        undecided = [c for c in causes if c != CAUSE_MISMATCH]
        if first_undecided is None and undecided:
            first_undecided = site(rec, undecided[0])
    return {
        "clusters_probed": n,
        "clusters_requested": clusters,
        "clusters_frozen": o4.FROZEN["g3_clusters"],
        "covers_frozen_stress_set": n == o4.FROZEN["g3_clusters"],
        "clusters_with_no_cause": clean,
        "first_undecided": first_undecided,
        "first_site_per_cause": [
            {"cause": c, **{k: v for k, v in firsts[c].items()
                            if k != "cause"}}
            for c in CAUSES if c in firsts],
        "counts": [{"cause": c, "clusters": per_cause[c],
                    "occurrences": per_cause_legs[c]} for c in CAUSES],
        "counts_are": ("diagnostic frequencies over ONE frozen stress "
                       "set reproduced from a retired stream; they are "
                       "not estimates of a rate and carry no interval. "
                       "`clusters` counts stress points carrying the "
                       "cause; `occurrences` counts undecided LEGS for "
                       "the two band causes and mismatching PROBES for "
                       "the third, since a mismatch is a property of "
                       "the probe rather than of one leg"),
        "rederivation_disagreements": disagreements,
        "rederivation_note": (
            "this module re-derives dt/t_min/err from the same "
            "flight_time call the predicate makes and checks the "
            "result against what causal_relation actually returned; a "
            "non-zero count would mean the readout is wrong, not the "
            "predicate"),
    }


# ------------------------------------------------------------ artifact

def assemble(clusters: int, progress=None) -> dict:
    surface = verify_replay_surface()
    tol = o4.FROZEN["tol"]
    t0 = time.perf_counter()
    found = walk(clusters, tol, progress)
    return {
        "run_kind": RUN_KIND,
        "replay_of": REPLAY_OF,
        "incident_record": _INCIDENT.relative_to(_REPO).as_posix(),
        "not_a_result": (
            "A replay is a reproduction, never an observation. Nothing "
            "here is a scientific quantity, no gate has a status, and "
            "no frozen rule changes as a consequence: the G3 redesign "
            "is a separate prereg re-opening."),
        "replayed_stream": {
            "name": _STREAM, "seed": _STREAM_SEED,
            "via": "probe_seed_ledger.replay_scalar",
            "state": "retired (observed, never re-entered)",
        },
        "streams_untouched": {
            "withdrawn_unspent": _WITHDRAWN_UNSPENT,
            "g2": ("o4_aborted_g2 is not replayed; G2 is not re-run at "
                   "all, since the G3 stress points come from G1"),
            "remote_reservation": ("refs/o4/reservation is retained; "
                                   "this module neither reads it as a "
                                   "claim nor writes it"),
        },
        "scope": {
            "reproduced": ("the G1 draw prefix, only as far as the "
                           "original G3 stress points required"),
            "not_computed": ["G1 statistics", "G2 leakage",
                             "any verdict"],
            "interleaving": (
                "each cluster is probed as soon as it is accepted; the "
                "four probe calls depend on that cluster's own "
                "(r, theta, T1, T2) alone, so the inputs to "
                "causal_relation are identical to the frozen order"),
            "does_not_stop_at_first_undecided": True,
        },
        "surface": surface,
        "frozen": {k: o4.FROZEN[k] for k in
                   ("tol", "g3_clusters", "g3_stress_offset")},
        "findings": found,
        "wall_s": time.perf_counter() - t0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "--clusters", type=int, default=o4.FROZEN["g3_clusters"],
        help="how many of the frozen stress clusters to probe "
             "(default: all of them)")
    parser.add_argument(
        "--out", metavar="PATH",
        help="publish the diagnostic artifact here, no-clobber; "
             "omitted, the run only prints")
    args = parser.parse_args()

    if not 1 <= args.clusters <= o4.FROZEN["g3_clusters"]:
        raise SystemExit(
            f"--clusters must lie in [1, {o4.FROZEN['g3_clusters']}]: "
            f"the stress set is frozen and this diagnostic may not "
            f"probe past it")

    def show(seen: int, found: int, secs: float) -> None:
        print(f"  prefix {seen:,} pts  clusters {found:,}/"
              f"{args.clusters:,}  t={secs:.0f}s", flush=True)

    art = assemble(args.clusters, progress=show)
    found = art["findings"]
    first = found["first_undecided"]
    if first is None:
        print(f"no undecided in the leading {found['clusters_probed']:,} "
              f"clusters")
    else:
        s = first["site"]
        print(f"first undecided: cluster {first['cluster_index']}  "
              f"{first['cause']}")
        print(f"  r={first['r']!r}  theta={first['theta']!r}  "
              f"L={first['window']['L']!r}")
        print(f"  probe={s['probe']}  leg={s['leg']}  "
              f"dt={s['dt']!r}  t_min={s['t_min']!r}")
        print(f"  err={s['err']!r}  |dt-t_min|="
              f"{s['abs_dt_minus_t_min']!r}  margin="
              f"{s['decision_margin']!r}")
    for row in found["counts"]:
        print(f"  {row['cause']}: {row['clusters']:,} clusters, "
              f"{row['occurrences']:,} occurrences")
    print(f"clean clusters: {found['clusters_with_no_cause']:,} / "
          f"{found['clusters_probed']:,}")
    if found["rederivation_disagreements"]:
        print(f"WARNING: {found['rederivation_disagreements']} legs "
              f"where the readout disagrees with causal_relation")

    if args.out:
        out = Path(args.out).resolve()
        o4._publish_write_once(
            out, json.dumps(art, ensure_ascii=False, indent=1) + "\n")
        print(f"artifact: {out}")


if __name__ == "__main__":
    main()
