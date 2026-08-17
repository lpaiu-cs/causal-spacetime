"""S6 rung M = 1.4: certified diamond volume at target ratio 0.005
-- the frozen runner.

EXECUTION IS NOT YET APPROVED. The S6 ruling approves implementing
this freeze; the run itself starts only after the M = 1.4 integrated
execution approval, from a clean exact checkout of the freeze branch
head (the merge commit's second parent), once. Until then this rung's
volume is known only as the DESIGN estimate (coarse survey,
non-binding, never a gate).

WHAT THIS RUNG IS. The mass ladder's first new rung: the certified
shell [10, 20], the anchors (12, 18) and the cap stay FIXED in
absolute coordinates; the mass is 1.4, so the compactness indicator
is MU = 2M/r_c = 0.18666666666666665 against the central rung's
0.13333... . The time window is the ladder convention, NOT a free
choice: dt = 8.5 * T_min(1.4) / T_min(1) = 9.070565190742672, pinned
in `s6_rungs.RUNG_CONSTANTS` and re-derived through the one exact
path at import. The import gate below also re-runs the rung's FULL
certified lemma table (exterior, K < 0, w-monotone, Q > 0, L2a, all
four L4 margins, the L5 winding cross-check) and refuses to load if
any row fails -- nothing is inherited silently from M = 1.

NO INTERSECTION BASE. Unlike O3', no prior certified interval exists
at this mass: the artifact publishes the STANDALONE certified
interval only. The design survey's coarse estimate is not a
certification and never becomes a gate; consistency across rungs is
a statement the S6 integration stage may make LATER, per verdict
separation (no joint primary).

The frozen choices and their epistemic grade:

- target_ratio = 0.005: the ladder-uniform target (S6 ruling), the
  same floor arithmetic that made tau = 2.5% feasible at M = 1.
- n_sub / k_micro / d_switch / init grids / max_depth: INHERITED
  from the O3/O3' freeze unchanged -- a freeze CHOICE, not a claim
  of optimality.
- max_calls = 900,000 and max_wall_s = 24 h: the O3' headroom kept.
  The O3' run cost 565,478 calls / 14.6 h at M = 1; the design
  survey converged FASTER per call at higher mass (smaller box), so
  the projection band ~467k-750k calls carries over as a projection,
  not a certification. The caps are frozen: RAISING THEM DURING OR
  AFTER THE RUN IS FORBIDDEN. If they fire, the result is the
  certified interval reached at that moment plus `target-not-met`
  plus the exact termination reason, published as is.

V_ref for this rung is REPORTED as the standalone midpoint,
recommendation only: adoption and every count-stage sizing (A,
acceptance, pilot n, powers -- all provisional in the design survey)
are re-derived at 96-bit from THIS artifact's actual endpoints, at
the rung's count freeze.

Deterministic certified integration: NO seed, NO sprinkling, no
stochastic draw anywhere.

Run (after the M = 1.4 integrated execution approval, from a clean
exact checkout of the approved SHA):

    python experiments/oracle/s6_m14_frozen_volume.py --preflight \
        --freeze-rev <approved 40-hex SHA>
    python experiments/oracle/s6_m14_frozen_volume.py \
        --freeze-rev <approved 40-hex SHA>

`--freeze-rev` is REQUIRED and must be the full 40-hex SHA: the
manifest cannot certify itself, so both `--preflight` and the real
execution refuse a missing, short, malformed, or
merely-digest-matching-but-different rev before anything runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gmpy2  # noqa: E402
import s6_rungs as s6  # noqa: E402
from volume_oracle import OracleConfig, assemble  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_MANIFEST = (_REPO / "docs" / "prereg"
             / "p14_s6_m14_freeze_manifest.json")
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s6_m14_volume.json"

_M = 1.4

#: The frozen configuration (S6 ruling). dt comes from the ladder
#: convention through the one exact path; caps are part of the
#: freeze: no auto-raise, ever.
FROZEN = {
    "r_in": 12.0, "r_out": 18.0,
    "dt": 9.070565190742672, "m": _M,
    "mu": 0.18666666666666665,
    "target_ratio": 0.005,
    "n_sub": 16, "k_micro": 4, "d_switch": 0.25,
    "init_rho": 12, "init_psi": 12,
    "max_depth": 18,
    "max_calls": 900_000,
    "max_wall_s": 86_400.0,
}


def _import_gate() -> dict:
    """The rung identity, refused at import if anything drifted:
    dt/mu must equal the ladder constants from the one exact path,
    and the rung's WHOLE certified lemma table must pass. The frozen
    table snapshot is embedded in the artifact at run time."""

    want = s6.RUNG_CONSTANTS[_M]
    if (FROZEN["dt"] != want["dt"] or FROZEN["mu"] != want["mu"]
            or FROZEN["m"] != _M):
        raise SystemExit(
            f"s6_m14: frozen (dt, mu) = ({FROZEN['dt']!r}, "
            f"{FROZEN['mu']!r}) is not the ladder constant "
            f"({want['dt']!r}, {want['mu']!r}) -- the one exact path "
            f"is the only source")
    table = s6.lemma_table(_M)
    if not table["all_pass"]:
        bad = [k for k, r in table["rows"].items()
               if not r["certified"]]
        raise SystemExit(
            f"s6_m14: certified lemma table FAILED on {bad} -- this "
            f"rung does not exist until every row certifies")
    return table


LEMMA_TABLE = _import_gate()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _environment() -> dict:
    """The instrument identity: the correct-rounding contract lives
    in MPFR, so its version is part of the freeze."""

    return {
        "python": ".".join(platform.python_version_tuple()[:2]),
        "gmpy2": gmpy2.version(),
        "mpfr": gmpy2.mpfr_version(),
        "gmp": gmpy2.mp_version(),
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
    env = _environment()
    locked = manifest["environment"]
    drift = {k: (locked[k], env.get(k)) for k in locked
             if env.get(k) != locked[k]}
    if drift:
        raise SystemExit(
            f"freeze verification ({stage}): environment drift "
            f"{drift} -- MPFR/gmpy2/python are part of the frozen "
            f"instrument")


def verify_freeze(stage: str) -> dict:
    manifest = verify_digests(stage)
    verify_environment(stage, manifest)
    return manifest


def _git_state() -> dict:
    rev = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=_REPO, check=True,
        capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], cwd=_REPO, check=True,
        capture_output=True, text=True).stdout.strip()
    return {"rev": rev, "dirty": bool(dirty)}


def require_full_sha(freeze_rev: str) -> str:
    """The approved execution SHA, exactly 40 hex characters -- the
    manifest cannot certify itself (the O4b/O3' rule)."""

    want = (freeze_rev or "").strip().lower()
    if len(want) != 40 or any(c not in "0123456789abcdef"
                              for c in want):
        raise SystemExit(
            f"--freeze-rev must be the full 40-hex SHA named in the "
            f"execution approval, got {freeze_rev!r} -- a short or "
            f"malformed rev is refused before anything runs")
    return want


def preflight(expected_rev: str) -> dict:
    """Everything the execution approval requires re-checked at the
    last moment: frozen digests, environment lock, clean tree AT THE
    EXACT APPROVED SHA, the certified lemma table, and the ABSENCE of
    any result for this rung -- observing nothing."""

    want = require_full_sha(expected_rev)
    manifest = verify_freeze("preflight")
    state = _git_state()
    if state["dirty"]:
        raise SystemExit(
            "preflight: working tree is dirty -- the campaign runs "
            "from a clean exact checkout only")
    if state["rev"] != want:
        raise SystemExit(
            f"preflight: HEAD is {state['rev']}, not the approved "
            f"freeze SHA {want} -- a clean tree with matching digests "
            f"is NOT enough; the campaign runs from the exact "
            f"approved commit only")
    if _ARTIFACT.exists():
        raise SystemExit(
            f"preflight: {_ARTIFACT.name} already exists -- this "
            f"rung's volume is write-once; a rerun may not overwrite "
            f"or relabel the first observation")
    if not s6.lemma_table(_M)["all_pass"]:
        raise SystemExit(
            "preflight: the certified lemma table no longer passes")
    return {"manifest": manifest, "git": state}


def v_ref_rung(v_lo: float, v_hi: float) -> dict:
    return {
        "value": 0.5 * (v_lo + v_hi),
        "definition": "midpoint of the STANDALONE M = 1.4 interval",
        "status": ("recommendation only -- adopted, and every "
                   "count-stage sizing (A, acceptance, pilot n, "
                   "powers) re-derived at 96-bit from the actual "
                   "endpoints, at this rung's count freeze"),
    }


def _serialize_result(res: dict) -> dict:
    """Outward binary64 endpoints, exactly as O3/O3'; status and
    termination reason published AS IS."""

    v = res["v"]
    return {
        "v_lo": v.lo_float(), "v_hi": v.hi_float(),
        "ratio": res["ratio"],
        "status": res["status"],
        "termination_reason": res["termination_reason"],
        "calls": res["calls"], "cells": res["cells"],
        "wall_s": res["wall_s"],
        "max_depth_reached": res["max_depth_reached"],
        "cells_at_max_depth": res["cells_at_max_depth"],
        "uncosted_cells": res["uncosted_cells"],
        "cell_counts_by_mode": res["modes"],
        "raw_width_by_mode": res["raw_width_by_mode"],
        "raw_total_before_intersection":
            res["raw_total_before_intersection"],
        "certified_total_after_intersection":
            res["certified_total_after_intersection"],
        "intersection_active": res["intersection_active"],
    }


def _publish_write_once(path: Path, payload: str) -> None:
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
                f"publish: {path.name} already exists -- this rung's "
                f"volume is write-once and the first observation "
                f"stands; this run's payload is NOT published"
            ) from None
    finally:
        tmp.unlink(missing_ok=True)


# ------------------------------------------------ the manifest

MANIFEST_NOTE = (
    "raw sha256; all paths .gitattributes eol=lf pinned. The rung "
    "inherits the certified integrator and flight-time contract "
    "unchanged; its identity (dt, mu, lemma table) comes from "
    "s6_rungs' one exact path and is re-derived at import. This "
    "manifest cannot certify itself: the run also demands the exact "
    "freeze branch head via --freeze-rev."
)

PROTOCOL_SURFACE = (
    "experiments/oracle/s6_m14_frozen_volume.py",
    "experiments/oracle/s6_rungs.py",
    "experiments/oracle/volume_oracle.py",
    "experiments/oracle/certified_flight_time.py",
    "experiments/oracle/certified_interval.py",
    "docs/theory/schwarzschild_volume_oracle_certification.md",
    "docs/prereg/p14_s6_m14_volume.md",
    "pyproject.toml",
)


def build_manifest() -> dict:
    return {
        "stage": ("S6 rung M = 1.4 oracle freeze manifest "
                  "(content-addressed protocol surface)"),
        "note": MANIFEST_NOTE,
        "files": {rel: _sha256(_REPO / rel)
                  for rel in PROTOCOL_SURFACE},
        "environment": _environment(),
    }


def write_manifest() -> Path:
    _MANIFEST.write_text(
        json.dumps(build_manifest(), ensure_ascii=False, indent=1)
        + "\n", encoding="utf-8", newline="\n")
    return _MANIFEST


def main() -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "--preflight", action="store_true",
        help="run every pre-execution check and observe nothing")
    parser.add_argument(
        "--freeze-rev", default="",
        help="REQUIRED for preflight and execution: the full 40-hex "
             "SHA named in the execution approval")
    parser.add_argument(
        "--write-manifest", action="store_true",
        help="regenerate the freeze manifest; a freeze action, never "
             "part of a run")
    args = parser.parse_args()

    if args.write_manifest:
        print(f"wrote {write_manifest()}")
        return

    checks = preflight(args.freeze_rev)
    if args.preflight:
        print("preflight PASS: freeze digests, environment lock, "
              f"clean tree at the approved "
              f"{checks['git']['rev'][:9]}, certified lemma table, "
              "no M=1.4 result artifact. Nothing was observed.")
        return

    # THE BASELINE IS THE PREFLIGHT'S VERIFIED STATE, NOT A RE-READ
    # (the PR #81 R3 rule)
    start = checks["git"]
    approved = start["rev"]
    cfg = OracleConfig(
        FROZEN["r_in"], FROZEN["r_out"], FROZEN["dt"], FROZEN["m"],
        target_ratio=FROZEN["target_ratio"], n_sub=FROZEN["n_sub"],
        k_micro=FROZEN["k_micro"], d_switch=FROZEN["d_switch"],
        max_calls=FROZEN["max_calls"],
        max_wall_s=FROZEN["max_wall_s"],
        max_depth=FROZEN["max_depth"],
        init_rho=FROZEN["init_rho"], init_psi=FROZEN["init_psi"])
    last = [0.0]

    def show(s: dict) -> None:
        if s["wall_s"] - last[0] < 300.0:
            return
        last[0] = s["wall_s"]
        print(f"  calls={s['calls']:7d} ratio={s['ratio']:.6f} "
              f"V=[{s['v_lo']:.5f}, {s['v_hi']:.5f}] "
              f"t={s['wall_s']:.0f}s", flush=True)

    t0 = time.perf_counter()
    res = assemble(cfg, progress=show)
    verify_freeze("exit")
    end = _git_state()
    # the exit rev is compared to the APPROVED SHA directly
    if end["rev"] != approved or end["dirty"]:
        raise SystemExit(
            f"exit lineage check: HEAD is {end['rev']} "
            f"(dirty={end['dirty']}), not the approved freeze SHA "
            f"{approved} -- the tree changed underneath the campaign; "
            f"refusing to write a result with broken lineage")

    result = _serialize_result(res)
    art = {
        "stage": ("S6 rung M = 1.4: certified diamond volume of the "
                  "fixed anchors (12, 18) at dt = 9.070565190742672, "
                  "target ratio 0.005"),
        "rung": {"m": _M, "mu": FROZEN["mu"],
                 "scale": s6.RUNG_CONSTANTS[_M]["scale"],
                 "ladder": list(s6.LADDER)},
        "lemma_table": {k: {kk: r[kk] for kk in
                            ("certified", "margin")}
                        for k, r in LEMMA_TABLE["rows"].items()},
        "frozen_config": FROZEN,
        "environment": _environment(),
        "host": {"machine": platform.machine(),
                 "system": platform.system()},
        "code": {"start": start, "end": end},
        "result": result,
        "v_ref_rung_recommendation": v_ref_rung(
            result["v_lo"], result["v_hi"]),
        "curve": res["curve"],
        "total_wall_s": time.perf_counter() - t0,
    }
    payload = json.dumps(art, ensure_ascii=False, indent=1) + "\n"
    _publish_write_once(_ARTIFACT, payload)
    r = art["result"]
    print(f"result: V=[{r['v_lo']:.6f}, {r['v_hi']:.6f}] "
          f"ratio={r['ratio']:.6f} {r['status']} "
          f"({r['termination_reason']})")
    print(f"artifact: {_ARTIFACT}")


if __name__ == "__main__":
    main()
