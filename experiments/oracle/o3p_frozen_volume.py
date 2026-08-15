"""O3': the re-certified diamond volume at target ratio 0.005 -- runner.

EXECUTION IS NOT YET APPROVED. The PI ruling on the O5 sizing audit
adopts the target r = 0.005 and approves PREPARING this freeze; the
campaign run itself needs a separate approval after the freeze PR
review converges and an exact-checkout preflight passes. Until then
the frozen configuration's volume stays observed only at the O3
precision -- running O3' early and deciding anything afterwards would
collapse the freeze boundary.

WHY O3' EXISTS. The O5 Poisson-count stage is infeasible against the
current O3 interval: its full width W = 1.135282 sets an equivalence
floor W/V_ref = 2.0%, and the exact sizing (three-engine reproduced)
puts even tau = 3.0% beyond the campaign call caps. Halving the
oracle's relative half-width to r = 0.005 lowers the floor to 1.0%
and makes tau_O5 = 2.5% feasible at ~5.4e7 expected predicate calls.
O3' is deterministic certified integration: NO seed, NO sprinkling,
no stochastic draw anywhere.

The frozen choices and their epistemic grade:

- target_ratio = 0.005 is the ADOPTED target (PI ruling): the
  back-computed minimal oracle width that strengthens the O5 allowed
  sentence within campaign caps.
- n_sub / k_micro / d_switch / init grids are INHERITED from the O3
  freeze unchanged -- a freeze CHOICE, not a claim of optimality.
- max_calls = 900,000 and max_wall_s = 24 h are HEADROOM above the
  projection RANGE ~467k-750k calls (~12.6-20.2 h), obtained by
  fitting the neighbor price ladder AND the executed O3 run's own
  1,078-sample curve separately over several fit windows
  (`o3p_projection.py`; +20% above the worst window). Projections,
  not certifications. The caps are frozen: RAISING THEM DURING OR
  AFTER THE RUN IS FORBIDDEN. If they fire, the result is the
  certified interval reached at that moment plus `target-not-met`
  plus the exact termination reason and full provenance, published
  as is.
- max_depth = 18 is the O3 value unchanged (O3 reached depth 12).

THE O3 ARTIFACTS ARE IMMUTABLE. This runner never writes them; it
READS the committed `p14_o3_volume.json` (digest-pinned in this
freeze's manifest) as the intersection base. The artifact reports the
O3' STANDALONE interval and, alongside it, the intersection with the
O3 interval. Both certifications enclose the same true volume, so an
EMPTY intersection is a CERTIFICATION INCONSISTENCY: it indicts at
least one of the two certifications, the artifact says so
machine-readably, and no downstream stage (O5 sizing, V_ref') may
consume either interval until the inconsistency is resolved. The
standalone O3' interval -- not the intersection -- is what a
consistent result hands to O5, with V_ref' = its midpoint as the
recommended reference (final adoption at the O5 freeze).

Freeze identity is content-addressed (the O3/S4/S5 pattern): a
manifest of raw SHA-256 digests plus the environment lock (python,
gmpy2, MPFR/GMP -- the correct-rounding contract lives in those
libraries, so a different MPFR is a different instrument).
`verify_freeze` runs at entry AND exit; `--preflight` additionally
checks a clean git tree and the ABSENCE of the O3' result artifact,
and observes nothing.

Run (after separate execution approval, from a clean exact
checkout):

    python experiments/oracle/o3p_frozen_volume.py --preflight
    python experiments/oracle/o3p_frozen_volume.py
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
from volume_oracle import OracleConfig, assemble  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_MANIFEST = _REPO / "docs" / "prereg" / "p14_o3p_freeze_manifest.json"
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o3p_volume.json"
#: The immutable O3 result this freeze reads as intersection base.
_O3_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o3_volume.json"

#: The frozen configuration (PI ruling on the O5 sizing audit). Caps
#: are part of the freeze: no auto-raise, ever.
FROZEN = {
    "r_in": 12.0, "r_out": 18.0, "dt": 8.5, "m": 1.0,
    "target_ratio": 0.005,
    "n_sub": 16, "k_micro": 4, "d_switch": 0.25,
    "init_rho": 12, "init_psi": 12,
    "max_depth": 18,
    "max_calls": 900_000,
    "max_wall_s": 86_400.0,
}


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
    """Content-addressed freeze identity: every pinned file matches
    its frozen digest (all LF-pinned, so the digests are
    checkout-independent). The pin list includes `p14_o3_volume.json`,
    so the intersection base cannot drift either."""

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


def preflight() -> dict:
    """Everything the execution approval requires re-checked at the
    last moment: frozen digests, environment lock, clean tree, and
    the ABSENCE of any O3' result -- observing nothing. The O3
    artifact must EXIST (it is the intersection base), and stays
    untouched."""

    manifest = verify_freeze("preflight")
    state = _git_state()
    if state["dirty"]:
        raise SystemExit(
            "preflight: working tree is dirty -- the campaign runs "
            "from a clean exact checkout only")
    if _ARTIFACT.exists():
        raise SystemExit(
            f"preflight: {_ARTIFACT.name} already exists -- the "
            f"O3' volume is write-once; a rerun may not overwrite "
            f"or relabel the first observation")
    if not _O3_ARTIFACT.exists():
        raise SystemExit(
            f"preflight: {_O3_ARTIFACT.name} is missing -- the O3 "
            f"result is this freeze's immutable intersection base")
    return {"manifest": manifest, "git": state}


def o3_interval() -> dict:
    """The immutable O3 base, read only. Its digest is pinned in this
    freeze's manifest, so what is read here is what was frozen."""

    o3 = json.loads(_O3_ARTIFACT.read_text(encoding="utf-8"))
    r = o3["result"]
    return {"v_lo": r["v_lo"], "v_hi": r["v_hi"],
            "ratio": r["ratio"], "status": r["status"],
            "artifact": _O3_ARTIFACT.name,
            "artifact_sha256": _sha256(_O3_ARTIFACT)}


def intersect_with_o3(v_lo: float, v_hi: float) -> dict:
    """The O3' standalone interval against the immutable O3 interval.

    Both certifications claim to enclose the same true volume, so a
    consistent pair overlaps. An empty intersection does NOT pick a
    winner: it is a CERTIFICATION INCONSISTENCY that indicts at least
    one of the two, and the record says so machine-readably. On
    consistency, the STANDALONE O3' interval is what downstream
    consumes; the intersection is the consistency check, not the
    estimator."""

    base = o3_interval()
    lo = max(v_lo, base["v_lo"])
    hi = min(v_hi, base["v_hi"])
    consistent = lo <= hi
    out = {
        "o3_base": base,
        "consistent": consistent,
        "certification_inconsistency": not consistent,
    }
    if consistent:
        out["v_lo"] = lo
        out["v_hi"] = hi
        out["note"] = (
            "consistency check only: downstream (O5, V_ref') consumes "
            "the STANDALONE O3' interval, not this intersection")
    else:
        out["note"] = (
            "EMPTY intersection: the O3 and O3' certifications cannot "
            "both enclose the true volume -- at least one is wrong. "
            "No downstream stage may consume either interval until "
            "this is resolved; this record publishes the fact, it "
            "does not adjudicate it")
    return out


def _serialize_result(res: dict) -> dict:
    """Outward binary64 endpoints (lo_float/hi_float, never plain
    float()), exactly as O3: a nearest-rounding conversion could move
    an endpoint inward and the stored pair would stop enclosing the
    true value. `status`/`termination_reason` are published AS IS --
    a `target-not-met` result is still the certified interval reached
    at that moment, with full provenance."""

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
    """Atomic no-clobber publication, exactly the O3 mechanism: fsync
    to a same-directory temp, then os.link, which fails atomically if
    the destination exists."""

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
                f"publish: {path.name} already exists -- the O3' "
                f"volume is write-once and the first observation "
                f"stands; this run's payload is NOT published"
            ) from None
    finally:
        tmp.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "--preflight", action="store_true",
        help="run every pre-execution check and observe nothing")
    args = parser.parse_args()

    checks = preflight()
    if args.preflight:
        print("preflight PASS: freeze digests, environment lock, "
              f"clean tree at {checks['git']['rev'][:9]}, no O3' "
              "result artifact, O3 base present. Nothing was "
              "observed.")
        return

    start = _git_state()
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
    if end["rev"] != start["rev"] or end["dirty"]:
        raise SystemExit(
            "the tree changed underneath the campaign -- refusing "
            "to write a result with broken lineage")

    result = _serialize_result(res)
    art = {
        "stage": ("O3': re-certified diamond volume of the frozen "
                  "anchors (12, 18, 8.5)M at target ratio 0.005"),
        "frozen_config": FROZEN,
        "environment": _environment(),
        "host": {"machine": platform.machine(),
                 "system": platform.system()},
        "code": {"start": start, "end": end},
        "result": result,
        "intersection_with_o3": intersect_with_o3(
            result["v_lo"], result["v_hi"]),
        "v_ref_prime_recommendation": {
            "value": 0.5 * (result["v_lo"] + result["v_hi"]),
            "definition": "midpoint of the STANDALONE O3' interval",
            "status": ("recommendation only -- adopted, and all O5 "
                       "sizing re-derived from the actual endpoints, "
                       "at the O5 freeze"),
        },
        "curve": res["curve"],
        "total_wall_s": time.perf_counter() - t0,
    }
    payload = json.dumps(art, ensure_ascii=False, indent=1) + "\n"
    _publish_write_once(_ARTIFACT, payload)
    r = art["result"]
    x = art["intersection_with_o3"]
    print(f"result: V=[{r['v_lo']:.6f}, {r['v_hi']:.6f}] "
          f"ratio={r['ratio']:.6f} {r['status']} "
          f"({r['termination_reason']})")
    verdict = ("consistent" if x["consistent"]
               else "EMPTY -- CERTIFICATION INCONSISTENCY")
    print(f"intersection with O3: {verdict}")
    print(f"artifact: {_ARTIFACT}")


if __name__ == "__main__":
    main()
