"""O3: the certified diamond volume of the FROZEN anchors -- runner.

EXECUTION IS NOT YET APPROVED. The freeze-review ruling (PR #66
aftermath) approves PREPARING this freeze; the campaign run itself
needs a separate approval after the freeze PR review converges and
an exact-checkout preflight passes. Until then the frozen
configuration's volume (12, 18, 8.5)M stays UNOBSERVED -- running it
early and deciding anything afterwards would collapse the freeze
boundary.

The frozen choices and their epistemic grade:

- n_sub = 16 is a freeze CHOICE based on the current measurements
  (the neighbor price ladder p14_oracle_price.json and the
  mode-width diagnostic p14_oracle_mode_width.json), NOT a claim of
  optimality.
- max_calls = 600,000 and max_wall_s = 24 h are HEADROOM above the
  fit-window extrapolation range 187k-237k calls (~2.53x its top,
  ~14.8 h at the measured call rate), sized for the frozen diamond
  being larger than the neighbor's and for the extrapolation being a
  model. They are frozen: RAISING THEM DURING OR AFTER THE RUN IS
  FORBIDDEN. If they fire, the result is the certified interval plus
  `target-not-met` plus the exact termination reason, reported as
  is.
- target_ratio = 0.01 is the frozen target from the certification
  document; unchanged.

Freeze identity is content-addressed (the S4/S5 pattern): an 8-file
manifest of raw SHA-256 digests plus a dependency/environment lock
(python, gmpy2, MPFR/GMP -- the correct-rounding contract lives in
those libraries, so a different MPFR is a different instrument).
`verify_freeze` runs at entry AND exit; `--preflight` additionally
checks a clean git tree and the ABSENCE of the result artifact, and
observes nothing.

Run (after separate execution approval, from a clean exact
checkout):

    python experiments/oracle/o3_frozen_volume.py --preflight
    python experiments/oracle/o3_frozen_volume.py
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
_MANIFEST = _REPO / "docs" / "prereg" / "p14_o3_freeze_manifest.json"
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o3_volume.json"

#: The frozen configuration (freeze-review ruling). Caps are part of
#: the freeze: no auto-raise, ever.
FROZEN = {
    "r_in": 12.0, "r_out": 18.0, "dt": 8.5, "m": 1.0,
    "target_ratio": 0.01,
    "n_sub": 16, "k_micro": 4, "d_switch": 0.25,
    "init_rho": 12, "init_psi": 12,
    "max_depth": 18,
    "max_calls": 600_000,
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
        # both libraries, by name: the correct-rounding contract is
        # MPFR's; GMP is the arithmetic underneath it
        "mpfr": gmpy2.mpfr_version(),
        "gmp": gmpy2.mp_version(),
    }


def verify_digests(stage: str) -> dict:
    """Content-addressed freeze identity, platform-independent part:
    every pinned file matches its frozen digest (all LF-pinned, so
    the digests are checkout-independent). Refuses on any mismatch
    -- a clean-but-later tree is a drifted protocol surface."""

    manifest = json.loads(_MANIFEST.read_text(encoding="utf-8"))
    bad = [rel for rel, want in manifest["files"].items()
           if _sha256(_REPO / rel) != want]
    if bad:
        raise SystemExit(
            f"freeze verification ({stage}): digest mismatch on "
            f"{bad} -- refusing to run on a drifted protocol surface")
    return manifest


def verify_environment(stage: str, manifest: dict) -> None:
    """Execution-host part of the freeze identity: the running
    python/gmpy2/MPFR/GMP must equal the lock. This binds the RUN,
    not CI -- other platforms legitimately carry other MPFR builds,
    which is exactly why they may not execute the campaign."""

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
    the ABSENCE of any result -- observing nothing."""

    manifest = verify_freeze("preflight")
    state = _git_state()
    if state["dirty"]:
        raise SystemExit(
            "preflight: working tree is dirty -- the campaign runs "
            "from a clean exact checkout only")
    if _ARTIFACT.exists():
        raise SystemExit(
            f"preflight: {_ARTIFACT.name} already exists -- the "
            f"frozen volume is write-once; a rerun may not overwrite "
            f"or relabel the first observation")
    return {"manifest": manifest, "git": state}


def _serialize_result(res: dict) -> dict:
    """The certified interval's endpoints serialized OUTWARD
    (lo_float/hi_float, never plain float()): a nearest-rounding
    conversion can move the lower endpoint up or the upper endpoint
    down, after which the stored numbers no longer contain the true
    value and the artifact stops being a certified interval (review
    R1). This artifact is the frozen run's only certified result, so
    the binary64 pair must still enclose the MPFR interval."""

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
    """Atomic no-clobber publication (review R1): the payload is
    fully written and fsynced to a per-process temp file in the SAME
    directory, then hard-linked to the destination -- os.link fails
    atomically if the destination exists, so a concurrent second run
    can never overwrite the first observation, and a crash mid-write
    leaves only the temp file, never a truncated destination that a
    later preflight would mistake for the write-once result."""

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
                f"publish: {path.name} already exists -- the frozen "
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
              f"clean tree at {checks['git']['rev'][:9]}, no result "
              "artifact. Nothing was observed.")
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

    art = {
        "stage": ("O3: certified diamond volume of the frozen "
                  "anchors (12, 18, 8.5)M"),
        "frozen_config": FROZEN,
        "environment": _environment(),
        "host": {"machine": platform.machine(),
                 "system": platform.system()},
        "code": {"start": start, "end": end},
        "result": _serialize_result(res),
        "curve": res["curve"],
        "total_wall_s": time.perf_counter() - t0,
    }
    payload = json.dumps(art, ensure_ascii=False, indent=1) + "\n"
    _publish_write_once(_ARTIFACT, payload)
    r = art["result"]
    print(f"result: V=[{r['v_lo']:.6f}, {r['v_hi']:.6f}] "
          f"ratio={res['ratio']:.6f} {res['status']} "
          f"({res['termination_reason']})")
    print(f"artifact: {_ARTIFACT}")


if __name__ == "__main__":
    main()
