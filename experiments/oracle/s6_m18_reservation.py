"""The S6 M=1.8 ambiguity pilot's seed reservation, in its own ref
namespace.

The reasoning is O4b's unchanged (its reviews R2/R3/R12/R23/R24): the
remote ref is the one authority that serialises stream opening across
every checkout; the claim object is unique per ATTEMPT by a 16-byte
nonce; a push whose reply is lost has still happened at the server, so
every uncertain outcome is treated as SEEDS POSSIBLY SPENT.

What is new here is only the namespace and the retained set:

  * this pilot claims `refs/s6m18pilot/reservation`;
  * ALL SIX prior reservations are retained permanently and verified
    before any claim -- refs/o4, refs/o4b, refs/o5pilot, refs/o5,
    refs/s6m14pilot and refs/s6m14 (every executed or in-flight
    claim on drawn streams). A cleared ref
    would make retired streams look free from a fresh clone, which is
    the one failure this program has already decided must never
    happen.

`ClaimUncertain`, `authority` and `_ls_remote` are IMPORTED from
`o4b_reservation`; `_make_commit`/`claim`/`verify_still_held` follow
`o5_pilot_reservation`'s corrected copies (the PR #85 pre-push
contract: every failure before the push is a clean SystemExit
refusal, never a raw subprocess exception).
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

_HERE = Path(__file__).resolve().parent
import sys  # noqa: E402

sys.path.insert(0, str(_HERE))

import o4b_reservation as base  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]

REMOTE = base.REMOTE
CANONICAL_AUTHORITY = base.CANONICAL_AUTHORITY
ClaimUncertain = base.ClaimUncertain

#: This campaign's ref.
REF = "refs/s6m18pilot/reservation"

#: Retained permanently, never written here, verified before any claim.
RETAINED = (
    ("refs/o4/reservation", "c4da1626463e6a6505374813cf3f56d6b429c209"),
    ("refs/o4b/reservation", "46acee340bc247511546964b2925953721d5bb59"),
    ("refs/o5pilot/reservation",
     "f5223636bc7c60603b18ecae3f474dd3dc37146f"),
    ("refs/o5/reservation",
     "b5d9d6131644e7c8c07477c0b676696aaf2d1ce1"),
    ("refs/s6m14pilot/reservation",
     "3fc1b28bdf8014245078f2ac3fd08faa17ea1f5b"),
    ("refs/s6m14/reservation",
     "a6c8fb74cbe7910d05f228368bef20843f3d8ff9"),
)

_IDENTITY = "s6m18pilot-reservation"


def held() -> str | None:
    """The object at THIS pilot's ref, or None."""

    return base._ls_remote(REF)


def verify_retained() -> None:
    """Every prior reservation is still exactly where its record left
    it -- checked, not assumed."""

    for ref, want in RETAINED:
        got = base._ls_remote(ref)
        if got is None:
            raise SystemExit(
                f"{ref} is GONE -- it is retained permanently and "
                f"holds a claim on retired streams; refusing to run "
                f"with it cleared")
        if got != want:
            raise SystemExit(
                f"{ref} holds {got}, not {want} -- it has been "
                f"rewritten, and the retired streams' claim is no "
                f"longer the one its record froze")


def _make_commit(message: str) -> str:
    """A commit over the empty tree under this campaign's identity."""

    env = dict(os.environ)
    for role in ("AUTHOR", "COMMITTER"):
        env[f"GIT_{role}_NAME"] = _IDENTITY
        env[f"GIT_{role}_EMAIL"] = f"{_IDENTITY}@invalid"
    # every failure here is PRE-PUSH by construction, so it must
    # surface as a clean SystemExit refusal, never as a raw
    # CalledProcessError/OSError that main()'s incident boundary
    # would file as a spent-nothing global incident (the PR #85
    # contract)
    try:
        hashed = subprocess.run(
            ["git", "hash-object", "-t", "tree", "--stdin"],
            cwd=_REPO, input="", capture_output=True, text=True)
    except OSError as exc:
        raise SystemExit(
            f"reservation: `git hash-object` could not run ({exc}) "
            f"-- nothing was pushed") from None
    if hashed.returncode != 0:
        raise SystemExit(
            f"reservation: `git hash-object` failed: "
            f"{hashed.stderr.strip()}")
    empty_tree = hashed.stdout.strip()
    try:
        made = subprocess.run(
            ["git", "commit-tree", empty_tree, "-m", message],
            cwd=_REPO, input="", capture_output=True, text=True,
            env=env)
    except OSError as exc:
        raise SystemExit(
            f"reservation: `git commit-tree` could not run ({exc}) "
            f"-- nothing was pushed") from None
    if made.returncode != 0:
        raise SystemExit(
            f"reservation: `git commit-tree` failed: "
            f"{made.stderr.strip()}")
    return made.stdout.strip()


def probe_namespace() -> None:
    """Prove at PREFLIGHT that `refs/s6m18pilot/` is writable: a NEW
    namespace, so more than a formality. The probe uses a different
    ref, is deleted immediately, and observes nothing."""

    base.authority()
    probe_ref = REF + "-preflight-probe"
    obj = _make_commit("s6m18pilot preflight probe")
    push = subprocess.run(
        ["git", "push", "--force", REMOTE, f"{obj}:{probe_ref}"],
        cwd=_REPO, capture_output=True, text=True)
    if push.returncode != 0:
        raise SystemExit(
            f"preflight: cannot write {REMOTE}/{probe_ref}: "
            f"{push.stderr.strip()} -- the pilot could not claim "
            f"its stream, so it must not start")
    drop = subprocess.run(["git", "push", REMOTE, f":{probe_ref}"],
                          cwd=_REPO, capture_output=True, text=True)
    if drop.returncode != 0:
        print(f"warning: left {probe_ref} behind on {REMOTE}; delete "
              f"it manually ({drop.stderr.strip()})")


def claim(payload: dict) -> str:
    """Claim the pilot stream globally, before the first draw, and
    RETURN the attempt object (the O4b abort's one-line lesson)."""

    verify_retained()
    if (already := held()) is not None:
        raise SystemExit(
            f"reservation: {REF} is already held by {already} -- a "
            f"pilot has already opened this stream from some "
            f"checkout; the seed is spent whether or not that run "
            f"published")
    attempt = {**payload, "attempt_nonce": os.urandom(16).hex()}
    obj = _make_commit(
        json.dumps(attempt, ensure_ascii=False, sort_keys=True))
    # FROM HERE ON THE OUTCOME IS UNCERTAIN, NOT MERELY UNSUCCESSFUL.
    try:
        push = subprocess.run(
            ["git", "push", f"--force-with-lease={REF}:", REMOTE,
             f"{obj}:{REF}"],
            cwd=_REPO, capture_output=True, text=True)
    except BaseException as exc:
        raise ClaimUncertain(
            f"the push to {REF} did not return",
            pushed=None, obj=obj,
            detail=f"{type(exc).__name__}: {exc}") from exc
    if push.returncode != 0:
        raise ClaimUncertain(
            f"could not claim {REF} -- another checkout may have won "
            f"the race, or the reply was lost after the ref was "
            f"created", pushed=False, obj=obj,
            detail=push.stderr.strip())
    try:
        confirmed = held()
    except BaseException as exc:
        raise ClaimUncertain(
            f"{REF} was pushed but could not be re-read",
            pushed=True, obj=obj, detail=str(exc)) from exc
    if confirmed != obj:
        raise ClaimUncertain(
            f"{REF} holds {confirmed}, not this attempt's object",
            pushed=True, obj=obj, detail=f"held={confirmed}")
    return obj


def verify_still_held(obj: str) -> str:
    """The ref still holds THIS attempt's object, re-read immediately
    before publication. A read failure here is a reason not to
    publish, never uncertain-and-continue."""

    try:
        got = base._ls_remote(REF)
    except BaseException as exc:
        raise ClaimUncertain(
            f"{REF} could not be re-read before publication",
            pushed=True, obj=obj,
            detail=f"{type(exc).__name__}: {exc}") from exc
    if got is None:
        raise ClaimUncertain(
            f"{REF} is GONE -- the claim this run drew under no "
            f"longer exists", pushed=True, obj=obj,
            detail="ref absent at exit")
    if got != obj:
        raise ClaimUncertain(
            f"{REF} holds {got}, not this run's {obj}",
            pushed=True, obj=obj, detail=f"held={got}")
    return got
