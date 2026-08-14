"""The O4b seed reservation, in its own ref namespace.

A file in the working tree is local state. `os.link` serialises only
processes sharing that directory, so two clones of the same approved
SHA would each create their own reservation and observe the same
streams, and a fail-closed abort followed by a discarded worktree
would leave the ledger reading `fresh` on seeds already drawn from.
The remote ref is one object shared by every checkout, and its update
is atomic at the server. That reasoning is O4's (its reviews R2 and
R3) and it is unchanged; what changes here is the namespace and the
one rule below.

`refs/o4/reservation` IS RETAINED AND IS NOT THIS RUN'S REF. It holds
`c4da162`, the O4 campaign that aborted, and it stays exactly as it
is: those streams were drawn from and are spent whether or not a
verdict came of them. O4b claims `refs/o4b/reservation` instead, and
`verify_o4_ref_retained()` checks the old one is still where the
incident says it is -- a run that had quietly cleared it would be
re-entering retired streams, which no digest check would notice.

WHY THE AUTHORITY IS THE CANONICAL REPOSITORY AND NOT `origin`.
`origin` is a local alias. A fork or mirror holding the same approved
SHA would find ITS OWN reservation ref empty and claim the same
streams on a different server, and the commit SHA would be identical,
so no rev check can see it.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]

#: This campaign's ref.
REMOTE = "origin"
REF = "refs/o4b/reservation"

#: O4's, retained permanently and never written by this module.
RETAINED_REF = "refs/o4/reservation"
RETAINED_OBJECT = "c4da1626463e6a6505374813cf3f56d6b429c209"

CANONICAL_AUTHORITY = "github.com/lpaiu-cs/causal-spacetime"

#: The identity the reservation commit is made under. Supplied
#: explicitly so the object does not depend on the ambient
#: `user.name`/`user.email` of whatever machine runs the campaign; a
#: checkout with none configured (CI is one) would otherwise fail at
#: `commit-tree`.
_IDENTITY = "o4b-reservation"


def normalise_remote(url: str) -> str:
    """`host/owner/repo`, from any of git's URL spellings."""

    u = url.strip().rstrip("/")
    for scheme in ("ssh://", "git+ssh://", "git://", "https://",
                   "http://"):
        u = u.removeprefix(scheme)
    u = u.removesuffix(".git").rstrip("/")
    head, _, rest = u.partition("/")
    if "@" in head:                       # user@host, incl. scp form
        head = head.split("@", 1)[1]
    if ":" in head:                       # scp form host:owner/repo
        host, _, tail = head.partition(":")
        head, rest = host, f"{tail}/{rest}" if rest else tail
    return f"{head}/{rest}".rstrip("/").lower()


def authority() -> str:
    """The remote's URL, once it is proved to BE the canonical
    repository and not a same-named local alias."""

    got = subprocess.run(["git", "remote", "get-url", REMOTE],
                         cwd=_REPO, input="", capture_output=True,
                         text=True)
    if got.returncode != 0:
        raise SystemExit(
            f"reservation: no `{REMOTE}` remote "
            f"({got.stderr.strip()}) -- the campaign cannot reach the "
            f"authority that serialises its streams")
    url = got.stdout.strip()
    if normalise_remote(url) != CANONICAL_AUTHORITY:
        raise SystemExit(
            f"reservation: `{REMOTE}` points at {url!r}, not "
            f"{CANONICAL_AUTHORITY} -- a fork or mirror has its own "
            f"empty reservation ref and would let the same streams be "
            f"drawn twice")
    return url


def _ls_remote(ref: str) -> str | None:
    """The object at `ref`, or None if it is free.

    A network or configuration failure is FATAL, never "assume free":
    the whole point of the ref is that no checkout may open the
    streams without consulting the one authority."""

    authority()
    probe = subprocess.run(["git", "ls-remote", REMOTE, ref],
                           cwd=_REPO, capture_output=True, text=True)
    if probe.returncode != 0:
        raise SystemExit(
            f"cannot read {REMOTE}/{ref}: {probe.stderr.strip()} -- "
            f"refusing to act without reaching the reservation "
            f"authority")
    line = probe.stdout.strip()
    return line.split()[0] if line else None


def held() -> str | None:
    """The object at THIS campaign's ref, or None."""

    return _ls_remote(REF)


def verify_o4_ref_retained() -> str:
    """O4's reservation is still exactly where the incident left it.

    Checked rather than assumed. Nothing in a digest manifest or a rev
    check can see a ref that was quietly deleted, and a cleared O4 ref
    would make the retired streams look free -- the one failure this
    program has already decided must never happen."""

    got = _ls_remote(RETAINED_REF)
    if got is None:
        raise SystemExit(
            f"{RETAINED_REF} is GONE -- it is retained permanently "
            f"and holds the O4 campaign's claim on retired streams; "
            f"refusing to run with it cleared")
    if got != RETAINED_OBJECT:
        raise SystemExit(
            f"{RETAINED_REF} holds {got}, not {RETAINED_OBJECT} -- it "
            f"has been rewritten, and the retired streams' claim is "
            f"no longer the one the incident recorded")
    return got


class ClaimUncertain(Exception):
    """The push was attempted and the outcome is not known to be clean.

    Raised instead of exiting, so the caller can record an incident
    (review R23). `claim()` pushes and then re-reads the ref, and a
    network failure BETWEEN those two steps leaves the ref created --
    the seeds spent by policy -- while the function never returns. An
    exit there produces the one thing this freeze exists to prevent:
    streams that are spent with no record saying so.

    `pushed` is what git reported. It is not a guarantee either way:
    a push that returns an error can still have been applied at the
    server. So the caller must treat any instance of this as SEEDS
    POSSIBLY SPENT -- the safe direction of the error is to over-
    record a claim, never to under-record one."""

    #: What git managed to tell us. THREE states, not two (review
    #: R24): the push can also fail to report at all -- an interrupt
    #: or an OS error inside `subprocess.run` itself, after the
    #: request has already reached the server. `None` is that case,
    #: and it is the one an earlier version left outside this
    #: exception entirely.
    def __init__(self, reason: str, pushed: bool | None,
                 obj: str | None, detail: str = "") -> None:
        super().__init__(reason)
        self.reason = reason
        self.pushed = pushed
        self.obj = obj
        self.detail = detail

    def as_record(self) -> dict:
        return {
            "reason": self.reason,
            "push_reported": {True: "success", False: "error",
                              None: "did not report"}[self.pushed],
            "attempt_object": self.obj,
            "detail": self.detail,
            "seeds_must_be_treated_as": "spent",
            "why": ("the push was attempted, so the ref may hold this "
                    "attempt whatever git reported -- including when "
                    "it reported nothing; recording the seeds as free "
                    "would be the unsafe direction"),
        }


def _make_commit(message: str) -> str:
    """A commit object over the empty tree, carrying `message`."""

    env = dict(os.environ)
    for role in ("AUTHOR", "COMMITTER"):
        env[f"GIT_{role}_NAME"] = _IDENTITY
        env[f"GIT_{role}_EMAIL"] = f"{_IDENTITY}@invalid"
    empty_tree = subprocess.run(
        ["git", "hash-object", "-t", "tree", "--stdin"], cwd=_REPO,
        input="", capture_output=True, text=True, check=True
    ).stdout.strip()
    made = subprocess.run(
        ["git", "commit-tree", empty_tree, "-m", message], cwd=_REPO,
        input="", capture_output=True, text=True, env=env)
    if made.returncode != 0:
        raise SystemExit(
            f"reservation: `git commit-tree` failed: "
            f"{made.stderr.strip()}")
    return made.stdout.strip()


def probe_namespace() -> None:
    """Prove at PREFLIGHT that the reservation can actually be made.

    Push permission and the server's ref-namespace policy are only
    exercised by a real push, so without this the first thing a
    campaign could learn -- after clearing every other check -- is
    that it cannot claim its streams. `refs/o4b/` is a NEW namespace,
    which makes this more than a formality here.

    The probe uses a DIFFERENT ref and is deleted immediately; it
    observes nothing, and it cannot stand in for the claim itself,
    which must be made by the run that draws."""

    authority()
    probe_ref = REF + "-preflight-probe"
    obj = _make_commit("o4b preflight probe")
    push = subprocess.run(
        ["git", "push", "--force", REMOTE, f"{obj}:{probe_ref}"],
        cwd=_REPO, capture_output=True, text=True)
    if push.returncode != 0:
        raise SystemExit(
            f"preflight: cannot write {REMOTE}/{probe_ref}: "
            f"{push.stderr.strip()} -- the campaign could not claim "
            f"its streams, so it must not start")
    drop = subprocess.run(["git", "push", REMOTE, f":{probe_ref}"],
                          cwd=_REPO, capture_output=True, text=True)
    if drop.returncode != 0:
        print(f"warning: left {probe_ref} behind on {REMOTE}; delete "
              f"it manually ({drop.stderr.strip()})")


def claim(payload: dict) -> str:
    """Claim the campaign streams globally, before the first draw.

    Create-only semantics need care: pushing the freeze commit itself
    would be a no-op "everything up-to-date" against an existing
    reservation at the same SHA, i.e. a silent second claim. The
    pushed object is therefore unique to this ATTEMPT, so
    `--force-with-lease=<ref>:` ("the ref must not exist") rejects
    every later attempt.

    UNIQUE BY A NONCE, NOT BY THE PAYLOAD (review R12). The payload is
    the freeze rev, the manifest digest and the seeds -- identical
    across attempts by construction. The identity is fixed on purpose.
    Git's commit timestamp has one-second resolution. So two runs on
    one host reaching this within the same second would build the same
    empty tree, message, identity and timestamp, hence the SAME commit
    SHA -- and then the second push is "everything up-to-date", the
    `held() == obj` check passes for both, and two processes draw the
    same campaign seeds while each believes it holds the claim. That
    is the very no-op this comment set out to avoid, arrived at from
    the other direction."""

    verify_o4_ref_retained()
    if (already := held()) is not None:
        raise SystemExit(
            f"reservation: {REF} is already held by {already} -- a "
            f"campaign has already opened these streams from some "
            f"checkout; the seeds are spent whether or not that run "
            f"published a result")
    attempt = {**payload, "attempt_nonce": os.urandom(16).hex()}
    obj = _make_commit(
        json.dumps(attempt, ensure_ascii=False, sort_keys=True))
    # FROM HERE ON THE OUTCOME IS UNCERTAIN, NOT MERELY UNSUCCESSFUL.
    # Everything above this line reads; this line writes, and a write
    # whose reply is lost has still happened at the server.
    try:
        push = subprocess.run(
            ["git", "push", f"--force-with-lease={REF}:", REMOTE,
             f"{obj}:{REF}"],
            cwd=_REPO, capture_output=True, text=True)
    except BaseException as exc:
        # The call did not come back -- an interrupt, an OS error --
        # and the request may already have reached the server. This
        # is the same failure class as a lost reply, so it takes the
        # same exit (review R24). `KeyboardInterrupt` included: a
        # deliberate stop here still cannot un-send the push.
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
            f"{REF} holds {confirmed}, not this attempt's object -- "
            f"refusing to draw on streams reserved by someone else",
            pushed=True, obj=obj, detail=f"held={confirmed}")
    return obj
