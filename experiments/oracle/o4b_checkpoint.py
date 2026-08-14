"""Progress that survives the run being killed.

O4 lost twelve hours of G1 and G2 statistics to an abort, and the
direct cause was not the abort -- it was that nothing intermediate had
been written down. The numbers existed only in a process that stopped.

So O4b writes a checkpoint at four points: the end of G3b, every G1
chunk, G1 completion, G2 completion. Each one carries enough to say
which freeze produced it and where the stream had got to.

TWO PROPERTIES MATTER, AND THEY ARE DIFFERENT.

ATOMIC. The write goes to a temporary file in the same directory and
then `os.replace`s the target. A checkpoint interrupted halfway must
not destroy the last good one -- the failure mode this exists to
prevent is exactly a run that dies at an awkward moment.

NOT A VERDICT. Every checkpoint is stamped `partial` and
`non_verdict`, and lives under a different name from the result
artifact. A checkpoint is a record of progress; reading one as a
result would publish a number the stopping rule never sanctioned.
That is why the stamps are written by this module rather than passed
in: a caller cannot forget them, and cannot claim otherwise.

WHY `os.replace` AND NOT `Path.rename`. On Windows `rename` fails when
the target exists, which would make every checkpoint after the first a
delete-then-write -- a window in which there is no checkpoint at all.
`os.replace` is atomic on both platforms.

Checkpoints are deliberately NOT write-once. They are meant to be
overwritten; the incident artifacts of the abort paths are the
write-once ones, and they are a different thing.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

#: The four points the freeze names. A checkpoint at any other stage
#: is a caller inventing a stage.
STAGES = ("g3b", "g1_chunk", "g1_complete", "g2_complete")

#: Keys every checkpoint must carry. Named here so a stage that
#: forgets one fails at the write rather than at the recovery, when
#: the run that could have supplied it is gone.
REQUIRED = (
    "freeze_sha",        # which commit
    "manifest_digest",   # which frozen configuration
    "seed",              # which stream
    "rng_position",      # how far into it
    "samples",           # how many points processed
    "statistics",        # the running accumulators
    "budget",            # calls and wall spent
)

#: Keys this module writes and a payload may not supply.
#:
#: Refused rather than overridden (review R3). Ordering the spread so
#: the stamps win protects `partial`, but a payload carrying `stage`
#: is a caller that thinks it is describing something else -- a reused
#: run-state dict whose `stage` is the RUN stage `g1`, not the
#: checkpoint point `g1_chunk`. Silently winning writes the right file
#: and hides the confusion; silently losing writes a `stage` that is
#: not even in `STAGES`, and a resume reading it continues from the
#: wrong place. Neither is a thing to discover during a recovery, so
#: the collision is an error at the write.
RESERVED = ("kind", "stage", "partial", "non_verdict", "why")


def write(path: Path, stage: str, payload: dict) -> Path:
    """Replace `path` with a checkpoint, atomically.

    Returns the path written, so a caller can record it in the
    incident without re-deriving the name."""

    if stage not in STAGES:
        raise ValueError(
            f"{stage!r} is not one of the frozen checkpoint stages "
            f"{STAGES}; the freeze names where progress is recorded")
    missing = [k for k in REQUIRED if k not in payload]
    if missing:
        raise ValueError(
            f"checkpoint at {stage!r} is missing {missing}: a "
            f"checkpoint that cannot say which freeze, which stream "
            f"and how far cannot be resumed from or audited")
    clashing = [k for k in RESERVED if k in payload]
    if clashing:
        raise ValueError(
            f"checkpoint payload supplies {clashing}, which this "
            f"module writes: `stage` is the checkpoint point, not the "
            f"run stage, and `partial`/`non_verdict` are what keep a "
            f"progress record from reading as a result")

    record = {
        "kind": "checkpoint",
        "stage": stage,
        **payload,
        "partial": True,
        "non_verdict": True,
        "why": ("progress record, not a result: the stopping rule has "
                "not fired, so no statistic here is a verdict"),
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent),
                               prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
            json.dump(record, f, ensure_ascii=False, indent=1,
                      sort_keys=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return path


def read(path: Path) -> dict:
    """The last checkpoint, or `{}` if there is none.

    A missing file is not an error: the run may have died before the
    first chunk, and the caller has to be able to tell that apart from
    a corrupt one, which raises."""

    path = Path(path)
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def is_verdict(record: dict) -> bool:
    """Never true for a checkpoint. Exists so the composition step can
    assert it against whatever it is handed, rather than trusting the
    filename it was handed it under."""

    return not (record.get("partial") and record.get("non_verdict"))
