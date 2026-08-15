"""Recover the O4b result from the preserved completed statistics.

The O4b campaign ran once, generated its data, and computed both gates
to their frozen sample sizes. What failed was a single line of exit
provenance wiring AFTER that -- `main`'s claim callback returned `None`,
so the pre-publish re-verification never received the object it was
supposed to check and the process raised before writing the result
(docs/prereg/p14_o4b_incident.md). The gate statistics survived, in the
incident and the checkpoint, as non-verdict partials.

A recovery is NOT a re-run and NOT a new observation. It re-applies the
SAME frozen decision functions to the SAME preserved sufficient
statistics and republishes the verdict the run computed but never wrote.
No seed is drawn, no solver is called, no point is resampled, and no
gate rule is touched. Reading a preserved artifact after the fact is not
contamination; changing an analysis rule after seeing the data would be,
and this does neither -- the gates, sample sizes and analysis formulae
were frozen before the run and are applied here byte-for-byte.

THE INTEGRITY GATE. The recovery recomputes G1 and G2 independently and
refuses to publish unless the recomputation reproduces the preserved
verdict BIT-FOR-BIT. A mismatch would mean the preserved record and the
frozen functions disagree -- a corrupted or tampered artifact -- and the
safe response is to publish nothing, not to paper over it.

WHAT IS RECOMPUTED, AND WHAT IS ONLY CARRIED. G1's interval is a pure
function of `(n, mean, var)`, all three preserved at full precision, so
it is reproduced exactly. G2's decisive quantity -- the Clopper-Pearson
upper leak bound -- is a pure function of `(leaks, n)`, also preserved,
so the `concordant` verdict (`upper <= budget`) is reproduced exactly.
G2's lower bound needs the G2 accumulator's own mean and variance, which
the checkpoint did not serialise (it carries the G1 stream's; incident
defect 5c); it does NOT enter the verdict, so its preserved value is
carried forward and labelled, never recomputed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import empirical_bernstein as eb  # noqa: E402
import exact_binomial as xb  # noqa: E402
import o4_sizing as sz  # noqa: E402
import o4b_reservation as reservation  # noqa: E402
import o4b_volume_audit as run  # noqa: E402
import probe_seed_ledger as ledger  # noqa: E402

_PREREG = _REPO / "docs" / "prereg"
_INCIDENT = _PREREG / "p14_o4b_incident.json"
_CHECKPOINT = _PREREG / "p14_o4b_checkpoint.json"
_EXECUTED_MANIFEST = _PREREG / "p14_o4b_executed_freeze_manifest.json"
_RESULT = _PREREG / "p14_o4b_results.json"

#: The approved freeze the run executed from (the freeze branch head,
#: not the merge commit) and the manifest that pins the surface it
#: verified against. Both are asserted against the preserved records.
FREEZE_SHA = "715865abc684224785e71a3130e17c50db35f947"
EXECUTED_DIGEST = (
    "cec650b9391af0fc11e4b6bb94455cdbc1a037c18ecec83149d0d4693e7d7be2")

#: The incident and checkpoint as committed in PR #76. Pinned so the
#: recovery inputs are authenticated against the immutable record, not
#: merely trusted (review PR #77 R1): without this, a self-consistent
#: edit to one file's G2 block -- its `n` and the matching `status`,
#: `leak_upper_abs`, `budget_abs` together -- would pass a recomputation
#: that reads both sides from that same file, and publish a wrong verdict
#: as a write-once result. With the blob pinned AND the two files
#: cross-checked, a tamper must edit the file, its independent second
#: copy, and this reviewed pin, all consistently and all in the diff.
INCIDENT_SHA256 = (
    "8106b16f0d03efe1acc81941f7ca3149cb8ddee0fbaaf00cfccb5c8376d4cc75")
CHECKPOINT_SHA256 = (
    "5dd6de7eea25c8384329b26ebf9f61c9dae51ee451c6d5c5717dbd463dc6679a")


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


class RecoveryError(Exception):
    """The preserved record and the frozen functions disagree, or the
    frozen surface drifted. No result is published."""


#: The files the recovery's correctness actually rests on, pinned
#: against the EXECUTED manifest. The verdict is a function of the
#: empirical-Bernstein interval, the Clopper-Pearson bound and the
#: o4_sizing constants; the publish-time reservation re-verify is the
#: reservation module. These are checked byte-for-byte -- NOT the whole
#: protocol surface (as an earlier version did). The campaign runner,
#: the stage machine and the checkpoint format do not enter the recovered
#: number, so ordinary maintenance to them (a fixed claim return, a stage
#: budget label, an explicit failure_point) must not break a recovery
#: whose inputs and decision functions are unchanged. The incident and
#: checkpoint DATA is authenticated separately, by blob pin and
#: cross-check in `_authenticate`.
_DECISION_FILES = (
    "experiments/oracle/empirical_bernstein.py",
    "experiments/oracle/exact_binomial.py",
    "experiments/oracle/o4_sizing.py",
    "experiments/oracle/o4b_reservation.py",
)


def verify_frozen_surface() -> dict:
    """The decision functions the recovery uses are the frozen ones.

    Verified against the EXECUTED manifest -- the surface the campaign
    actually ran against -- not the live one, which re-pins as the
    protocol surface is maintained. Only the files the recovered verdict
    and its reservation re-verify depend on are pinned (`_DECISION_FILES`);
    the ledger's one permitted change (the two O4b seeds now OBSERVED) is
    checked below."""

    manifest = _load(_EXECUTED_MANIFEST)
    if run._sha256(_EXECUTED_MANIFEST) != EXECUTED_DIGEST:
        raise RecoveryError(
            f"executed manifest digest drifted from {EXECUTED_DIGEST}")

    drifted = []
    for rel in _DECISION_FILES:
        want = manifest["files"].get(rel)
        if want is None:
            raise RecoveryError(
                f"{rel} is not pinned in the executed manifest")
        if run._sha256(_REPO / rel) != want:
            drifted.append(rel)
    if drifted:
        raise RecoveryError(
            f"frozen decision surface drifted on {drifted} -- refusing "
            f"to recover with functions that are not the frozen ones")

    # the ledger's permitted change, made explicit: the TWO O4b seeds
    # are retired. Checked per-seed, not by requiring the whole fresh
    # pool empty (review PR #77 R2): FRESH_PROBE_SCALARS is the
    # program-wide active ledger, so a later campaign freezing its own
    # new scalar is a normal operating state that must not fail this
    # recovery -- only the O4b seeds' own status is the recovery's
    # business.
    for name, seed in (("o4b_g1_audit", 40_000_401),
                       ("o4b_g2_leakage", 40_000_411)):
        if name in ledger.FRESH_PROBE_SCALARS:
            raise RecoveryError(
                f"{name} is still FRESH; the O4b seeds must be retired "
                f"before their result is recovered")
        if ledger.OBSERVED_PROBE_SCALARS.get(name) != seed:
            raise RecoveryError(
                f"{name} is not retired at {seed} in OBSERVED")

    env, locked = run.environment(), manifest["environment"]
    drift = {k: (locked[k], env.get(k)) for k in locked
             if env.get(k) != locked[k]}
    if drift:
        raise RecoveryError(
            f"environment drift {drift} -- the rounding instrument is "
            f"part of the frozen apparatus the functions ran under")
    return manifest


def recompute_g1(stats: dict) -> dict:
    """G1's verdict from the preserved `(n, mean_z, var_z)`.

    `eb.interval` uses only those three; `_statistics` stored `acc.mean`
    and `acc.var` at full precision, so the empirical-Bernstein interval,
    its rescaling, the identified discrepancy and the concordance test
    reproduce the run's numbers exactly."""

    n, mean, var = stats["n"], stats["mean_z"], stats["var_z"]
    delta = sz.DELTA_G1
    interval = eb.EBInterval(
        n=n, mean=mean, var=var,
        half_width=eb.half_width(n, var, delta), delta=delta)
    lo, hi = interval.rescaled(sz.SCALE)
    gap = (lo - sz.V_HI, hi - sz.V_LO)
    band = sz.TAU * sz.V_REF
    if gap[0] >= -band and gap[1] <= band:
        status = "concordant"
    elif gap[1] < -band or gap[0] > band:
        status = "discordant"
    else:
        status = "inconclusive"
    return {"status": status, "n": n, "v_s1_lo": lo, "v_s1_hi": hi,
            "identified_discrepancy": list(gap), "band_abs": band}


def recompute_g2(preserved: dict) -> dict:
    """G2's verdict from the preserved `(leaking_points, n)`.

    The Clopper-Pearson upper bound and the budget are pure functions of
    the preserved counts and the frozen config, so the `concordant`
    decision (`upper <= budget`) is reproduced exactly. The lower bound
    is not recomputable from what was serialised (incident defect 5c) and
    does not enter the verdict; the preserved value is carried, labelled,
    never recomputed."""

    leaks, n = preserved["leaking_points"], preserved["n"]
    alpha = sz.ALPHA_G2
    upper = sz.DT * sz.B_OUT * xb.cp_upper(leaks, n, alpha)
    allowed = sz.LEAK_BUDGET * sz.V_REF
    lower_carried = preserved["leak_lower_abs"]
    if upper <= allowed:
        status = "concordant"
    elif lower_carried > allowed:
        status = "discordant"
    else:
        status = "inconclusive"
    return {
        "status": status, "n": n, "leaking_points": leaks,
        "leak_upper_abs": upper, "budget_abs": allowed,
        "leak_lower_abs": lower_carried,
        "leak_lower_abs_provenance": (
            "carried from the preserved record, NOT independently "
            "recomputed: the g2_complete checkpoint serialised the G1 "
            "stream's statistics, not G2's (incident defect 5c). It does "
            "not enter the concordant verdict, which rests on "
            "leak_upper_abs <= budget_abs."),
    }


def _assert_reproduces(kind: str, got: dict, preserved: dict,
                       keys: tuple[str, ...]) -> None:
    """The integrity gate: bit-for-bit, or no result is published."""

    for k in keys:
        g, p = got[k], preserved[k]
        if isinstance(g, list):
            ok = list(g) == list(p)
        else:
            ok = g == p
        if not ok:
            raise RecoveryError(
                f"{kind}.{k} recomputed as {g!r} but the preserved "
                f"record holds {p!r} -- the frozen functions and the "
                f"preserved statistics disagree; publishing nothing")


def _authenticate() -> tuple[dict, dict]:
    """The preserved inputs are the committed record, and the two files
    agree (review PR #77 R1).

    Two independent defences against a tampered input:

      * BLOB PIN. Each file's sha256 is pinned to what PR #76 committed,
        so any edit to either changes a hash checked here.

      * CROSS-CHECK. The incident and the checkpoint were written at
        different points in the run and carry independent copies of the
        same facts -- the G2 block, the G1 sufficient statistics, the
        lineage. They must agree, so a single-file edit that keeps one
        file self-consistent still diverges from the other and is caught
        here, before any recomputation reads from it."""

    for path, want in ((_INCIDENT, INCIDENT_SHA256),
                       (_CHECKPOINT, CHECKPOINT_SHA256)):
        got = run._sha256(path)
        if got != want:
            raise RecoveryError(
                f"{path.name} sha256 {got} != committed {want} -- the "
                f"preserved input is not the record PR #76 committed")

    incident, checkpoint = _load(_INCIDENT), _load(_CHECKPOINT)

    for label, rec in (("incident", incident), ("checkpoint", checkpoint)):
        if rec["freeze_sha"] != FREEZE_SHA:
            raise RecoveryError(
                f"{label}.freeze_sha {rec['freeze_sha']} != {FREEZE_SHA}")
        if rec["manifest_digest"] != EXECUTED_DIGEST:
            raise RecoveryError(
                f"{label}.manifest_digest != executed manifest digest")

    # the two independent copies of the same facts must agree
    g1_partial = incident["preserved"]["g1_partial"]
    stats = checkpoint["statistics"]
    if (g1_partial["n"], g1_partial["mean_z"], g1_partial["var_z"]) != (
            stats["n"], stats["mean_z"], stats["var_z"]):
        raise RecoveryError(
            "the incident's G1 sufficient statistics and the checkpoint's "
            "disagree -- one of the preserved files was altered")
    if incident["preserved"]["completed_gates"]["g2"] != stats.get("g2"):
        raise RecoveryError(
            "the incident's G2 block and the checkpoint's copy disagree "
            "-- one of the preserved files was altered")
    return incident, checkpoint


def recover() -> dict:
    """Rebuild the O4b result. Raises RecoveryError rather than publish
    anything that does not reproduce the preserved verdict exactly."""

    incident, checkpoint = _authenticate()
    obj = incident["reservation_object"]
    if obj is None or incident["reservation_claimed"] is not True:
        raise RecoveryError(
            "the incident does not record a claimed reservation object")

    verify_frozen_surface()

    preserved_gates = incident["preserved"]["completed_gates"]
    # G1: recompute from the CHECKPOINT's statistics, compare against the
    # INCIDENT's gate -- two independently-written files, so the input
    # and the target cannot both be moved by editing one of them.
    g1 = recompute_g1(checkpoint["statistics"])
    _assert_reproduces("g1", g1, preserved_gates["g1"],
                       ("status", "n", "v_s1_lo", "v_s1_hi",
                        "identified_discrepancy", "band_abs"))
    # G2: recompute from the CHECKPOINT's g2 counts, compare against the
    # INCIDENT's g2 gate. `_authenticate` has already asserted the two
    # g2 blocks are identical, so this is a genuine cross-file check
    # rather than a file agreeing with itself.
    g2 = recompute_g2(checkpoint["statistics"]["g2"])
    _assert_reproduces("g2", g2, preserved_gates["g2"],
                       ("status", "n", "leaking_points",
                        "leak_upper_abs", "budget_abs"))

    return {
        "kind": "results",
        "run_kind": "recovered_completed_campaign",
        "freeze_sha": FREEZE_SHA,
        "manifest_digest": EXECUTED_DIGEST,
        "seeds": incident["seeds"],
        "recovery": {
            "source": [_INCIDENT.name, _CHECKPOINT.name],
            "method": ("the frozen decision functions (empirical "
                       "Bernstein, Clopper-Pearson, o4_sizing constants) "
                       "re-applied to the preserved sufficient "
                       "statistics"),
            "reproduced_preserved_verdicts_bit_exact": True,
            "no_new_seed": True,
            "no_solver_call": True,
            "no_resampling": True,
            "no_gate_change": True,
            "seeds_status": (
                "OBSERVED -- the two streams were drawn and are retired "
                "(PR #76); a recovery is a reproduction of that observed "
                "run, never a new observation"),
            "why": (
                "the campaign completed both gates; only the exit "
                "provenance wiring failed, so the verdict the run "
                "computed is recovered from its own preserved statistics "
                "rather than re-run"),
        },
        "reservation": {
            "ref": reservation.REF,
            "authority": reservation.CANONICAL_AUTHORITY,
            "object": obj,
            "seeds_spent": True,
            "verified_at_recovery": True,
            "why": ("the streams were opened by this attempt and are "
                    "retired; the object is the one the ref held when "
                    "the claim was made and when this recovery published"),
        },
        "order": list(run.stages.STAGES),
        "g1": g1,
        "g2": g2,
        "availability": incident["preserved"]["availability"],
        "budget": incident["preserved"]["budget"],
        "environment": run.environment(),
    }


def main() -> None:
    result = recover()
    # Item 3 of the recovery contract: re-verify the remote reservation
    # immediately before publishing, with the object actually in hand --
    # the exact check the campaign's wiring defect skipped. A lost or
    # replaced ref raises here, and nothing is published over it.
    reservation.verify_still_held(result["reservation"]["object"])
    published = run.publish_write_once(_RESULT, result)
    if published:
        print(f"recovered O4b result published to {_RESULT.name}")
        print(f"  G1 {result['g1']['status']} "
              f"[{result['g1']['v_s1_lo']}, {result['g1']['v_s1_hi']}]")
        print(f"  G2 {result['g2']['status']} "
              f"leak_upper {result['g2']['leak_upper_abs']} "
              f"<= budget {result['g2']['budget_abs']}")
    else:
        print(f"{_RESULT.name} already exists -- the first recovery "
              f"stands; nothing overwritten")


if __name__ == "__main__":
    main()
