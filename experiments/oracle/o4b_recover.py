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

#: The one protocol-surface file that legitimately changed after the
#: run: the ledger moved the two O4b seeds FRESH -> OBSERVED. Every
#: OTHER file in the executed manifest must still match byte-for-byte.
_LEDGER_REL = "experiments/positive_control/probe_seed_ledger.py"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


class RecoveryError(Exception):
    """The preserved record and the frozen functions disagree, or the
    frozen surface drifted. No result is published."""


def verify_frozen_surface() -> dict:
    """The decision functions the recovery uses are the frozen ones.

    Verified against the EXECUTED manifest -- the surface the campaign
    actually ran against -- not the live one, which re-pinned the ledger
    when the seeds retired. Every file matches byte-for-byte except the
    ledger, whose only permitted change is that the two O4b seeds are now
    OBSERVED rather than FRESH."""

    manifest = _load(_EXECUTED_MANIFEST)
    if run._sha256(_EXECUTED_MANIFEST) != EXECUTED_DIGEST:
        raise RecoveryError(
            f"executed manifest digest drifted from {EXECUTED_DIGEST}")

    drifted = []
    for rel, want in manifest["files"].items():
        got = run._sha256(_REPO / rel)
        if got == want:
            continue
        if rel == _LEDGER_REL:
            continue                        # the intended, verified delta
        drifted.append(rel)
    if drifted:
        raise RecoveryError(
            f"frozen decision surface drifted on {drifted} -- refusing "
            f"to recover with functions that are not the frozen ones")

    # the ledger's permitted change, made explicit: the seeds are spent
    if ledger.FRESH_PROBE_SCALARS:
        raise RecoveryError(
            f"the ledger still lists fresh allocations "
            f"{ledger.FRESH_PROBE_SCALARS}; the O4b seeds must be retired")
    for name, seed in (("o4b_g1_audit", 40_000_401),
                       ("o4b_g2_leakage", 40_000_411)):
        if ledger.OBSERVED_PROBE_SCALARS.get(name) != seed:
            raise RecoveryError(
                f"{name} is not retired at {seed} in the ledger")

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


def recover() -> dict:
    """Rebuild the O4b result. Raises RecoveryError rather than publish
    anything that does not reproduce the preserved verdict exactly."""

    incident = _load(_INCIDENT)
    checkpoint = _load(_CHECKPOINT)

    # the preserved records agree on which run this is
    for label, rec in (("incident", incident), ("checkpoint", checkpoint)):
        if rec["freeze_sha"] != FREEZE_SHA:
            raise RecoveryError(
                f"{label}.freeze_sha {rec['freeze_sha']} != {FREEZE_SHA}")
        if rec["manifest_digest"] != EXECUTED_DIGEST:
            raise RecoveryError(
                f"{label}.manifest_digest != executed manifest digest")
    obj = incident["reservation_object"]
    if obj is None or incident["reservation_claimed"] is not True:
        raise RecoveryError(
            "the incident does not record a claimed reservation object")

    verify_frozen_surface()

    preserved_gates = incident["preserved"]["completed_gates"]
    g1 = recompute_g1(checkpoint["statistics"])
    _assert_reproduces("g1", g1, preserved_gates["g1"],
                       ("status", "n", "v_s1_lo", "v_s1_hi",
                        "identified_discrepancy", "band_abs"))
    g2 = recompute_g2(preserved_gates["g2"])
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
