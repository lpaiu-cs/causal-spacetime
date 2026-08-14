"""O4b sizing: the same O4 gate, re-costed for the G3 redesign.

This re-run repairs an instrumentation contract; it does not sharpen
anything. The oracle is NOT re-frozen: the existing O3 certified
interval and tau = 3.0% are used unchanged, and G1 keeps `N_G1` with
its certified worst-case power. Tightening the oracle would be a
different question -- new interval, new sample size, new cost -- and is
left to a separate precision arc.

`o4_sizing` is digest-pinned by the O4 freeze manifest and by the
executed snapshot the replay diagnostic verifies against, so it is
imported, never edited. Everything here is additive.

WHAT CHANGED, AND WHY IT COSTS MORE. The frozen G3 spent four solver
calls per cluster: two probes, two legs each. G3b spends more, for two
reasons that are both consequences of the redesign rather than
overheads bolted onto it:

  * the predicate's own `(t1, err1, t2, err2)` must be computed at
    `dpsi`, not at `theta`, before any probe can be placed -- two calls
    per CANDIDATE examined, not per cluster accepted;
  * there are three probes now, not two, each with two legs -- six
    calls per fully-testable cluster.

So the cost is

    G3b calls = 2 * scanned_positive_candidates
              + 6 * fully_testable_clusters

and the first term is charged on every candidate the scan looks at,
including the ones it rejects. That is what makes the total depend on
availability, which is the thing no freeze may assume.

THE CAP IS A COMPLETION BUDGET, NOT A POWER STATEMENT. If availability
were poor the scan could run to the end of the G1 sample, and the call
cap would then fire before the scan does -- the hard worst case is
about 107.5M calls against an 80M cap. The freeze does not resize the
caps to cover that: reserving for 26.2M rejected candidates would be
wildly out of proportion to what the census actually saw. Instead the
caps stay where they are, the freeze STATES that they may fire first,
and firing yields an incident and `INCONCLUSIVE`. Raising a cap during
or after a run remains forbidden.

Run:  python experiments/oracle/o4b_sizing.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "positive_control"))

import o4_g3_redesign as g3  # noqa: E402
import o4_sizing as base  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o4b_sizing.json"

#: Unchanged from O4, and unchanged deliberately. The O3 certified
#: interval is reused as it stands; the full-width floor it implies is
#: 1.999426%, so tau = 3.0% remains valid, and `N_G1` keeps the
#: worst-case power certified on the continuum at that size.
TAU = base.TAU
N_G1 = base.N_G1

#: The G3 unit costs (design section 3, prereg re-opening section 2.4).
STATE_CALLS_PER_CANDIDATE = 2      # (t1, err1), (t2, err2) at dpsi
PROBE_CALLS_PER_CLUSTER = 6        # three probes, two legs each

#: The scan's bound, both as the reference the runner reads and as the
#: number that reference resolves to today. The prereg re-opening asks
#: for both: the reference so a re-derived G1 size moves the bound with
#: it, the number so the freeze manifest records what was actually in
#: force.
SCAN_CAP_REALIZED = 26_200_000

#: Frozen caps, carried over from O4 unchanged. They are a COMPLETION
#: budget: enough for the run the census makes likely, not enough for
#: the worst case the scan rule permits. See `cap_contract()`.
MAX_CALLS = 80_000_000
MAX_WALL_S = 86_400.0

#: The committed census artifact, the only lineage the projection is
#: allowed to lean on. An earlier draft quoted a positive rate taken
#: from a run log that was never committed and that no total used --
#: a number with no provenance and no purpose, so it is gone.
_CENSUS = _REPO / "docs" / "prereg" / "p14_o4_replay_diagnostic.json"


def scan_cap() -> int:
    """The end of the frozen G1 sample. Delegates, never transcribes."""

    return g3.scan_cap()


def g3a_calls(tol: float = s1.DEFAULT_TOL) -> dict:
    """The G3a preflight's solver cost, MEASURED END TO END.

    An earlier draft reckoned this: it counted the table's bisections
    for real but added the row calls arithmetically, because `run_g3a`
    did not exist yet. It came out 72 calls low -- `_family_of` costs
    one call per case that the arithmetic forgot, the two
    construction-unavailable rows cost nothing rather than one each,
    and the determinism probe costs four. Small, and wrong in the
    direction a reckoning is always wrong in: it counted what someone
    remembered.

    So the real preflight is now executed under the meter and the
    number is whatever it actually spends."""

    import o4b_g3a as g3a
    import o4b_meter as meter
    from o4b_budget import Budget

    probe = Budget(max_calls=10 ** 9, max_wall_s=10 ** 9)
    with meter.metered(probe):
        preflight = g3a.run_preflight(tol)
    if not preflight["passed"]:
        raise SystemExit(
            f"G3a preflight failed {preflight['failed_conditions']}: "
            f"refusing to size a stage whose wrapper contract does not "
            f"hold -- the run would end INVALID before any seed")
    result = preflight["table"]
    return {
        "cases_resolved": result["cases"],
        "cases_unreachable": len(result["unreachable"]),
        "cases_unreachable_are_the_named_ones": (
            not result["unexpected_unreachable"]
            and not result["unexpectedly_reachable"]),
        "calls": probe.reserved,
        "completed": probe.completed,
        "construction_unavailable_rows":
            result["construction_unavailable"],
        "ulp_distance_total": result["ulp_distance_total"],
        "search_comparisons_total":
            result["search_comparisons_total"],
        "distance_is_not_cost": result["distance_is_not_cost"],
        "host_dependent": ["ulp_distance_total",
                           "search_comparisons_total"],
        "why_host_dependent": (
            "the solver reaches libm for sin/cos/acos, so `t_min` and "
            "`err` differ in their last bits across platforms and the "
            "ulp distance to a satisfying placement moves with them. "
            "CI measured 8,111 where this host measured 8,120. The "
            "CALL counts do not move, so the cost model and the "
            "contract are unaffected -- only these two diagnostics "
            "are, and they must never be pinned as constants"),
        "preflight_conditions": preflight["conditions"],
        "measures": ("run_preflight, the single authority; the stage "
                     "and this measurement read the same verdict"),
        "families_covered": result["families_covered"],
        "basis": ("end-to-end measurement of run_preflight under the "
                  "meter"),
        "includes": ["the table's branch-edge bisections",
                     "one family probe per case",
                     "rows A/B/C where constructible",
                     "the solver determinism probe"],
        "row_d_costs_nothing": ("row D proves the predicate never "
                                "consults the solver, so it spends no "
                                "calls by construction"),
    }


def census_eligibility() -> dict:
    """The committed census row the projection rests on.

    The projection assumes the scan stops at `K_G3B` -- i.e. that
    essentially every positive candidate is fully testable. That
    assumption comes from the census's eligibility count at the frozen
    eta, and it is READ from the artifact rather than transcribed, so a
    changed artifact changes the projection instead of silently
    contradicting it.

    It remains a design input, and a partial one. The census measured
    eligibility and the lower probe's reach; it did NOT measure the
    realized-margin or nudge conditions that `fully-testable` also
    requires, because no probe was constructed. So the projection it
    supports is optimistic as well as non-probabilistic -- which is why
    the budget says `census_projected`, never `expected`."""

    art = json.loads(_CENSUS.read_text(encoding="utf-8"))
    rows = art["findings"]["census"]["eligibility_by_eta"]
    row = next(r for r in rows if r["eta"] == g3.ETA)
    return {
        "source": _CENSUS.relative_to(_REPO).as_posix(),
        "eta": row["eta"],
        "eligible_clusters": row["eligible_clusters"],
        "clusters_probed": art["findings"]["clusters_probed"],
        "lower_probe_unreachable_among_eligible":
            row["eligible_and_lower_probe_t_x_negative"],
        "measured_here": ["eligibility at the frozen eta",
                          "the lower probe's reach among eligible"],
        "not_measured_here": [
            "the realized-margin check", "the nudge search",
            "therefore the fully-testable rate itself"],
        "is": ("a design input from a fixed stress set on a retired "
               "stream; NOT a probabilistic expectation about the new "
               "campaign's stream and NOT a passing condition"),
    }


def g3b_calls(scanned_positive: int, fully_testable: int) -> int:
    """`2 * scanned + 6 * fully_testable`, the redesign's G3 cost."""

    return (STATE_CALLS_PER_CANDIDATE * scanned_positive
            + PROBE_CALLS_PER_CLUSTER * fully_testable)


def budget(tol: float = s1.DEFAULT_TOL) -> dict:
    """Expected and hard-worst call budgets, itemised."""

    n_g2 = base.g2_points()
    g3a = g3a_calls(tol)
    g1_calls, g2_calls = 2 * N_G1, 2 * n_g2

    # census-projected: the census saw every candidate eligible, so
    # the scan is projected to stop at K. A projection, not an
    # expectation -- see `census_eligibility`.
    projected_scanned = g3.K_G3B
    projected = g3b_calls(projected_scanned, g3.K_G3B)

    # hard worst: the scan reaches the end of the G1 sample. Only
    # POSITIVE candidates are charged the state calls, but nothing
    # committed bounds that fraction, so the bound charges every G1
    # point. Conservative in the direction a cap has to be.
    worst = g3b_calls(scan_cap(), g3.K_G3B)

    total_projected = g1_calls + g2_calls + g3a["calls"] + projected
    total_worst = g1_calls + g2_calls + g3a["calls"] + worst
    #: `MS_PER_POINT` is named for milliseconds but carries SECONDS per
    #: point (the frozen module divides it by 3600 to get hours); one
    #: call is half a point.
    s_per_call = base.MS_PER_POINT / 2.0
    return {
        "items": {
            "g1": g1_calls,
            "g2": g2_calls,
            "g3a": g3a["calls"],
            "g3b_state_census_projected":
                STATE_CALLS_PER_CANDIDATE * projected_scanned,
            "g3b_probes": PROBE_CALLS_PER_CLUSTER * g3.K_G3B,
            "g3b_state_worst":
                STATE_CALLS_PER_CANDIDATE * scan_cap(),
        },
        "g3a_detail": g3a,
        "worst_case_is_an_upper_bound_because": (
            "the state calls are charged only on positive candidates, "
            "but nothing committed bounds that fraction and no rate is "
            "a guarantee about a new stream, so the bound charges "
            "every G1 point"),
        "total_census_projected": total_projected,
        "census_projected_assumes": (
            "fully-testable rate = 1, so the scan stops at K_G3B. The "
            "census measured ELIGIBILITY and the lower probe's reach, "
            "not the realized-margin and nudge conditions that "
            "fully-testable also requires, so this is a "
            "census-informed OPTIMISTIC projection and not a measured "
            "rate. If fully-testable falls short of eligible, the scan "
            "runs longer and the total rises toward the hard bound"),
        "total_hard_worst": total_worst,
        "o4_total_for_comparison": (2 * base.N_G1 + 2 * n_g2
                                    + 4 * base.G3_CLUSTERS),
        "hours_census_projected": total_projected * s_per_call / 3600.0,
        "hours_hard_worst": total_worst * s_per_call / 3600.0,
        "caps": {"calls": MAX_CALLS, "wall_s": MAX_WALL_S},
        "cap_ratio_census_projected": MAX_CALLS / total_projected,
        "cap_ratio_hard_worst": MAX_CALLS / total_worst,
    }


#: The frozen outcome table. `INVALID` means the instrument is not
#: behaving as its contract says, wherever that is found; it is NOT
#: G3a's alone, because the redesign already declares a G3b probe
#: mismatch a contract failure. `INCONCLUSIVE` means the run could not
#: reach a decision with the resources or the sample it had.
OUTCOMES = (
    {"condition": "G3a case fails (wrong tri-state, or the recovery "
                  "or determinism check fails)",
     "outcome": "INVALID",
     "when": "before any fresh seed is touched"},
    {"condition": "a fully-testable G3b probe returns None, or "
                  "bool(a and b) != want",
     "outcome": "INVALID",
     "when": "the redesign calls this a contract failure and requires "
             "mismatch 0"},
    {"condition": "a single candidate is construction-unavailable",
     "outcome": "no stage outcome",
     "when": "tallied by reason and the scan CONTINUES; one candidate "
             "that cannot be probed says nothing about the stage"},
    {"condition": "fewer than K_G3B fully-testable clusters by the end "
                  "of the scan",
     "outcome": "INCONCLUSIVE",
     "when": "the sample was short, not the instrument wrong"},
    {"condition": "the call or wall cap binds before K_G3B is reached",
     "outcome": "INCONCLUSIVE",
     "when": "a resource limit, not an instrument violation"},
    {"condition": "provenance failure, or a solver exception",
     "outcome": "no scientific verdict",
     "when": "write-once incident only"},
)


def cap_contract() -> dict:
    """What happens when the completion budget runs out.

    Frozen as an outcome, not as a licence to spend more."""

    b = budget()
    return {
        "kind": "completion budget, not a power statement",
        "may_fire_before_the_scan_ends": b["cap_ratio_hard_worst"] < 1.0,
        "on_firing": "write-once incident, then INCONCLUSIVE",
        "not_invalid": ("a cap is a resource limit, not an instrument "
                        "violation. INVALID belongs to contract "
                        "failures wherever they occur -- G3a's cases "
                        "AND a G3b probe that comes back undecided or "
                        "boolean-mismatched"),
        "auto_raise": "forbidden, during or after the run",
        "charged_before_spending": (
            "the counter and the clock start before G3a and every "
            "call or bounded batch is reserved first; a cap checked "
            "after the calls were made is not a cap"),
        "why_not_resized": (
            "covering the hard worst case would reserve for 26.2M "
            "rejected candidates, wildly out of proportion to what the "
            "census saw; the honest move is to state that the cap can "
            "bind and to end INCONCLUSIVE when it does"),
        "outcomes": list(OUTCOMES),
    }


def summary() -> dict:
    b = budget()
    return {
        "stage": ("O4b sizing: the O4 gate re-costed for the G3 "
                  "redesign, with the oracle and tau unchanged"),
        "not_a_precision_change": (
            "This re-run uses the existing O3 certified interval and "
            "tau = 3.0% unchanged, with no O3' re-freeze. Its purpose "
            "is to re-run the same O4 gate with a repaired G3 "
            "instrumentation contract, not to improve precision."),
        "inherited": {
            "tau": TAU, "n_g1": N_G1,
            "delta_g1_per_side": base.DELTA_G1,
            "alpha_g2_per_end": base.ALPHA_G2,
            "leak_budget_frac": base.LEAK_BUDGET,
            "n_g2": base.g2_points(),
            "oracle": {"v_lo": base.V_LO, "v_hi": base.V_HI,
                       "v_ref": base.V_REF},
            "tau_floor_full_width": (base.V_HI - base.V_LO) / base.V_REF,
            "power_lower_bound":
                base.certify_power(N_G1)["power_lower_bound"],
        },
        "g3": {
            "eta": g3.ETA,
            "outward_search_cap": ("none: the search ends on "
                                   "representability, not a budget"),
            "n_avail": g3.N_AVAIL, "k_g3b": g3.K_G3B,
            "alpha_g3b": g3.ALPHA_G3B,
            "cp_upper_zero_mismatch": base.g3_upper(g3.K_G3B,
                                                    g3.ALPHA_G3B),
            "cp_is": "characterisation only, never an accuracy gate",
            "scan_cap_reference": "o4_g3_redesign.scan_cap()",
            "scan_cap_realized": SCAN_CAP_REALIZED,
            "cost_model": ("2 * scanned_positive_candidates + 6 * "
                           "fully_testable_clusters"),
            "state_calls_per_candidate": STATE_CALLS_PER_CANDIDATE,
            "probe_calls_per_cluster": PROBE_CALLS_PER_CLUSTER,
        },
        "census_design_input": census_eligibility(),
        "budget": b,
        "cap_contract": cap_contract(),
    }


def main() -> None:
    s = summary()
    b, g = s["budget"], s["g3"]
    print(f"tau = {TAU}  (floor "
          f"{s['inherited']['tau_floor_full_width']:.6%}), "
          f"n_g1 = {N_G1:,}, power >= "
          f"{s['inherited']['power_lower_bound']:.5f}")
    print(f"scan cap = {g['scan_cap_realized']:,} "
          f"(via {g['scan_cap_reference']} = {scan_cap():,})")
    print(f"\nG3a: {b['g3a_detail']['cases_resolved']} cases, "
          f"{b['g3a_detail']['cases_unreachable']} unreachable, "
          f"{b['items']['g3a']:,} calls")
    print("G3b cost = 2 * scanned + 6 * fully_testable")
    for name in ("g1", "g2", "g3a", "g3b_state_census_projected",
                 "g3b_probes"):
        print(f"  {name:<22} {b['items'][name]:>14,}")
    print(f"  {'TOTAL census-projected':<26} "
          f"{b['total_census_projected']:>14,}"
          f"   ({b['hours_census_projected']:.2f} h)")
    print(f"  {'(O4 was)':<22} "
          f"{b['o4_total_for_comparison']:>14,}")
    print(f"  {'g3b_state_worst':<22} "
          f"{b['items']['g3b_state_worst']:>14,}")
    print(f"  {'TOTAL hard worst':<22} {b['total_hard_worst']:>14,}"
          f"   ({b['hours_hard_worst']:.2f} h)")
    print(f"\ncaps {MAX_CALLS:,} calls / {MAX_WALL_S:.0f} s "
          f"({MAX_WALL_S / 3600:.0f} h)")
    print(f"  vs projected  {b['cap_ratio_census_projected']:.3f}x")
    print(f"  vs hard worst {b['cap_ratio_hard_worst']:.3f}x"
          f"  -> cap can bind first: "
          f"{s['cap_contract']['may_fire_before_the_scan_ends']}")
    print(f"  on firing: {s['cap_contract']['on_firing']}")
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(s, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"\nartifact: {_ARTIFACT}")


if __name__ == "__main__":
    main()
