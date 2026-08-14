"""O4b sizing and the completion budget.

Two things are under test. The cost model, which must charge the G3
redesign honestly -- including the state calls the scan pays on
candidates it rejects. And the budget, which must refuse a reservation
BEFORE the calls are made, because a cap noticed afterwards is not a
cap.

Nothing here is a result. O4 has no verdict.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_g3_redesign as g3  # noqa: E402
import o4_sizing as base  # noqa: E402
import o4b_budget as bud  # noqa: E402
import o4b_sizing as sizing  # noqa: E402

# ------------------------------------------------- the frozen choices

def test_the_oracle_and_tau_are_inherited_not_re_frozen():
    assert sizing.TAU == base.TAU == 0.03
    assert sizing.N_G1 == base.N_G1 == 26_200_000
    floor = (base.V_HI - base.V_LO) / base.V_REF
    assert floor == pytest.approx(0.01999426, rel=1e-6)
    assert sizing.TAU > floor


def test_the_frozen_sizing_module_is_not_edited():
    """It is digest-pinned by the O4 manifest and by the executed
    snapshot the replay diagnostic checks against."""

    snapshot = json.loads(
        (_REPO / "docs" / "prereg"
         / "p14_o4_executed_freeze_manifest.json").read_text(
            encoding="utf-8"))
    import hashlib
    rel = "experiments/oracle/o4_sizing.py"
    got = hashlib.sha256((_REPO / rel).read_bytes()).hexdigest()
    assert got == snapshot["files"][rel]


# ----------------------------------------------------- the cost model

def test_the_cost_model_charges_rejected_candidates_too():
    """The state calls are per CANDIDATE, not per accepted cluster --
    which is what makes the total depend on availability."""

    assert sizing.STATE_CALLS_PER_CANDIDATE == 2
    assert sizing.PROBE_CALLS_PER_CLUSTER == 6
    assert sizing.g3b_calls(1_000, 10) == 2 * 1_000 + 6 * 10
    assert sizing.g3b_calls(100_000, 100_000) == 800_000


def test_the_new_g3_costs_more_than_the_old_one():
    b = sizing.budget()
    assert b["items"]["g3b_probes"] == 6 * g3.K_G3B
    old_g3 = 4 * base.G3_CLUSTERS
    projected_g3 = (b["items"]["g3b_state_census_projected"]
                    + b["items"]["g3b_probes"])
    assert projected_g3 == 2 * old_g3


def test_the_totals_and_cap_ratios():
    b = sizing.budget()
    assert b["total_census_projected"] == 55_351_840
    assert b["total_hard_worst"] == 107_551_840
    assert b["cap_ratio_census_projected"] == pytest.approx(1.445,
                                                            abs=5e-4)
    assert b["cap_ratio_hard_worst"] == pytest.approx(0.744, abs=5e-4)
    assert b["hours_census_projected"] == pytest.approx(11.808,
                                                        abs=5e-3)


def test_the_hard_worst_case_is_declared_an_upper_bound():
    b = sizing.budget()
    assert b["items"]["g3b_state_worst"] == 2 * sizing.scan_cap()
    assert "nothing committed bounds that fraction" in (
        b["worst_case_is_an_upper_bound_because"])


def test_the_scan_cap_is_delegated_and_also_recorded():
    assert sizing.scan_cap() == g3.scan_cap() == base.N_G1
    assert sizing.SCAN_CAP_REALIZED == 26_200_000


def test_the_g3a_cost_is_now_measured_end_to_end():
    """The reckoning came out 72 low. The preflight is executed under
    the meter instead, so the number is what it spends."""

    detail = sizing.g3a_calls()
    assert detail["basis"] == (
        "end-to-end measurement of run_preflight under the meter")
    assert "single authority" in detail["measures"]
    # named rather than counted, so a condition that gets renamed or
    # dropped fails here instead of being replaced by a new one
    assert detail["preflight_conditions"] == {
        "tri_state_rows": True,
        "recovered_dpsi_bit_identical": True,
        "families_covered": True,
        "ab_realized_margin": True,
        "only_expected_unreachable": True,
        "every_case_built": True,
        "row_d_no_solver_call": True,
        "solver_determinism": True,
    }
    # ulp distances are HOST-DEPENDENT: the solver reaches libm, so
    # `t_min` and `err` differ in their last bits and the distance to a
    # satisfying placement moves with them (CI 8,111 vs 8,120 here).
    # The call counts do not move, which is why they can be pinned and
    # these cannot.
    assert detail["ulp_distance_total"] >= 8_000
    assert detail["search_comparisons_total"] < 1_000
    assert "is the work actually done" in (
        detail["distance_is_not_cost"])
    assert set(detail["host_dependent"]) == {"ulp_distance_total",
                                             "search_comparisons_total"}
    assert "never be pinned as constants" in detail["why_host_dependent"]
    assert detail["calls"] == detail["completed"] == 6_448
    assert detail["cases_resolved"] == 70
    assert detail["cases_unreachable"] == 7
    assert detail["construction_unavailable_rows"] == 2
    assert set(detail["families_covered"]) == set(g3.FAMILIES)
    assert "no calls by construction" in detail["row_d_costs_nothing"]


# ------------------------------------------------- the census lineage

def test_the_projection_reads_its_input_from_the_committed_census():
    """The earlier draft quoted a rate from an uncommitted run log
    that no total used. Provenance or nothing."""

    row = sizing.census_eligibility()
    assert row["source"] == "docs/prereg/p14_o4_replay_diagnostic.json"
    assert (_REPO / row["source"]).exists()
    assert row["eta"] == g3.ETA
    assert row["eligible_clusters"] == row["clusters_probed"] == 100_000
    assert row["lower_probe_unreachable_among_eligible"] == 0
    assert "NOT a probabilistic expectation" in row["is"]


def test_it_is_called_projected_not_expected():
    b = sizing.budget()
    assert "total_census_projected" in b
    assert "total_expected" not in b
    source = (_REPO / "experiments" / "oracle"
              / "o4b_sizing.py").read_text(encoding="utf-8")
    assert "CENSUS_PREFIX_POINTS" not in source


# --------------------------------------------------- the outcome table

def test_invalid_is_not_g3a_alone():
    """A G3b probe that comes back undecided or boolean-mismatched is a
    contract failure too -- the redesign requires mismatch 0."""

    outcomes = {o["outcome"] for o in sizing.OUTCOMES}
    assert outcomes == {"INVALID", "INCONCLUSIVE",
                        "no scientific verdict",
                        "no stage outcome"}
    invalid = [o["condition"] for o in sizing.OUTCOMES
               if o["outcome"] == "INVALID"]
    assert any("G3a" in c for c in invalid)
    assert any("fully-testable G3b probe" in c for c in invalid)
    assert "INVALID belongs to contract failures wherever they occur" in (
        sizing.cap_contract()["not_invalid"])


def test_shortfalls_and_caps_are_inconclusive():
    inconclusive = [o["condition"] for o in sizing.OUTCOMES
                    if o["outcome"] == "INCONCLUSIVE"]
    assert any("fewer than K_G3B" in c for c in inconclusive)
    assert any("cap binds" in c for c in inconclusive)


def test_one_unavailable_candidate_is_not_a_stage_outcome():
    """It is tallied and the scan carries on. Only failing to reach
    K_G3B by the end of the scan is INCONCLUSIVE."""

    row = next(o for o in sizing.OUTCOMES
               if "a single candidate" in o["condition"])
    assert row["outcome"] == "no stage outcome"
    assert "scan CONTINUES" in row["when"]


def test_the_cap_may_bind_and_says_so():
    contract = sizing.cap_contract()
    assert contract["may_fire_before_the_scan_ends"] is True
    assert contract["on_firing"].endswith("INCONCLUSIVE")
    assert contract["auto_raise"].startswith("forbidden")


# ------------------------------------------------------- the budget

def test_a_reservation_is_refused_before_the_calls_are_made():
    """The O4 runner incremented and then compared, so a batch could
    carry the total past the cap. Here it cannot."""

    b = bud.Budget(max_calls=10, max_wall_s=1e9)
    b.reserve(6)
    assert b.reserved == 6
    with pytest.raises(bud.CapReached) as caught:
        b.reserve(5)
    assert b.reserved == 6                   # refused, not charged
    assert caught.value.reason == "max-calls"
    assert caught.value.wanted == 5


def test_the_boundary_is_inclusive():
    b = bud.Budget(max_calls=10, max_wall_s=1e9)
    b.reserve(10)
    assert b.reserved == 10
    with pytest.raises(bud.CapReached):
        b.reserve(1)


def test_reserved_and_completed_are_kept_apart():
    """A batch that dies part way through completes fewer calls than
    it reserved, and the incident must be able to say which."""

    b = bud.Budget(max_calls=100, max_wall_s=1e9)
    b.reserve(8)
    b.complete(3)                     # the solver raised on the fourth
    state = b.state()
    assert (state["reserved"], state["completed"]) == (8, 3)
    assert state["unfinished"] == 5
    assert "never under-charges" in state["cap_is_charged_against"]


def test_the_wall_clock_also_refuses_before_spending():
    ticks = iter([0.0, 0.0, 5.0, 5.0, 5.0])
    b = bud.Budget(max_calls=10 ** 9, max_wall_s=4.0,
                   clock=lambda: next(ticks))
    b.reserve(1)
    with pytest.raises(bud.CapReached) as caught:
        b.reserve(1)
    assert caught.value.reason == "max-wall"
    assert b.reserved == 1


def test_the_wall_cap_is_a_start_rule_with_no_bound_in_seconds():
    """It does not stop a batch already in flight. The overrun is one
    batch in CALLS, but the 768 us/pair cost is a host-dependent
    average and not a per-call maximum, so there is no rigorous bound
    in seconds and the module must not claim one."""

    b = bud.Budget(max_calls=10 ** 9, max_wall_s=100.0,
                   max_batch_calls=8_192)
    says = b.state()["wall_cap_is"]
    assert "a start rule" in says
    assert "NO rigorous bound in seconds" in says
    assert "host-dependent average" in says
    assert not hasattr(b, "max_overrun_s")
    assert b.nominal_overrun_s == pytest.approx(8_192 * 0.768e-3)
    assert "projection" in b.state()["nominal_overrun_is"]
    with pytest.raises(ValueError, match="exceeds max_batch_calls"):
        b.reserve(8_193)


def test_the_call_cap_is_never_crossed_at_all():
    """Unlike the wall cap, this one has no overshoot: a batch that
    does not fit is refused outright."""

    b = bud.Budget(max_calls=10, max_wall_s=1e9, max_batch_calls=8)
    b.reserve(8)
    with pytest.raises(bud.CapReached):
        b.reserve(3)
    assert b.reserved == 8 <= b.max_calls
    assert "never crossed" in b.state()["call_cap_overshoot"]


def test_the_actual_overrun_is_recorded_not_projected():
    now = [0.0]
    b = bud.Budget(max_calls=10, max_wall_s=4.0,
                   clock=lambda: now[0])
    now[0] = 5.0
    assert b.wall_overrun_s == pytest.approx(1.0)
    assert b.state()["wall_overrun_s"] == pytest.approx(1.0)


def test_completed_may_never_exceed_reserved():
    """A completion without a reservation is a call the budget never
    charged for -- the very thing pre-charging exists to prevent."""

    b = bud.Budget(max_calls=100, max_wall_s=1e9)
    with pytest.raises(ValueError, match="must have been charged"):
        b.complete(3)                       # nothing reserved at all
    b.reserve(2)
    b.complete(2)                           # exactly the reservation
    assert b.completed == 2 and b.unfinished == 0
    with pytest.raises(ValueError, match="must have been charged"):
        b.complete(1)                       # one past the reservation
    assert b.completed == 2
    b.reserve(4)
    b.complete(1)                           # partial, allowed
    assert (b.reserved, b.completed, b.unfinished) == (6, 3, 3)


def test_the_counters_are_read_only_too():
    b = bud.Budget(max_calls=10, max_wall_s=1e9)
    for name in ("reserved", "completed", "unfinished"):
        with pytest.raises(AttributeError):
            setattr(b, name, 99)


def test_the_clock_and_counter_start_before_g3a():
    b = bud.Budget(max_calls=10, max_wall_s=1e9)
    assert b.stage == "g3a"
    b.reserve(3)
    assert b.per_stage == {"g3a": 3}
    b.enter("g1")
    b.reserve(2)
    assert b.per_stage == {"g3a": 3, "g1": 2}
    assert b.state()["clock_started_at"] == "before G3a"


def test_the_limits_are_read_only_not_merely_unraised():
    b = bud.Budget(max_calls=10, max_wall_s=1e9)
    for name in ("max_calls", "max_wall_s", "max_batch_calls"):
        with pytest.raises(AttributeError):
            setattr(b, name, 10 ** 12)
    for forbidden in ("raise_cap", "extend", "reset", "set_max_calls"):
        assert not hasattr(b, forbidden)


def test_a_cap_breach_carries_the_state_the_incident_must_report():
    b = bud.Budget(max_calls=3, max_wall_s=1e9)
    b.reserve(3)
    b.complete(2)
    with pytest.raises(bud.CapReached) as caught:
        b.reserve(1)
    exc = caught.value
    assert (exc.reserved, exc.completed, exc.wanted) == (3, 2, 1)
    assert "refused a reservation" in str(exc)
    assert exc.reserved == b.state()["reserved"]
    assert exc.completed == b.state()["completed"]


# ------------------------------------------- the projection's honesty

def test_the_projection_declares_its_optimistic_assumption():
    """The census never measured the realized-margin or nudge
    conditions, so `fully-testable = 1` is an assumption, not a rate."""

    says = sizing.budget()["census_projected_assumes"]
    assert "fully-testable rate = 1" in says
    assert "OPTIMISTIC projection" in says
    row = sizing.census_eligibility()
    assert "the realized-margin check" in row["not_measured_here"]
    assert "the nudge search" in row["not_measured_here"]


def test_no_deleted_lineage_is_quoted_anywhere():
    source = (_REPO / "experiments" / "oracle"
              / "o4b_sizing.py").read_text(encoding="utf-8")
    for gone in ("2.4%", "99,391", "4,149,248", "positive_rate"):
        assert gone not in source


def test_the_committed_artifact_matches_the_module():
    """A stale artifact is a freeze that says something the code does
    not."""

    artifact = json.loads(
        (_REPO / "docs" / "prereg"
         / "p14_o4b_sizing.json").read_text(encoding="utf-8"))
    fresh = json.loads(json.dumps(sizing.summary()))

    # the host-dependent diagnostics are allowed to differ between the
    # machine that published and the machine that checks; everything
    # the freeze actually rests on is not
    volatile = set(artifact["budget"]["g3a_detail"]["host_dependent"])
    for side in (artifact, fresh):
        for key in volatile:
            side["budget"]["g3a_detail"].pop(key, None)
    assert artifact == fresh
