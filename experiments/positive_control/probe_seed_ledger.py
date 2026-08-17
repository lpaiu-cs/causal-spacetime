"""Aggregated scalar-seed ledger for post-stage exploration probes.

ONE SOURCE for "every scalar seed the program has spent": the frozen
chains keep their own seed constants, and this module UNIONS them
instead of hand-copying values -- the hand-copied set in the first
S3/W1 cut missed 20260831/32 and 20260841/42, and the first union
missed S1's `BENCH_SEED + 1` reference-bench stream (PR review).

Fresh vs replay (PR review R2): a seed is SPENT the moment results
were observed on it, and observation is a ledger fact, not something
an API may silently wave through. The first cut's
`spent_scalars(excluding=...)` excluded the asking probe's own seed,
which certified the S3 pilot seed as "fresh" for the official rerun.
The two operations are now distinct:

- `assert_fresh_scalar(name)`: a NEW allocation -- the seed must be
  outside every observed scalar, every other active allocation, and
  every spent range. Used exactly once per experiment.
- `replay_scalar(name)`: a deterministic REPRODUCTION of an already
  observed run -- returns the recorded observed seed, and the caller
  must label its output as a replay, never as a new observation.
"""

from __future__ import annotations

import p14_prereg
import p14_probe_p2 as p2
import p14_probe_p3c as p3c
import p14_probe_p3e as p3e
import s1_schwarzschild_cost as s1
from seed_windows import (
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_point_seeds_fresh,
)

#: Probe scalars whose results HAVE BEEN OBSERVED -- spent forever.
#: s3_pilot: the fixed-N dirty-tree run preserved as pilot only.
#: w1_exploration: the approved W1 artifact (its official regeneration
#: was a deterministic replay of the same stream, not a new draw).
#: s3_smoke: first allocated as the official S3 seed, then
#: smoke-observed during review-R2 validation (twice, 3 readings at
#: E_N = 40) -- demoted to the dedicated smoke/validation stream by
#: the review-R3 principle that observing ANY output spends the seed.
#: s3_exploration: the official S3 artifact was observed 2026-08-10
#: (moved from FRESH in the same commit that added the artifact, per
#: the results-commit obligation below).
#: s4_campaign: the S4 preregistered confirmation observed 2026-08-10
#: (verdict CONFIRMED; docs/prereg/p14_s4_results.json).
#: s5_curved / s5_flat: the S5 C2-unpaired campaign arms observed
#: 2026-08-11 (outcome DETECTED; docs/prereg/p14_s5_results.json).
#: oracle_mc_diagnostic: the volume-oracle Monte Carlo cross-check
#: stream, first observed 2026-08-11 (PR-O2 development). DIAGNOSTIC
#: ONLY -- no statistical claim ever attaches to it (the oracle is
#: deterministic; MC disjointness raises an investigation flag, never
#: a verdict), but the observing-spends principle still applies: this
#: stream is permanently spent for any claim-bearing use.
#: o4_smoke: the O4 validation stream. Spent from the moment it was
#: allocated, because the contract tests run smoke assemblies on it
#: and observing ANY output spends the seed (the S3 review-R3
#: principle); it exists precisely so smoke may never draw a campaign
#: seed. Its two consecutive successors are reserved with it: the
#: smoke path derives its G2 stream as SEED+1, and `o4_smoke_g3`
#: (SEED+2) stays retired although the G3 seed-inheritance fix left it
#: unread -- a seed once listed as spent is never un-spent, since the
#: safe direction of that error is to retire one value too many.
#: o4_aborted_g1 / o4_aborted_g2: the O4 campaign strata, drawn on
#: 2026-08-12 from the freeze merge `1eb9461`. The run completed G1's
#: 26.2M points and G2's call, entered G3, and stopped on the frozen
#: fail-closed path (`causal_relation` undecided at a stress point), so
#: it published NO result and NO gate has a status
#: (docs/prereg/p14_o4_incident.json). The names say `aborted` rather
#: than `campaign` precisely so nothing downstream can read these as
#: the provenance of a verdict. They are spent all the same: the points
#: were observed, and the streams are retired whether or not a verdict
#: came of them.
#: o4b_g1_audit / o4b_g2_leakage: the O4b campaign strata, drawn on
#: 2026-08-15 from the clean checkout of freeze SHA `715865a` (the
#: freeze branch head, merged in PR #75). The run passed G3a's 6,537
#: metered wrapper-contract calls, claimed `refs/o4b/reservation`
#: (`46acee34`), and ran G1 (26,200,000 points, concordant) and G2
#: (1,072,696 points, 0 leaks, concordant) to their frozen sample
#: sizes. It then stopped at the pre-publication provenance check --
#: NOT on a scientific fail-closed path and NOT on a cap, but on a
#: pure implementation defect: `main`'s claim callback stored the
#: reservation object in an external dict yet returned `None`, so the
#: exit re-read received `None` and raised before any result was
#: published (docs/prereg/p14_o4b_incident.json). No verdict; the gate
#: statistics are preserved there and in the checkpoint as non-verdict
#: partials. The names stay FUNCTIONAL rather than `aborted` because
#: the defect is in publication wiring, not the science -- both gates
#: completed WITH a status, and whether the preserved statistics become
#: the provenance of a verdict is decided by a separate recovery audit,
#: not prejudged here. They are spent regardless: the points were
#: observed, and `refs/o4b/reservation` is RETAINED, not deleted.
#: o5_amb_pilot: the O5 ambiguity pilot's fixed-n indicator stream,
#: drawn 2026-08-17 from the clean exact checkout of freeze SHA
#: `e9e7e15` (the freeze branch head, merged in PR #83). The run
#: scanned all 10,736,965 events in one pass (21,473,930 predicate
#: calls -- exactly 2n -- after G3a's 6,537 frozen wrapper-preflight
#: calls; k_ambiguous = 0), published FEASIBLE under the frozen
#: general-in-k rule (docs/prereg/p14_o5_amb_pilot.json), and
#: re-verified `refs/o5pilot/reservation` (`f522363`) at exit; that
#: ref is RETAINED, not deleted. One stream, one run, never extended:
#: a re-run needs a NEW freeze with a NEW scalar, and replaying this
#: one is a reproduction, never an observation (`replay_scalar`).
#: o5_campaign: the O5 Poisson-count campaign's sprinkling stream,
#: drawn 2026-08-17 from the clean exact checkout of freeze SHA
#: `c6eb85e` (the freeze branch head, merged in PR #87) under the
#: integrated execution approval. The run drew N = 26,831,117
#: proposals (Poisson at A*SCALE), scanned them all (29,335,028
#: predicate calls after G3a's 6,537; K_certain = 53,285, zero
#: ambiguous), published CONCORDANT under the frozen general rule
#: (docs/prereg/p14_o5_count.json: D = [-0.855, +0.677] contained in
#: [-1.4194, +1.4194]), and re-verified `refs/o5/reservation`
#: (`b5d9d61`) at exit; that ref is RETAINED, not deleted. One
#: stream, one run, never extended: a re-run needs a NEW freeze with
#: a NEW scalar, and replaying this one is a reproduction, never an
#: observation (`replay_scalar`).
OBSERVED_PROBE_SCALARS = {
    "s3_pilot": 40_000_201,
    "w1_exploration": 40_000_211,
    "s3_smoke": 40_000_221,
    "s3_exploration": 40_000_231,
    "s4_campaign": 40_000_241,
    "s5_curved": 40_000_251,
    "s5_flat": 40_000_261,
    "oracle_mc_diagnostic": 40_000_271,
    "o4_smoke": 40_000_311,
    "o4_smoke_g2": 40_000_312,
    "o4_smoke_g3": 40_000_313,
    "o4_aborted_g1": 40_000_281,
    "o4_aborted_g2": 40_000_291,
    "o4b_g1_audit": 40_000_401,
    "o4b_g2_leakage": 40_000_411,
    "o5_amb_pilot": 40_000_441,
    "o5_campaign": 40_000_451,
}

#: Active fresh allocations, not yet observed when allocated. A
#: results commit MUST move the scalar to OBSERVED_PROBE_SCALARS in
#: the same change that adds the observed artifact -- and so must an
#: ABORT record, since spending does not depend on the run succeeding.
#:
#: Empty: the O4 strata `g1_audit`/`g2_leakage` were drawn on
#: 2026-08-12 and are now retired as `o4_aborted_g1`/`o4_aborted_g2`.
#: A re-run of O4 needs a NEW freeze with NEW scalars; the retired
#: streams may not be re-entered, and replaying them is a reproduction,
#: never an observation (`replay_scalar`).
#:
#: There is deliberately no G3 scalar. G3's unit is the boundary-stress
#: cluster drawn from `G1 measure | L_S1 > 0`, and it inherits G1's own
#: accepted points; the v1 O4 draft allocated `g3_wrapper = 40_000_301`
#: and then never read it, which would have recorded a spent seed as G3
#: provenance for samples G1 produced (O4 review R1). The value is
#: withdrawn unspent -- it was never drawn from, so it is not moved to
#: OBSERVED; `40_000_301` simply returns to the unallocated pool, and
#: the abort does not change that (nothing ever drew from it).
#: O4b: the re-run under the redesigned G3. NEW scalars, because the
#: O4 streams were drawn from and are retired whether or not a verdict
#: came of them. `40_000_301` stays where the O4 review left it --
#: withdrawn unspent, back in the unallocated pool -- and is
#: deliberately NOT reused here: reaching for it would make the
#: withdrawal look like a deferral.
#:
#: The O4b strata `o4b_g1_audit`/`o4b_g2_leakage` were drawn on
#: 2026-08-15 and are retired in OBSERVED under those same names; the
#: recovered result is a REPRODUCTION of that observed run, and a
#: re-run of O4b needs a NEW freeze with NEW scalars.
#:
#: `o5_amb_pilot` (allocated 2026-08-16 under the delegated pilot
#: freeze) was drawn and is retired in OBSERVED under the same name,
#: in the same commit that added the observed FEASIBLE artifact --
#: the results-commit obligation above, honoured for an outcome as it
#: would have been for an abort.
#:
#: Empty again since 2026-08-17: `o5_campaign` (allocated 2026-08-17
#: under the campaign freeze, docs/prereg/p14_o5_count.md) was drawn
#: under the integrated execution approval and is retired in OBSERVED
#: under the same name, in the same commit that added the observed
#: CONCORDANT artifact -- the results-commit obligation above. The
#: next allocation belongs to whatever stage the PI rules next (S6
#: mass ladder is the stated candidate); nothing here pre-allocates
#: it.
FRESH_PROBE_SCALARS: dict[str, int] = {}

S3_PILOT_SEED = OBSERVED_PROBE_SCALARS["s3_pilot"]
W1_SEED = OBSERVED_PROBE_SCALARS["w1_exploration"]
S3_SMOKE_SEED = OBSERVED_PROBE_SCALARS["s3_smoke"]
S3_SEED = OBSERVED_PROBE_SCALARS["s3_exploration"]
S4_SEED = OBSERVED_PROBE_SCALARS["s4_campaign"]
O4_SMOKE_SEED = OBSERVED_PROBE_SCALARS["o4_smoke"]

SPENT_RANGES = P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,)


def spent_scalars() -> frozenset[int]:
    """Every scalar seed with observed results, program-wide."""

    return frozenset(
        set(p3c.BURNED_SEEDS) | set(p3c.CAMPAIGN_SEEDS.values())
        | set(p3e.BURNED_SEEDS) | set(p3e.SEEDS.values())
        | set(p2.BURNED_SEEDS) | set(p2.CAMPAIGN_SEEDS.values())
        | {p2.MC_SEED, p2.DESIGN_MC_SEED}
        | set(p14_prereg.CAMPAIGN_SEEDS.values())
        | {s1.BENCH_SEED, s1.BENCH_SEED + 1}
        | set(OBSERVED_PROBE_SCALARS.values()))


def assert_fresh_scalar(name: str) -> int:
    """Abort unless `name`'s allocation is fresh against the observed
    scalars, the other active allocations, and the spent ranges.
    Returns the seed."""

    if name not in FRESH_PROBE_SCALARS:
        retired = sorted(k for k in OBSERVED_PROBE_SCALARS
                         if k.endswith(name.split("_")[0])
                         or name.split("_")[0] in k)
        raise KeyError(
            f"{name!r} has no fresh allocation. A retired stream is "
            f"never re-entered: allocate a NEW scalar under a new "
            f"freeze, or use replay_scalar() for a reproduction that "
            f"is labelled as one. Related retired entries: {retired}")
    seed = FRESH_PROBE_SCALARS[name]
    others = {n: s for n, s in FRESH_PROBE_SCALARS.items() if n != name}
    assert_point_seeds_fresh(
        {name: seed, **others},
        spent_scalars(), SPENT_RANGES, name)
    return seed


def replay_scalar(name: str) -> int:
    """The recorded observed seed for a deterministic reproduction.
    The caller must label the output as a replay of the observed run."""

    return OBSERVED_PROBE_SCALARS[name]
