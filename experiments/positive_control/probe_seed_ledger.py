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
OBSERVED_PROBE_SCALARS = {
    "s3_pilot": 40_000_201,
    "w1_exploration": 40_000_211,
    "s3_smoke": 40_000_221,
    "s3_exploration": 40_000_231,
    "s4_campaign": 40_000_241,
    "s5_curved": 40_000_251,
    "s5_flat": 40_000_261,
}

#: Active fresh allocations, not yet observed when allocated. A
#: results commit MUST move the scalar to OBSERVED_PROBE_SCALARS in
#: the same change that adds the observed artifact.
FRESH_PROBE_SCALARS: dict[str, int] = {}

S3_PILOT_SEED = OBSERVED_PROBE_SCALARS["s3_pilot"]
W1_SEED = OBSERVED_PROBE_SCALARS["w1_exploration"]
S3_SMOKE_SEED = OBSERVED_PROBE_SCALARS["s3_smoke"]
S3_SEED = OBSERVED_PROBE_SCALARS["s3_exploration"]
S4_SEED = OBSERVED_PROBE_SCALARS["s4_campaign"]

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
