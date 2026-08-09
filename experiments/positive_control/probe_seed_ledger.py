"""Aggregated scalar-seed ledger for post-stage exploration probes.

ONE SOURCE for "every scalar seed the program has spent": the frozen
chains keep their own seed constants, and this module UNIONS them
instead of hand-copying values -- the hand-copied set in the first
S3/W1 cut missed 20260831/32 and 20260841/42 (PR review). A new probe
registers its scalar in PROBE_SEEDS and asserts freshness through
`assert_probe_seed_fresh`, which excludes only the probe's own seed
from the spent set.

`seed_windows` is digest-pinned (frozen), which is why the aggregation
lives here and not there.
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

S3_SEED = 40_000_201  # S3 Schwarzschild paired exploration
W1_SEED = 40_000_211  # W1 channel-decomposition exploration

PROBE_SEEDS = {"s3_exploration": S3_SEED, "w1_exploration": W1_SEED}

SPENT_RANGES = P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,)


def spent_scalars(excluding: str) -> frozenset[int]:
    """Every spent scalar seed except the named probe's own."""

    spent = (set(p3c.BURNED_SEEDS) | set(p3c.CAMPAIGN_SEEDS.values())
             | set(p3e.BURNED_SEEDS) | set(p3e.SEEDS.values())
             | set(p2.BURNED_SEEDS) | set(p2.CAMPAIGN_SEEDS.values())
             | {p2.MC_SEED, p2.DESIGN_MC_SEED}
             | set(p14_prereg.CAMPAIGN_SEEDS.values())
             | {s1.BENCH_SEED}
             | set(PROBE_SEEDS.values()))
    spent.discard(PROBE_SEEDS[excluding])
    return frozenset(spent)


def assert_probe_seed_fresh(name: str) -> None:
    """Abort unless the named probe's seed is fresh against the full
    aggregated ledger (scalars and ranges)."""

    assert_point_seeds_fresh({name: PROBE_SEEDS[name]},
                             spent_scalars(name), SPENT_RANGES, name)
