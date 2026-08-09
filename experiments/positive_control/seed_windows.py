"""Seed-freshness checks, extracted from P12 Stage B as pure functions.

`p14_weyl_curvature.md` section 8.2 told the eventual preregistration
to inherit P12's disjointness helper rather than write a third
variant; the P14 preregistration review (R4) added that the P12
function is hardcoded to P12's module constants, so the logic lives
here as parameterized PURE functions and `p12_stage_b` keeps a
wrapper with its original name and behavior.

Two kinds of streams exist and each needs its own check (prereg
v0.4): integer scalar seeds (this module -- window/range blocks and
point seeds) and `SeedSequence` entropy/spawn-key paths (structural
key uniqueness, asserted where the layout is defined, because their
identity is the `(entropy, spawn_key)` pair and not an integer that
ranges can describe).

Failure is `SystemExit`, matching the P12 original: a seed collision
is a protocol violation, not a condition to handle.
"""

from __future__ import annotations

Range = tuple[str, int, int]


def assert_windows_disjoint_and_fresh(
        blocks: list[Range],
        spent_ranges: tuple[Range, ...],
        window_floor: int,
        window_ceil: int,
        stage: str) -> list[Range]:
    """Pairwise disjointness, window containment, and freshness for
    (label, lo, hi) INCLUSIVE integer blocks. Aborts on violation."""

    for label, lo, hi in blocks:
        if not (window_floor <= lo and hi <= window_ceil):
            raise SystemExit(
                f"seed block {label} [{lo}..{hi}] leaves {stage}'s "
                f"window [{window_floor}..{window_ceil}]")
        for name, slo, shi in spent_ranges:
            if not (hi < slo or shi < lo):
                raise SystemExit(
                    f"seed block {label} [{lo}..{hi}] overlaps the "
                    f"already-spent {name} [{slo}..{shi}]")
    for i, (la, loa, hia) in enumerate(blocks):
        for lb, lob, hib in blocks[i + 1:]:
            if not (hia < lob or hib < loa):
                raise SystemExit(
                    f"seed blocks overlap: {la} [{loa}..{hia}] vs "
                    f"{lb} [{lob}..{hib}]")
    return blocks


def assert_point_seeds_fresh(
        seeds: dict[str, int],
        spent_scalars: frozenset[int],
        spent_ranges: tuple[Range, ...],
        stage: str) -> None:
    """Freshness for INDIVIDUAL scalar seeds: pairwise distinct, not
    in the spent scalar set, and inside no spent range. Aborts."""

    values = list(seeds.values())
    if len(set(values)) != len(values):
        raise SystemExit(f"{stage} seeds repeat: {sorted(values)}")
    for label, s in seeds.items():
        if s in spent_scalars:
            raise SystemExit(
                f"{stage} seed {label}={s} was already spent")
        for name, slo, shi in spent_ranges:
            if slo <= s <= shi:
                raise SystemExit(
                    f"{stage} seed {label}={s} sits inside the "
                    f"already-spent {name} [{slo}..{shi}]")
