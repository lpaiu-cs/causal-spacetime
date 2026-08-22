"""Paper A <-> artifact integration contracts for Sections 4 through 6.7.

Until now the manuscript's own numbers were machine-checked in exactly
two places: `test_paper_a_count_integration.py` re-derives the Section
6.8 per-rung table, and `test_o4b_downstream.py` pins the O4b audit's
figures verbatim. `test_paper_a_manifest.py` recomputes FILE digests
only -- it never compares a manuscript number to an artifact -- so
editing a value like Section 4.2's volume-estimator RMSE passed every
gate. This file closes that window over the rest of the results.

The shape follows the Section 6.8 contract. Every claim below states
WHERE it is printed, WHICH artifact produces it, and HOW to get from
the artifact to the printed characters; the expected value appears
nowhere in this file. `derive()` returns the formatted string(s), and
the assertion is that the manuscript sentence -- with those strings
substituted in -- occurs in the section that is supposed to carry it.
Pinning the surrounding sentence, not the bare number, is what makes
the check sensitive to a value that silently moves: a whole-document
substring search would still pass if 4.2's two RMSE figures swapped.

Numbers that no artifact can produce (a mathematical constant, a row
selector, a construction setting) are NOT skipped: each is listed in
`ACCEPTED_EXCLUSIONS` bound to the exact PHRASE that justifies it, in
the section it occurs in, with a one-line reason.
`test_every_number_is_derived_or_explicitly_excluded` strikes out
every sentence a claim matched and every excused phrase, then requires
nothing numeric to be left. An excuse therefore covers the sentence it
was written about and nothing else -- not the same literal elsewhere,
and not a rewrite of its own sentence into a new result. That guard is
the actual promise; the claim list alone could always be out-run by a
sentence added later.
"""

from __future__ import annotations

import csv
import functools
import importlib.util
import json
import math
import re
import statistics
from collections.abc import Callable
from pathlib import Path
from typing import NamedTuple

import pytest

REPO = Path(__file__).resolve().parents[1]
_PREREG = REPO / "docs" / "prereg"
_PAPER = REPO / "docs" / "paper" / "paper_a"
_FIG_DATA = _PAPER / "figures" / "data"

MANUSCRIPT = _PAPER / "manuscript.md"
CLAIM_BOUNDARY = _PAPER / "claim_boundary.md"


# --------------------------------------------------------- artifacts

def _csv(name: str) -> list[dict[str, str]]:
    """One of the 19 committed legacy summary tables."""

    with (_FIG_DATA / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _json(name: str) -> dict:
    """One of the p14_* preregistered artifacts."""

    return json.loads((_PREREG / name).read_text(encoding="utf-8"))


def _text(rel: str) -> str:
    return (REPO / rel).read_text(encoding="utf-8")


def _rows(name: str, **eq: float) -> list[dict[str, str]]:
    """Rows whose numeric columns equal the given values."""

    return [r for r in _csv(name)
            if all(float(r[k]) == v for k, v in eq.items())]


def _one(name: str, **eq: float) -> dict[str, str]:
    matched = _rows(name, **eq)
    assert len(matched) == 1, (name, eq, len(matched))
    return matched[0]


def _col(rows: list[dict[str, str]], key: str) -> list[float]:
    return [float(r[key]) for r in rows]


def _num(value: float) -> str:
    """Shortest exact decimal, the way the manuscript prints a count or
    an exact ratio: 496.0 -> '496', 14.75 -> '14.75', 2.25 -> '2.25'."""

    return str(int(value)) if float(value).is_integer() else repr(float(value))


@functools.cache
def _volume_ratio(wT: float) -> float:
    """`p14_checks/p14_interval_volume_constant_a.py` is a committed
    design check rather than an importable package module; it is
    digest-locked in the artifact manifest like every other artifact
    read here, and it is pure numpy quadrature (no RNG, ~50 ms)."""

    return _interval_volume_check().volume_ratio(wT)


@functools.cache
def _interval_volume_check():
    path = _PREREG / "p14_checks" / "p14_interval_volume_constant_a.py"
    spec = importlib.util.spec_from_file_location("_p14_interval_volume", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------- manuscript sections

_HEADING = re.compile(r"^#{2,3} (\d+(?:\.\d+)?)[. ]")


def _sections() -> dict[str, str]:
    """The manuscript split by its own numbered headings, so a claim is
    checked against the section that is supposed to carry it. Heading
    lines are dropped: '### 4.1 R0' would otherwise donate a '4.1' to
    the numeric scan below."""

    out: dict[str, list[str]] = {}
    current: list[str] | None = None
    for line in MANUSCRIPT.read_text(encoding="utf-8").splitlines():
        heading = _HEADING.match(line)
        if heading:
            current = out.setdefault(heading.group(1), [])
            continue
        if line.startswith("#"):
            current = None
            continue
        if current is not None:
            current.append(line)
    return {k: " ".join(" ".join(v).split()) for k, v in out.items()}


SECTIONS = _sections()


def _section(number: str) -> str:
    """Section text, whitespace-flattened so markdown line wrapping
    cannot hide a needle. A parent number carries its subsections."""

    parts = [v for k, v in SECTIONS.items()
             if k == number or k.startswith(number + ".")]
    assert parts, number
    return " ".join(parts)


# ------------------------------------------------------- claim table

class Claim(NamedTuple):
    """One quantitative statement in the manuscript.

    `context` is the sentence as printed, with `{}` where each derived
    value belongs; `derive` returns those values already formatted the
    way the manuscript formats them.
    """

    section: str
    source: str
    context: str
    derive: Callable[[], tuple[str, ...]]


CLAIMS: list[Claim] = []


def claim(section: str, source: str, context: str) -> Callable:
    def wrap(fn: Callable[[], tuple[str, ...]]) -> Claim:
        made = Claim(section, source, context, fn)
        CLAIMS.append(made)
        return made
    return wrap


# ============================================ 4.1 R0 -- order alone

_DIM_N = (300.0, 600.0, 1200.0, 2400.0)


def _dim_track(dim: float) -> str:
    return " -> ".join(
        "{:.2f}".format(float(_one("dimension_reconstruction_summary.csv",
                                   spacetime_dim=dim,
                                   N=n)["mean_estimated_dim"]))
        for n in _DIM_N)


@claim("4.1", "exp10 dimension_reconstruction_summary.csv",
       "for true dimension 2 the estimate moves {} at N = 300, 600, 1200, 2400")
def _dim2() -> tuple[str, ...]:
    return (_dim_track(2.0),)


@claim("4.1", "exp10 dimension_reconstruction_summary.csv",
       "for dimension 3, {};")
def _dim3() -> tuple[str, ...]:
    return (_dim_track(3.0),)


@claim("4.1", "exp10 dimension_reconstruction_summary.csv",
       "for dimension 4, {}.")
def _dim4() -> tuple[str, ...]:
    return (_dim_track(4.0),)


def _chain_rows() -> list[dict[str, str]]:
    return _rows("longest_chain_calibration_summary.csv", rho=300.0)


@claim("4.1", "exp09 longest_chain_calibration_summary.csv",
       "at rho = 300 the normalized length L/(sqrt(rho) tau) is {} with "
       "endpoints included and {} with endpoints removed, versus the "
       "asymptotic sqrt(2) ~ {}")
def _chain_normalized() -> tuple[str, ...]:
    rows = _chain_rows()

    def normalized(key: str) -> str:
        values = [float(r[key]) / (math.sqrt(float(r["rho"])) * float(r["true_tau"]))
                  for r in rows]
        return f"{statistics.fmean(values):.3f}"

    return (normalized("chain_length_including_endpoints"),
            normalized("effective_chain_length_minus_endpoints"),
            f"{math.sqrt(2):.3f}")


@claim("4.1", "exp09 longest_chain_calibration_summary.csv",
       "its mean error there is {};")
def _chain_bias() -> tuple[str, ...]:
    return (f'{statistics.fmean(_col(_chain_rows(), "tau_chain_error")):.3f}',)


# ============================== 4.2 R1 -- order + measure: proper time

def _pair_rmse(key: str) -> tuple[str, str]:
    return tuple(  # type: ignore[return-value]
        "{:.3f}".format(float(_one("timelike_pair_reconstruction_summary.csv",
                                   N=n)[key]))
        for n in (300.0, 2400.0))


@claim("4.2", "exp07 timelike_pair_reconstruction_summary.csv",
       "The volume-estimator relative RMSE falls from {} at N = 300 to {} at "
       "N = 2400.")
def _volume_rmse() -> tuple[str, ...]:
    return _pair_rmse("tau_volume_relative_rmse")


@claim("4.2", "exp07 timelike_pair_reconstruction_summary.csv",
       "below the chain estimator's ({} -> {}) at every N")
def _chain_rmse() -> tuple[str, ...]:
    return _pair_rmse("tau_chain_relative_rmse")


@claim("4.2", "exp07 timelike_pair_reconstruction_summary.csv",
       "the volume statistic uses {} pairs per N whereas the chain statistic "
       "uses {}")
def _pair_counts() -> tuple[str, ...]:
    rows = _csv("timelike_pair_reconstruction_summary.csv")
    pair = {float(r["pair_count"]) for r in rows}
    chain = {float(r["chain_pair_count"]) for r in rows}
    assert len(pair) == len(chain) == 1, (pair, chain)
    return _num(pair.pop()), _num(chain.pop())


@claim("4.2", "exp03 timelike_reconstruction_summary.csv",
       "the interval-cardinality formula returns tau = {} exactly")
def _identity_tau() -> tuple[str, ...]:
    values = {float(r["interval_tau_estimate"])
              for r in _csv("timelike_reconstruction_summary.csv")}
    assert len(values) == 1, values
    return (f"{values.pop():.1f}",)


@claim("4.2", "exp03 timelike_reconstruction_summary.csv",
       "({} at N = 200 -> {} at N = 2000)")
def _chain_convergence() -> tuple[str, ...]:
    return tuple(
        "{:.{p}f}".format(
            float(_one("timelike_reconstruction_summary.csv",
                       n_events=n)["chain_tau_abs_error"]), p=places)
        for n, places in ((200.0, 3), (2000.0, 4)))


@claim("4.2", "exp08 probe_pair_statistical_calibration_summary.csv",
       "the ratio of reconstruction RMSE to the predicted Poisson standard "
       "deviation is {}-{} across N (binomial {}-{})")
def _noise_ratios() -> tuple[str, ...]:
    rows = _csv("probe_pair_statistical_calibration_summary.csv")
    out: list[str] = []
    for key in ("rmse_to_poisson_abs_std_ratio", "rmse_to_binomial_abs_std_ratio"):
        values = _col(rows, key)
        out.append(f"{min(values):.2f}")
        out.append(_range_top(max(values)))
    return tuple(out)


def _range_top(value: float, places: int = 2) -> str:
    """The printed top of a stated RANGE, where Section 4.2 is not
    self-consistent: it rounds one of its two range tops to nearest
    (1.14031 -> 1.14) and the other outward (1.06492 -> 1.07, the
    honest bound for a range). Rather than pick a rule and edit prose
    that is defensible either way, EXACTLY these two candidates are
    accepted and exactly one of them must be printed.

    Stated precisely, because this file exists to stop verification
    claims from outrunning what is enforced: this is the one figure
    per range whose last digit may sit either side, so 1.07 -> 1.06
    would pass here. Nothing else does -- 1.08 and 1.05 both fail, as
    does any change to a lower bound, which uses no such tolerance."""

    step = 10.0 ** -places
    nearest = f"{value:.{places}f}"
    outward = f"{math.ceil(value / step) * step:.{places}f}"
    section = _section("4.2")
    printed = [c for c in dict.fromkeys((nearest, outward)) if c in section]
    assert len(printed) == 1, (value, nearest, outward, printed)
    return printed[0]


# =================== 4.3 R2 -- order + observer: radar time / distance

def _radar(ticks: float, key: str, places: int) -> str:
    row = _one("discrete_radar_reconstruction_summary.csv",
               N=300.0, tick_count=ticks)
    return f"{float(row[key]):.{places}f}"


@claim("4.3", "exp11 discrete_radar_reconstruction_summary.csv",
       "radar-time RMSE from {} at 16 ticks to {} at 128 ticks")
def _radar_time() -> tuple[str, ...]:
    return (_radar(16.0, "radar_time_rmse", 3),
            _radar(128.0, "radar_time_rmse", 4))


@claim("4.3", "exp11 discrete_radar_reconstruction_summary.csv",
       "radar-distance RMSE from {} to {}")
def _radar_distance() -> tuple[str, ...]:
    return (_radar(16.0, "radar_distance_rmse", 3),
            _radar(128.0, "radar_distance_rmse", 4))


@claim("4.3", "exp11 discrete_radar_reconstruction_summary.csv",
       "with the accessible fraction {} throughout")
def _radar_accessible() -> tuple[str, ...]:
    values = {float(r["accessible_fraction"])
              for r in _rows("discrete_radar_reconstruction_summary.csv", N=300.0)}
    assert len(values) == 1, values
    return (f"{values.pop():.1f}",)


# ================== 4.4 R3 -- + orientation: signed coords, Lorentz map

def _beta_rmse(beta: str, ticks: str, places: int) -> str:
    """Figure 2 panel (d) averages over N; the prose quotes that mean."""

    rows = [r for r in _csv("oriented_radar_lorentz_summary.csv")
            if r["beta"] == beta and r["tick_count"] == ticks]
    assert len(rows) == 3, (beta, ticks, len(rows))
    return f'{statistics.fmean(_col(rows, "fitted_beta_rmse")):.{places}f}'


@claim("4.4", "exp13 oriented_radar_lorentz_summary.csv",
       "it falls from {} at 32 ticks to {} at 128 ticks for beta = 0.3")
def _beta03() -> tuple[str, ...]:
    return (_beta_rmse("0.3", "32.0", 5), _beta_rmse("0.3", "128.0", 6))


@claim("4.4", "exp13 oriented_radar_lorentz_summary.csv",
       "and from {} to {} for beta = 0.6")
def _beta06() -> tuple[str, ...]:
    return (_beta_rmse("0.6", "32.0", 5), _beta_rmse("0.6", "128.0", 5))


# ================= 4.5 R4 -- + oriented atlas: transition-map consistency

def _atlas_mean(name: str, ticks: float, key: str, places: int,
                absolute: bool = False) -> str:
    values = _col(_rows(name, tick_count=ticks), key)
    if absolute:
        values = [abs(v) for v in values]
    return f"{statistics.fmean(values):.{places}f}"


@claim("4.5", "exp14 observer_atlas_transition_summary.csv",
       "the mean transition-map beta error falls {} -> {}")
def _atlas_beta() -> tuple[str, ...]:
    return tuple(
        _atlas_mean("observer_atlas_transition_summary.csv", t,
                    "fitted_beta_error", 4, absolute=True)
        for t in (32.0, 128.0))


@claim("4.5", "exp14 observer_atlas_transition_summary.csv",
       "and the invariant-interval RMSE {} -> {} (ticks 32 -> 128)")
def _atlas_invariant() -> tuple[str, ...]:
    return tuple(
        _atlas_mean("observer_atlas_transition_summary.csv", t,
                    "invariant_interval_rmse", 3)
        for t in (32.0, 128.0))


@claim("4.5", "exp14 observer_atlas_loop_summary.csv",
       "has a beta-composition error of {} -> {}")
def _atlas_loop() -> tuple[str, ...]:
    return tuple(
        _atlas_mean("observer_atlas_loop_summary.csv", t,
                    "beta_composition_error", 4, absolute=True)
        for t in (32.0, 128.0))


@claim("4.5", "exp15 exact_poincare_map_sanity.csv",
       "recovers the maps to machine precision (beta error ~{}, RMSE ~{})")
def _poincare_exact() -> tuple[str, ...]:
    """The quoted pair is the lab -> moving_pos transition, the check's
    first row. Its neighbour B -> C is NOT at machine precision -- that
    protocol's exact beta falls off the fit grid -- so the row is named
    rather than aggregated over."""

    row = next(r for r in _csv("exact_poincare_map_sanity.csv")
               if r["kind"] == "transition"
               and r["source_protocol"] == "A_lab"
               and r["target_protocol"] == "B_moving_pos")
    return f'{float(row["fitted_beta_error"]):.1e}', f'{float(row["rmse"]):.0e}'


# ================ 4.6 R5 -- + conformal measure: volume, coarse-graining

def _conformal_constant() -> list[dict[str, str]]:
    rows = [r for r in _csv("conformal_order_ambiguity_summary.csv")
            if r["profile"].startswith("constant_")]
    assert len(rows) == 3, len(rows)
    return sorted(rows, key=lambda r: float(r["scale"]))


@claim("4.6", "exp18 conformal_order_ambiguity_summary.csv",
       "under constant ({}) and sinusoidal rescalings the causal matrix is "
       "unchanged and the reconstructed dimension is identical ({}), while "
       "the proper-time ratio tracks {} and the volume ratio {}")
def _conformal() -> tuple[str, ...]:
    rows = _conformal_constant()
    dims = {float(r["estimated_dimension"])
            for r in _csv("conformal_order_ambiguity_summary.csv")}
    assert len(dims) == 1, dims

    def joined(key: str) -> str:
        return "/".join(f"{float(r[key]):g}" if float(r[key]) % 1 else
                        f"{float(r[key]):.1f}" for r in rows)

    return (joined("scale"), f"{dims.pop():.3f}",
            joined("proper_time_ratio"), joined("volume_ratio"))


def _profile_rows(profile: str) -> list[dict[str, str]]:
    """exp19's rows for one profile, N-sorted.

    The R5 evidence was re-pointed from `constant_1.5` to
    `sinusoidal_0.3` after a referee showed the constant profile's
    weighted relative RMSE is the flat profile's identically (a constant
    weight rescales estimate and truth alike), so that run demonstrates
    the global density -- ingredient M -- and nothing local. The
    sinusoidal profile is the package's only genuinely
    position-dependent weight, and its story is scaling, not bias: the
    profile is odd over the t-symmetric diamond, so its bias is forced
    toward zero by symmetry while its per-pair error is not."""

    rows = sorted((r for r in _csv("weighted_conformal_volume_summary.csv")
                   if r["profile"] == profile),
                  key=lambda r: float(r["N"]))
    assert len(rows) == 3, (profile, len(rows))
    return rows


def _track(profile: str, key: str, places: int = 3) -> str:
    return " -> ".join(f"{float(r[key]):.{places}f}"
                       for r in _profile_rows(profile))


@claim("4.6", "exp18 conformal_order_ambiguity_summary.csv",
       "the diamond volume moves by {}% (exp18's sinusoidal row)")
def _sin_volume_effect() -> tuple[str, ...]:
    row = [r for r in _csv("conformal_order_ambiguity_summary.csv")
           if r["profile"] == "sinusoidal_0.3"]
    assert len(row) == 1, len(row)
    return (f'{(float(row[0]["volume_ratio"]) - 1) * 100:.1f}',)


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "(unweighted volume bias {} to {})")
def _sin_unweighted_bias() -> tuple[str, ...]:
    rows = _profile_rows("sinusoidal_0.3")
    return (f'{float(rows[0]["unweighted_volume_bias"]):.4f}',
            f'{float(rows[-1]["unweighted_volume_bias"]):.4f}')


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "its relative RMSE stalls at {} across the tested N")
def _sin_unweighted_track() -> tuple[str, ...]:
    return (_track("sinusoidal_0.3", "unweighted_relative_rmse"),)


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "the weighted relative RMSE falls {}, tracking the flat-profile "
       "sampling floor ({}) to within {}%")
def _sin_weighted_track() -> tuple[str, ...]:
    weighted = [float(r["weighted_relative_rmse"])
                for r in _profile_rows("sinusoidal_0.3")]
    flat = [float(r["weighted_relative_rmse"]) for r in _profile_rows("flat")]
    worst = max(abs(a / b - 1) for a, b in zip(weighted, flat, strict=True))
    return (_track("sinusoidal_0.3", "weighted_relative_rmse"),
            _track("flat", "weighted_relative_rmse"),
            f"{worst * 100:.1f}")


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "volume-RMSE ratio grows {} with N")
def _sin_ratio_track() -> tuple[str, ...]:
    rows = _profile_rows("sinusoidal_0.3")
    ratios = [float(r["unweighted_volume_rmse"])
              / float(r["weighted_volume_rmse"]) for r in rows]
    assert ratios == sorted(ratios), ratios          # "grows" is asserted
    return (" -> ".join(f"{v:.1f}" for v in ratios),)


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "its unweighted bias ({}) is exactly the missing global factor")
def _const_identity() -> tuple[str, ...]:
    """The constant profile IS the flat computation: pinned as an
    identity, so the prose sentence stating it cannot outlive the data."""

    const = _profile_rows("constant_1.5")
    flat = _profile_rows("flat")
    for c_row, f_row in zip(const, flat, strict=True):
        assert abs(float(c_row["weighted_relative_rmse"])
                   - float(f_row["weighted_relative_rmse"])) < 1e-15
    bias = statistics.fmean(_col(const, "unweighted_volume_bias"))
    return (f"~ {bias:.2f}",)


@claim("4.6", "exp20 conformal_volume_exact_sanity.csv",
       "The analytic volume/proper-time formulas are verified to ~{}")
def _exact_sanity() -> tuple[str, ...]:
    worst = max(abs(v) for v in
                _col(_csv("conformal_volume_exact_sanity.csv"), "absolute_error"))
    return (f"1e-{-math.floor(math.log10(worst))}",)


def _thinning(keep: float, key: str, places: int) -> str:
    row = _one("thinning_coarse_graining_summary.csv", keep_probability=keep)
    return f"{float(row[key]):.{places}f}"


@claim("4.6", "exp23 thinning_coarse_graining_summary.csv",
       "reconstruction is stable (volume RMSE {} -> {}, dimension steady ~{})")
def _thinning_corrected() -> tuple[str, ...]:
    dims = _col(_csv("thinning_coarse_graining_summary.csv"), "dimension_mean")
    return (_thinning(0.25, "corrected_volume_rmse", 3),
            _thinning(1.0, "corrected_volume_rmse", 3),
            f"{statistics.fmean(dims):.1f}")


@claim("4.6", "exp23 thinning_coarse_graining_summary.csv",
       "blows up to RMSE {} (bias {}) at 25% retention")
def _thinning_uncorrected() -> tuple[str, ...]:
    return (_thinning(0.25, "uncorrected_volume_rmse", 3),
            _thinning(0.25, "uncorrected_volume_bias", 3))


# ============ 6.9 the certified oracle and the auxiliary audit

@claim("6.9", "p14_o3_volume.json",
       "`V \u2208 [{}, {}]`, relative half-width {}")
def _o3_certified() -> tuple[str, ...]:
    """The O3 certification the O4b audit consumes; Section 6.9 (formerly
    Section 9) prints its enclosure and half-width from the artifact."""

    r = _json("p14_o3_volume.json")["result"]
    return (f'{r["v_lo"]:.6f}', f'{r["v_hi"]:.6f}', f'{r["ratio"]:.6f}')


@claim("6.9", "p14_o3p_volume.json",
       "giving `V \u2208 [{}, {}]` — strictly inside O3")
def _o3p_recertified() -> tuple[str, ...]:
    """The half-width re-certification (O3') the mu = 0.1333 rung is
    gated against. "Strictly inside" is asserted, not assumed."""

    o3 = _json("p14_o3_volume.json")["result"]
    o3p = _json("p14_o3p_volume.json")["result"]
    assert o3["v_lo"] < o3p["v_lo"] and o3p["v_hi"] < o3["v_hi"]
    return (f'{o3p["v_lo"]:.6f}', f'{o3p["v_hi"]:.6f}')


# ================== 4.7 Horizon analogue -- Rindler inaccessibility

def _variance_structure() -> dict[str, float]:
    """S4's paired arms and S5's unpaired arms, from the stored
    per-reading arrays.

    Section 6.7 attributed its C2 grade to type D versus type N until a
    referee showed the difference is variance structure: C1 pairs on one
    point set and sees only the relation change, C2 does not and so also
    carries the between-reading event-count spread. Derived here rather
    than stored, so the identity is pinned and not the digits."""

    s4 = _json("p14_s4_results.json")
    s5 = _json("p14_s5_results.json")
    curved = s4["f_schwarzschild_lower"]["per_reading"]
    flat = s4["f_flat"]["per_reading"]
    paired = [c - f for c, f in zip(curved, flat, strict=True)]
    paired_sd = statistics.pstdev(paired)
    arm_c = s5["curved_f_lower"]["per_reading"]
    arm_f = s5["flat_f"]["per_reading"]
    pooled_sd = math.sqrt(
        (statistics.pstdev(arm_c) ** 2 + statistics.pstdev(arm_f) ** 2) / 2)
    unpaired_effect = abs(statistics.fmean(arm_c) - statistics.fmean(arm_f))
    return {
        "corr": statistics.correlation(curved, flat),
        "paired_sd": paired_sd,
        "pooled_sd": pooled_sd,
        "variance_removed": 1.0 - (paired_sd / pooled_sd) ** 2,
        "paired_sigma": abs(statistics.fmean(paired)) / paired_sd,
        "unpaired_sigma": unpaired_effect / pooled_sd,
    }


@claim("6.4", "p14_prereg_results.json c1.raw",
       "the relation fraction runs {} flat against {} curved, a factor of {}")
def _capstone_fractions() -> tuple[str, ...]:
    """The paper reported only the paired difference and its ratio to an
    operational margin, which hides the effect size; a referee asked for
    the two fractions themselves. Derived from the paired arrays."""

    raw = _json("p14_prereg_results.json")["c1"]["raw"]
    flat = statistics.fmean(raw["f_flat"])
    curved = statistics.fmean(raw["f_curved"])
    return (f"{flat:.6f}", f"{curved:.6f}", f"{curved / flat:.2f}")


@claim("6.4", "p14_prereg_results.json c2.metrics",
       "so at BA = {} it reaches {}, outside the range the quantity can take")
def _ba_out_of_range() -> tuple[str, ...]:
    """The frozen BA interval's upper end lies above 1. Derived, not
    typed, so the sentence cannot outlive the artifact it describes."""

    metrics = _json("p14_prereg_results.json")["c2"]["metrics"]
    assert metrics["ci_ba"][1] > 1.0, metrics["ci_ba"]
    return (f'{metrics["ba"]:.1f}', f'{metrics["ci_ba"][1]:.3f}')


@claim("6.4", "p14_prereg_results.json c2.raw",
       "Bonferroni-combined — gives [{}, 1.0] on the held-out halves")
def _ba_reference_interval() -> tuple[str, ...]:
    """A NON-frozen companion to the frozen BA interval, which is an
    unclipped Wald construction and reaches 1.014.

    Referee C established the boundary this sits inside: a different
    INTERVAL CONSTRUCTION on already-frozen data is descriptive and
    carries no verdict, whereas a new STATISTIC chosen after seeing the
    result would be a retrofitted gate. So this is derived here, from
    the stored arrays, by the same split the frozen code uses -- train
    on the first half of each arm, threshold at the midpoint of the
    train means -- and the interval is the one S5 later uses for this
    same quantity: exact Clopper-Pearson per arm at one-sided
    alpha/2 = 0.0125, Bonferroni-combined."""

    raw = _json("p14_prereg_results.json")["c2"]["raw"]
    curved, flat = raw["f_curved"], raw["f_flat"]
    half = len(curved) // 2
    threshold = (statistics.fmean(curved[:half])
                 + statistics.fmean(flat[:half])) / 2
    hits = (sum(1 for x in curved[half:] if x > threshold)
            + sum(1 for x in flat[half:] if x <= threshold))
    assert hits == 2 * half, hits          # complete separation on the test halves
    return (f"{0.0125 ** (1 / half):.6f}",)


@claim("6.7", "p14_s4_results.json + p14_s5_results.json",
       "the paired readings correlate at {} and pairing removes {}% of the "
       "spread, leaving a paired SD of {}")
def _s4_pairing() -> tuple[str, ...]:
    v = _variance_structure()
    return (f'{v["corr"]:.4f}', f'{100 * v["variance_removed"]:.1f}',
            f'{v["paired_sd"]:.6f}')


@claim("6.7", "p14_s4_results.json event_counts + per-reading arrays",
       "that count correlates with the paired difference at {}")
def _s4_count_is_not_the_nuisance() -> tuple[str, ...]:
    """The first draft of this passage blamed the event count, following
    one referee; a second referee recomputed and refuted it. The count
    explains ~1% of the variance in f and essentially none of the
    residual in delta, so the sentence had to name the shared point
    configuration instead. Pinned here so the refuted mechanism cannot
    creep back."""

    s4 = _json("p14_s4_results.json")
    delta = [c - f for c, f in zip(s4["f_schwarzschild_lower"]["per_reading"],
                                   s4["f_flat"]["per_reading"], strict=True)]
    r = statistics.correlation(delta, s4["event_counts"]["per_reading"])
    assert abs(r) < 0.10, r
    return (f"{r:.3f}",)


@claim("6.7", "p14_prereg_results.json + p14_s4_results.json",
       "the flat-arm relation fraction averages {} in the plane wave against "
       "{} here")
def _base_rates() -> tuple[str, ...]:
    """Why the unpaired SDs differ by an order of magnitude: a relation
    fraction near one half is near its most variable, and these two
    constructions sit an order of magnitude apart in it."""

    plane = _json("p14_prereg_results.json")["c1"]["raw"]["f_flat"]
    schwarzschild = _json("p14_s4_results.json")["f_flat"]["per_reading"]
    return (f"{statistics.fmean(plane):.4f}",
            f"{statistics.fmean(schwarzschild):.4f}")


@claim("6.7", "p14_s5_results.json + p14_s4_results.json",
       "its pooled SD of {} is an order of magnitude larger, and against "
       "that the same geometric effect is worth {} standard deviations "
       "rather than {}")
def _s5_pooled() -> tuple[str, ...]:
    v = _variance_structure()
    return (f'{v["pooled_sd"]:.6f}', f'{v["unpaired_sigma"]:.2f}',
            f'{v["paired_sigma"]:.1f}')


@claim("4.1", "exp09 longest_chain_calibration_summary.csv",
       "the same chains scored under the inclusive convention give {}, the "
       "two separated by the exact offset 2/sqrt(2 rho) = {}")
def _chain_convention_offset() -> tuple[str, ...]:
    """Section 3 declared the INCLUSIVE convention while the estimator
    reads the exclusive length; a referee caught the mismatch. The two
    residuals differ by a deterministic offset, so the inclusive figure
    is derived from the committed exclusive one rather than stored --
    pin the identity, not the digits."""

    rows = _chain_rows()
    exclusive = statistics.fmean(_col(rows, "tau_chain_error"))
    rho = {r["rho"] for r in rows}
    assert rho == {"300.0"}, rho
    offset = 2.0 / math.sqrt(2 * 300.0)
    return (f"{exclusive + offset:.3f}", f"{offset:.3f}")


@claim("4.1", "exp09 longest_chain_calibration_summary.csv",
       "the median chain carries {} elements including its endpoints")
def _chain_median_length() -> tuple[str, ...]:
    lengths = _col(_chain_rows(), "chain_length_including_endpoints")
    return (f"{statistics.median(lengths):.1f}",)


@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "zero false positives and zero false negatives, precision = recall = {}")
def _rindler_exact() -> tuple[str, ...]:
    rows = _csv("rindler_horizon_reconstruction_summary.csv")
    values = set(_col(rows, "precision")) | set(_col(rows, "recall"))
    assert len(values) == 1, values
    return (f"{values.pop():.1f}",)


@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "recovering a fraction {} to {} of it")
def _rindler_wedge_recall() -> tuple[str, ...]:
    """Recall against the IDEAL wedge, which is NOT what the table's
    `recall` column scores. That column is scored against the events
    finite clock coverage reaches, so it reads 1.0 while part of the
    wedge is never reconstructed; the honest figure is the ratio of the
    two accessible fractions. Section 4.7 said "the ideal-wedge
    classification is exact" until a referee caught the mismatch."""

    rows = _csv("rindler_horizon_reconstruction_summary.csv")
    ratios = [c / w for c, w in
              zip(_col(rows, "finite_coverage_accessible_fraction"),
                  _col(rows, "wedge_accessible_fraction"), strict=True)]
    return (f"{min(ratios):.2f}", f"{max(ratios):.2f}")


@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "(accessible fraction ~{})")
def _rindler_fraction() -> tuple[str, ...]:
    rows = _csv("rindler_horizon_reconstruction_summary.csv")
    return (f'{statistics.fmean(_col(rows, "wedge_accessible_fraction")):.2f}',)


@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "radar-time RMSE falls with tick resolution (e.g. {} -> {})")
def _rindler_rmse() -> tuple[str, ...]:
    """The illustration is the table's first configuration: the lowest
    acceleration at the smallest N."""

    return tuple(
        "{:.4f}".format(float(
            _one("rindler_horizon_reconstruction_summary.csv",
                 acceleration=1.5, N=600.0,
                 tick_count=t)["radar_time_rmse_finite_coverage"]))
        for t in (32.0, 128.0))


@claim("4.7", "exp17 inertial_vs_rindler_accessibility.csv",
       "only ~{} (ideal wedge) / ~{} (finite clock coverage) are accessible")
def _accessibility() -> tuple[str, ...]:
    rows = _csv("inertial_vs_rindler_accessibility.csv")
    return tuple(f"{statistics.fmean(_col(rows, key)):.2f}" for key in
                 ("rindler_wedge_accessible",
                  "rindler_finite_coverage_accessible"))


# ===================== 5. Negative results that bound the ladder

@claim("5", "exp12 single_observer_reflection_degeneracy.csv",
       "two targets at x = {} and x = {} return the identical single-observer "
       "distance {}, while the two-chain oriented protocol recovers the "
       "signed positions {} and {}")
def _reflection() -> tuple[str, ...]:
    rows = [r for r in _csv("single_observer_reflection_degeneracy.csv")
            if abs(float(r["true_x"])) == 0.1]
    assert len(rows) == 2, len(rows)
    pair = sorted(rows, key=lambda r: -float(r["true_x"]))
    distances = {f'{float(r["single_observer_radar_distance"]):.1f}' for r in pair}
    assert len(distances) == 1, distances
    return (*(f'{float(r["true_x"]):+.1f}' for r in pair),
            distances.pop(),
            *(f'{float(r["two_chain_signed_position"]):+.1f}' for r in pair))


@claim("5", "exp05 finite_speed_lattice_growth.csv",
       "finite-t counts differ before the calibration time (t = 5: {} vs {}) "
       "and agree at it by construction (t = 30: {} vs {})")
def _lattice_counts() -> tuple[str, ...]:
    out: list[str] = []
    for time in (5.0, 30.0):
        row = _one("finite_speed_lattice_growth.csv", time=time)
        out.append(_num(float(row["lattice_cumulative_count"])))
        out.append(_num(float(row["continuum_expected_count"])))
    return tuple(out)


@claim("5", "exp05 producer (lattice edge-direction census)",
       "its edges lie only along the two lightcone diagonals ({} each)")
def _lattice_edges() -> tuple[str, ...]:
    """The one legacy figure whose committed table cannot carry it:
    exp05 prints its edge-direction census to stdout rather than into
    `finite_speed_lattice_growth.csv`. It is re-derived from the same
    RNG-free producer call the script makes, at the script's own
    `T_STEPS`, so the number stays tied to the producer rather than
    being excused as unverifiable."""

    import numpy as np

    from causal_spacetime_lab.lattice import (
        edge_displacements,
        regular_lattice_causal_graph_1p1,
    )

    steps = int(re.search(
        r"^T_STEPS = (\d+)$",
        _text("experiments/exp05_finite_speed_lattice_counterexample.py"),
        re.M).group(1))
    graph = regular_lattice_causal_graph_1p1(steps)
    _, counts = np.unique(edge_displacements(graph), axis=0, return_counts=True)
    assert len(counts) == 2 and len(set(counts.tolist())) == 1, counts
    return (str(int(counts[0])),)


# =============================== 6.2 Pure-Weyl control construction

@claim("6.2", "p14_checks/p14_interval_volume_constant_a.py",
       "diamond volume sits {}% (`wT = 1`) to {}% (`wT = 2`) above flat")
def _diamond_excess() -> tuple[str, ...]:
    return tuple(f"{(_volume_ratio(wT) - 1.0) * 100:.1f}" for wT in (1.0, 2.0))


# ============================== 6.3 Finite-density detection design

@claim("6.3", "p14_prereg.md; p14_prereg_results.json",
       "slab `({})`, `w = {}`, expected {} events per sprinkling")
def _operating_point() -> tuple[str, ...]:
    rule = _text("docs/prereg/p14_prereg.md")
    slab = re.search(r"slab \(Δu,Δv,Δx,Δy\)=\(([\d., ]+)\)", rule).group(1)
    w = re.search(r"aniso-a1\.0 — w=([\d.]+),", rule).group(1)
    return (", ".join(part.strip() for part in slab.split(",")), w,
            _num(float(_json("p14_prereg_results.json")["e_n"])))


def _eps_delta_string() -> str:
    return f'{_json("p14_prereg_results.json")["eps_delta"] * 1e4:.3f}e-4'


@claim("6.3", "p14_prereg_results.json",
       "a frozen margin `epsilon_Delta = {}`")
def _eps_delta() -> tuple[str, ...]:
    return (_eps_delta_string(),)


@claim("6.3", "p14_prereg_results.json",
       "mean paired difference (n = {} sprinklings)")
def _n_c1() -> tuple[str, ...]:
    return (str(_json("p14_prereg_results.json")["n_c1"]),)


@claim("6.3", "p14_prereg_results.json",
       "flat ensembles (n = {} per unpaired arm)")
def _n_c2() -> tuple[str, ...]:
    return (str(_json("p14_prereg_results.json")["n_c2"]),)


@claim("6.3", "p14_prereg_preflight.json",
       "({}/{} joint-effect replicates;")
def _joint_certification() -> tuple[str, ...]:
    branch = _json("p14_prereg_preflight.json")["branches"]["joint_effect"]
    return (str(branch["counts"]["joint_confirmed"]), str(branch["reps"]))


@claim("6.3", "p14_prereg_preflight.json",
       "an exact Clopper-Pearson 95% lower bound of at least {};")
def _null_certification() -> tuple[str, ...]:
    """The stated floor is the weaker of the two null branches, rounded
    DOWN -- a floor that rounded up would not be a floor."""

    branches = _json("p14_prereg_preflight.json")["branches"]
    floors = [branches[b]["pass_ci95_exact"][0] for b in ("c1_null", "c2_null")]
    return (f"{math.floor(min(floors) * 100) / 100:.2f}",)


# ==================================== 6.4 Preregistered result

@claim("6.4", "p14_prereg_results.json",
       "| C1 paired ensemble mean, n = {} | confirmed | mean {} [{}, {}]; "
       "lower end {}x the margin {} |")
def _c1_row() -> tuple[str, ...]:
    art = _json("p14_prereg_results.json")
    metrics = art["c1"]["metrics"]
    return (str(art["n_c1"]), f'{metrics["mean"]:.7f}',
            *(f"{v:.7f}" for v in metrics["ci"]),
            f'{metrics["ci"][0] / art["eps_delta"]:.0f}', _eps_delta_string())


@claim("6.4", "p14_prereg_results.json",
       "| C2 classifier replication, n = {}/arm | confirmed | separation "
       "s = {} [{}, {}]; AUC = {} [{}, {}]; balanced accuracy = {} [{}, {}] |")
def _c2_row() -> tuple[str, ...]:
    art = _json("p14_prereg_results.json")
    metrics = art["c2"]["metrics"]

    def bound(value: float) -> str:
        """A bound pinned at the exact boundary prints as the boundary;
        the interior bound carries the frozen six decimals."""

        return f"{value:.1f}" if value in (0.0, 1.0) else f"{value:.6f}"

    return (str(art["n_c2"]),
            f'{metrics["s"]:.3f}', *(f"{v:.3f}" for v in metrics["ci_s"]),
            f'{metrics["auc"]:.1f}', *(bound(v) for v in metrics["ci_auc"]),
            f'{metrics["ba"]:.1f}', *(f"{v:.3f}" for v in metrics["ci_ba"]))


@claim("6.4", "p14_prereg_results.json",
       "(the minimum curved-arm value exceeds the maximum flat-arm value in "
       "{} draws per arm)")
def _complete_separation() -> tuple[str, ...]:
    art = _json("p14_prereg_results.json")
    raw = art["c2"]["raw"]
    assert min(raw["f_curved"]) > max(raw["f_flat"])
    assert len(raw["f_curved"]) == len(raw["f_flat"]) == art["n_c2"]
    return (str(art["n_c2"]),)


# ========== 6.6 What the plane-wave result alone does not establish

@claim("6.7", "p14_s1_cost.json",
       "about {} ms per pair on the tested solver, patch, and tolerance, "
       "roughly {}x the plane-wave predicate")
def _s1_price() -> tuple[str, ...]:
    art = _json("p14_s1_cost.json")
    rung = next(e for e in art["ladder"] if e["tol"] == art["default_tol"])
    return (f'{rung["us_per_pair"] / 1000:.2f}',
            f'{round(art["price_ratio_at_default_tol"], -1):.0f}')


# ================= 6.7 Type-D extension: S4 confirmation, S5 detection

@claim("6.7", "p14_s1_cost.json; p14_s4_results.json",
       "(`M = {}`, exterior shell `r` in `[{}, {}]`, polar cap, "
       "coordinate-time extent {}), with `N ~ Poisson({})` events per reading")
def _s4_domain() -> tuple[str, ...]:
    domain = _json("p14_s1_cost.json")["domain"]
    return (_num(domain["m"]), *(_num(v) for v in domain["r_shell"]),
            _num(domain["t_extent"]),
            _num(float(_json("p14_s4_results.json")["params"]["e_n"])))


@claim("6.7", "p14_s1_cost.json",
       "the S1 predicate at tolerance {} with escalation to {}")
def _s4_tolerance() -> tuple[str, ...]:
    art = _json("p14_s1_cost.json")
    finest = min(e["tol"] for e in art["ladder"])
    return tuple(f"1e-{-round(math.log10(t))}"
                 for t in (art["default_tol"], finest))


@claim("6.7", "p14_s4_results.json; p14_s3_probe_results.json",
       "the frozen threshold `eps_det = {}`, about {}% of the exploration "
       "anchor")
def _eps_det() -> tuple[str, ...]:
    eps = _json("p14_s4_results.json")["margins"]["eps_det"]
    anchor = abs(_json("p14_s3_probe_results.json")["delta_lower"]["mean"])
    return (_num(eps), f"{eps / anchor * 100:.0f}")


@claim("6.7", "p14_s4_results.json",
       "inside `+-eps_rep = +-{}`, one exploration reading-SD")
def _eps_rep() -> tuple[str, ...]:
    return (_num(_json("p14_s4_results.json")["margins"]["eps_rep"]),)


@claim("6.7", "p14_s4_schwarzschild_c1.md",
       "({}/{} on every branch, exact Clopper-Pearson 95% lower bound {})")
def _s4_power() -> tuple[str, ...]:
    rule = _text("docs/prereg/p14_s4_schwarzschild_c1.md")
    rows = re.findall(r"\| (\d+)/(\d+) = 1\.00000 \| (0\.\d+) \|", rule)
    assert len(rows) == 5 and len(set(rows)) == 1, rows
    hit, total, bound = rows[0]
    return (hit, total, f"{float(bound):.4f}")


@claim("6.7", "p14_s4_freeze_manifest.json",
       "with an {}-file content-addressed manifest verified at entry and at "
       "exit")
def _s4_manifest() -> tuple[str, ...]:
    return (str(len(_json("p14_s4_freeze_manifest.json")["files"])),)


@claim("6.7", "p14_s4_results.json",
       "| A: C1 detection (primary) | identified CI95 top < {} | "
       "CI95 [{}, {}] | pass ({}x margin) |")
def _gate_a() -> tuple[str, ...]:
    art = _json("p14_s4_results.json")
    lo, hi = art["identified_ci95"]
    eps = art["margins"]["eps_det"]
    assert art["gate_a"] is True and hi < -eps
    return (_num(-eps), f"{lo:.6f}", f"{hi:.6f}", f"{abs(hi) / eps:.0f}")


@claim("6.7", "p14_s4_results.json",
       "| B: replication (secondary) | Welch CI95 of S4-S3 inside +-{} | "
       "[{}, {}] | REPLICATED |")
def _gate_b() -> tuple[str, ...]:
    art = _json("p14_s4_results.json")
    lo, hi = art["welch_identified_ci95"]
    assert art["gate_b"] == "REPLICATED"
    return (_num(art["margins"]["eps_rep"]), f"{lo:.6f}", f"{hi:+.6f}")


def _s5_floor() -> str:
    """The S5 rule states its floor as chance plus the frozen margin,
    the same number for the primary and the secondary gate."""

    margins = _json("p14_s5_results.json")["margins"]
    assert margins["eps_auc"] == margins["eps_ba"], margins
    return f'{0.5 + margins["eps_auc"]:.2f}'


@claim("6.7", "p14_s5_results.json (margin); p14_s5_schwarzschild_c2.md",
       "AUC {} as the minimum practically useful single-poset discrimination")
def _s5_declared_floor() -> tuple[str, ...]:
    floor = _s5_floor()
    assert f"AUC {floor}" in _text("docs/prereg/p14_s5_schwarzschild_c2.md")
    return (floor,)


@claim("6.7", "p14_s5_results.json",
       "; {} readings each by a pre-declared minimum-n rule")
def _s5_readings() -> tuple[str, ...]:
    return (_num(float(_json("p14_s5_results.json")["params"]["n_arm"])),)


@claim("6.7", "p14_s5_results.json",
       "| Primary: AUC (DeLong) | CI95 lower > {} | AUC {}, CI95 [{}, {}] | "
       "**DETECTED** |")
def _s5_primary() -> tuple[str, ...]:
    art = _json("p14_s5_results.json")
    auc, lo, hi = art["auc"]["auc_lower_series"]
    assert art["outcome"] == "DETECTED"
    return (_s5_floor(), f"{auc:.4f}", f"{lo:.4f}", f"{hi:.4f}")


@claim("6.7", "p14_s5_results.json",
       "| Secondary: out-of-sample BA | joint CI95 lower > {} | BA {}, "
       "CI95 [{}, {}] | pass |")
def _s5_secondary() -> tuple[str, ...]:
    ba = _json("p14_s5_results.json")["ba"]
    assert ba["pass"] is True
    return (_s5_floor(), f'{ba["ba"]:.4f}',
            *(f"{v:.4f}" for v in ba["ci95_cp_bonferroni"]))


@claim("6.7", "p14_s5_results.json",
       "does not (AUC ≈ {};")
def _s5_prose_auc() -> tuple[str, ...]:
    return (f'{_json("p14_s5_results.json")["auc"]["auc_lower_series"][0]:.3f}',)


# ---------------------------- the claim boundary quotes the same figures

#: The claim boundary restates the Sections 6.4/6.7 figures for the
#: capstone and the type-D extension. Nothing kept the two documents in
#: step, which is how one drifted point estimate came to be printed in
#: both; these hold them to the SAME artifacts as the manuscript.
BOUNDARY_CLAIMS: list[Claim] = []


def boundary(source: str, context: str) -> Callable:
    def wrap(fn: Callable[[], tuple[str, ...]]) -> Claim:
        made = Claim("claim_boundary.md", source, context, fn)
        BOUNDARY_CLAIMS.append(made)
        return made
    return wrap


@boundary("p14_prereg_results.json",
          "{} [{}, {}] vs epsilon_Delta = {}")
def _boundary_c1() -> tuple[str, ...]:
    metrics = _json("p14_prereg_results.json")["c1"]["metrics"]
    return (f'{metrics["mean"]:.7f}', *(f"{v:.7f}" for v in metrics["ci"]),
            _eps_delta_string())


@boundary("p14_prereg_results.json",
          "(s = {} [{}, {}]; AUC = {} [{}, {}]; BA = {} [{}, {}])")
def _boundary_c2() -> tuple[str, ...]:
    n, *rest = _c2_row.derive()
    assert n == str(_json("p14_prereg_results.json")["n_c2"])
    return tuple(rest)


@boundary("p14_s4_results.json",
          "(identified CI95 [{}, {}] vs threshold {})")
def _boundary_s4_gate_a() -> tuple[str, ...]:
    threshold, lo, hi, _margin = _gate_a.derive()
    return (lo, hi, threshold)


@boundary("p14_s4_results.json",
          "(Welch CI95 of the difference inside +-{})")
def _boundary_s4_gate_b() -> tuple[str, ...]:
    return (_gate_b.derive()[0],)


@boundary("p14_s5_results.json",
          "the independently declared {} threshold: AUC {}, DeLong CI95 "
          "[{}, {}]")
def _boundary_s5_auc() -> tuple[str, ...]:
    return _s5_primary.derive()


@boundary("p14_s5_results.json",
          "out-of-sample BA {}, CP-Bonferroni CI95 [{}, {}]")
def _boundary_s5_ba() -> tuple[str, ...]:
    _floor, *rest = _s5_secondary.derive()
    return tuple(rest)


# ------------------------------------------------------ the contracts

@pytest.mark.parametrize(
    "claim_", CLAIMS, ids=lambda c: f"{c.section}{c.derive.__name__}")
def test_the_manuscript_sentence_is_the_artifact_value(claim_: Claim):
    """Every quantitative claim in Sections 4-6.7: the value is read out
    of its artifact, formatted the way the manuscript formats it, and
    the resulting SENTENCE must occur in the section that carries it."""

    values = claim_.derive()
    assert claim_.context.count("{}") == len(values), (claim_.section, values)
    expected = claim_.context.format(*values)
    assert expected in _section(claim_.section), (
        claim_.section, claim_.source, expected)


def test_each_claim_is_pinned_to_a_single_section():
    """A sentence that also occurs elsewhere would let a value drift
    into the wrong section and still pass the check above. Every
    section is searched, including the unnumbered-subsection ones (4,
    5, 6): an earlier draft skipped those, which silently exempted
    Section 5's claims from this guard entirely."""

    for claim_ in CLAIMS:
        expected = claim_.context.format(*claim_.derive())
        hits = [s for s in SECTIONS if expected in SECTIONS[s]]
        assert hits == [claim_.section], (claim_.section, hits)


# ------------------------------------------- claims that are not numbers

def test_the_qualitative_claims_hold_in_the_artifacts():
    """Sections 4-5 also make statements ABOUT the tables that carry no
    printed figure. They are contracts too: each is the reason a number
    above is allowed to be read the way it is read."""

    dim = _csv("dimension_reconstruction_summary.csv")
    for spacetime_dim in (2.0, 3.0, 4.0):
        rmse = {float(r["N"]): float(r["rmse"]) for r in dim
                if float(r["spacetime_dim"]) == spacetime_dim}
        assert rmse[2400.0] < rmse[300.0], spacetime_dim
    # "non-monotonic finite-sample fluctuations for dimensions 3 and 4"
    for spacetime_dim, monotone in ((2.0, True), (3.0, False), (4.0, False)):
        series = [float(r["rmse"]) for r in dim
                  if float(r["spacetime_dim"]) == spacetime_dim]
        assert (series == sorted(series, reverse=True)) is monotone, spacetime_dim

    # 4.2: the fixed-interval sanity check is a normalization identity
    for row in _csv("timelike_reconstruction_summary.csv"):
        assert float(row["interval_count"]) == float(row["n_events"])
        assert float(row["rho"]) == (float(row["n_events"])
                                     / float(row["diamond_volume"]))

    # 4.2: the volume estimator sits below the chain estimator at EVERY N
    for row in _csv("timelike_pair_reconstruction_summary.csv"):
        assert (float(row["tau_volume_relative_rmse"])
                < float(row["tau_chain_relative_rmse"]))

    # 4.6: positive rescaling leaves the causal matrix alone
    for row in _csv("conformal_order_ambiguity_summary.csv"):
        assert float(row["causal_matrix_changed"]) == 0.0

    # 4.6: the R5 story is scaling, and the constant profile is M-only.
    # (i) the sinusoidal unweighted error is a floor: its last two values
    # agree within 1% and sit far above the weighted error; (ii) the
    # weighted error recovers the flat sampling floor within 2%; (iii)
    # the constant profile's weighted relative RMSE is the flat one
    # identically, and its unweighted bias is the missing global factor.
    def _p46(profile, key):
        return [float(r[key])
                for r in sorted((r for r in
                                 _csv("weighted_conformal_volume_summary.csv")
                                 if r["profile"] == profile),
                                key=lambda r: float(r["N"]))]
    u_sin = _p46("sinusoidal_0.3", "unweighted_relative_rmse")
    w_sin = _p46("sinusoidal_0.3", "weighted_relative_rmse")
    w_flat = _p46("flat", "weighted_relative_rmse")
    assert abs(u_sin[1] - u_sin[2]) / u_sin[2] < 0.01, u_sin
    assert u_sin[2] > 2 * w_sin[2], (u_sin[2], w_sin[2])
    assert all(abs(a / b - 1) < 0.02
               for a, b in zip(w_sin, w_flat, strict=True)), (w_sin, w_flat)
    w_const = _p46("constant_1.5", "weighted_relative_rmse")
    assert all(abs(a - b) < 1e-15
               for a, b in zip(w_const, w_flat, strict=True))
    assert all(b < -0.25
               for b in _p46("constant_1.5", "unweighted_volume_bias"))

    # 4.7: exact wedge classification, and Rindler access a strict subset
    for row in _csv("rindler_horizon_reconstruction_summary.csv"):
        assert float(row["false_positive"]) == float(row["false_negative"]) == 0.0
    access = _csv("inertial_vs_rindler_accessibility.csv")
    assert {float(r["inertial_accessible"]) for r in access} == {1.0}
    assert not [r for r in access if float(r["rindler_accessible"]) == 1.0
                and float(r["inertial_accessible"]) == 0.0]

    # 5: the lattice shares the continuum's leading quadratic growth
    last = _one("finite_speed_lattice_growth.csv", time=30.0)
    assert (float(last["lattice_cumulative_count"])
            == float(last["continuum_expected_count"]))


def test_the_preregistered_verdicts_are_the_artifacts_verdicts():
    """Sections 6.4 and 6.7 state verdicts, censuses and frozen
    sentences, not only figures."""

    plane = _json("p14_prereg_results.json")
    assert plane["stage_positive"] is True
    for key in ("c1", "c2"):
        assert plane[key]["verdict"] == "confirmed"
        for arm in plane[key]["ambiguity"].values():
            assert arm["ambiguous"] == arm["escalated"] == 0

    s4 = _json("p14_s4_results.json")
    assert s4["verdict"] == "CONFIRMED"
    assert s4["ambiguity"] == {"ambiguous": 0, "escalated": 1}

    s5 = _json("p14_s5_results.json")
    assert s5["outcome"] == "DETECTED"
    assert s5["ambiguity"] == {"ambiguous": 0, "escalated": 0}
    assert s5["auc"]["auc_lower_series"] == s5["auc"]["auc_upper_series"]

    section = _section("6.7")
    assert "Stage verdict: CONFIRMED." in section
    assert "zero ambiguous pairs and one escalated pair" in section
    assert "Zero ambiguous and zero escalated pairs" in section


def test_each_frozen_sentence_reaches_the_paper_as_its_rendering():
    """The guard that survived dropping the Korean quotes.

    Both halves matter and neither implies the other. The first pins
    every original against the artifact that holds it, so re-freezing a
    sentence cannot slip past while its translation goes stale. The
    second requires that translation to be printed, in the section that
    carries the verdict -- not merely somewhere in the paper, which
    would let a claim drift into a section whose evidence does not
    support it."""

    originals = _frozen_sentences()
    assert len(originals) == len(FROZEN_RENDERINGS), (
        len(originals), len(FROZEN_RENDERINGS))

    for original, rendering in zip(originals, FROZEN_RENDERINGS, strict=True):
        assert original == rendering.original, original[:40]
        assert rendering.english in _section(rendering.section), (
            rendering.section, rendering.english[:60])


def test_the_latex_rendition_prints_each_frozen_rendering():
    """FROZEN_RENDERINGS holds the manuscript to the artifacts; this
    holds the LaTeX rendition to the same six English renderings, so the
    two renderings cannot drift apart on the paper's load-bearing
    sentences. The LaTeX escapes nothing inside these sentences except
    dashes, so an em-dash-normalized, whitespace-flattened containment
    check is exact enough to bite."""

    sections_dir = _PAPER / "latex" / "sections"
    flat = " ".join(
        " ".join(path.read_text(encoding="utf-8").split())
        for path in sorted(sections_dir.glob("*.tex")))
    flat = flat.replace("---", " — ").replace("--", "–")
    # the two escapes the renderings carry in the LaTeX
    flat = flat.replace("$\\varepsilon_\\Delta$", "epsilon_Delta")
    flat = flat.replace("\\%", "%")
    flat = " ".join(flat.split())
    for rendering in FROZEN_RENDERINGS:
        needle = " ".join(rendering.english.split())
        assert needle in flat, (rendering.section, needle[:60])


def test_the_manuscript_prints_no_korean():
    """Appendix B says the paper prints none; this is what makes that
    true tomorrow. The frozen originals belong in their artifacts and in
    `FROZEN_RENDERINGS`, where they are checked -- not in a manuscript
    whose readers cannot check them, and whose target journal sets its
    body and references in roman characters."""

    hangul = re.findall(r"[가-힣]+",
                        MANUSCRIPT.read_text(encoding="utf-8"))
    assert not hangul, hangul[:5]


# ----------------------------------------------- the completeness guard

class Exclusion(NamedTuple):
    """One numeral occurrence no artifact produces as a result, bound
    to the WORDING that justifies it rather than to the bare literal.

    `phrase` is the text as printed, quoted from the section residue;
    `count` is how many times that wording occurs there. The reason has
    to hold for the phrase, so if the phrase goes the excuse goes.
    """

    section: str
    phrase: str
    count: int
    reason: str


#: The honest edge of this contract, arrived at over two review rounds.
#:
#: R1 killed the first form, a global `set[str]`: one justified `300`
#: licensed every other `300` in Sections 4-6.7. R2 killed the second,
#: `(section, literal) -> count`: the count survives a change of
#: MEANING inside one section, so 4.3's `at N = 300 ...` row selector
#: could become `the run took 300 ms` -- an unverified new result
#: wearing a permitted occurrence's clothes -- and still pass.
#:
#: Binding to the phrase closes both. An excuse now covers exactly the
#: sentence it was written about; rewrite that sentence and the excuse
#: stops applying, whatever the literal or the count.
ACCEPTED_EXCLUSIONS: tuple[Exclusion, ...] = (
    Exclusion("4.1", "Endpoint RMSE is lower at N = 2400 than at N = 300", 1,
              "N selectors of a comparison asserted as an inequality, not printed"),
    Exclusion("4.1", "fluctuations for dimensions 3 and 4", 1,
              "dimension labels of a claim asserted qualitatively, not printed"),
    Exclusion("4.1", "scaling L ~ sqrt(2 rho) tau", 1,
              "the Brightwell-Gregory scaling law, a formula constant"),
    Exclusion("4.2", "tau_est = sqrt(2K/rho)", 1,
              "the estimator's defining formula, a constant in it"),
    Exclusion("4.3", "at N = 300 the error falls roughly by half per doubling", 1,
              "N selector naming which exp11 rows the tick series is read from"),
    Exclusion("4.4", "dimension estimates lie near D = 2, 3, 4", 1,
              "Figure 2 caption: the panel's dimension labels"),
    Exclusion("4.4", "endpoint RMSE is lower at N = 2400 than N = 300", 1,
              "Figure 2 caption restating 4.1's endpoint-RMSE inequality"),
    Exclusion("4.6", "The constant-1.5 profile in the same table", 1,
              "the profile's label, a row selector"),
    Exclusion("4.6", "`Omega(t) = 1 + 0.3 sin(pi t/T)`", 1,
              "the tested profile's definition, a construction setting"),
    Exclusion("4.6", "by a factor of 3.45 across the diamond (0.49 to 1.69)", 1,
              "arithmetic of the declared amplitude: (1 +- 0.3)^2 = 0.49 "
              "and 1.69, ratio 3.45"),
    Exclusion("4.6", "(constant-1.5 conformal profile)", 0,
              "Figure 3 caption: the exp19 profile label"),
    Exclusion("6.2", "profile `A(u)(x^2 - y^2)`", 1,
              "the construction's profile, its exponents a setting not a result"),
    Exclusion("6.2", "`det g = -1`", 2,
              "an exact property of the Brinkmann construction"),
    Exclusion("6.2", "read with `A = 0`", 1,
              "the control reading of the paired design, a construction setting"),
    Exclusion("6.2", "`V_A/V_0 = 1 + (wT)^4/252 + O((wT)^8)`", 1,
              "the analytic volume response of the construction, an exact "
              "closed form and not a measurement"),
    Exclusion("6.2", "`1/pi = 0.318` of the way to the first conjugate point", 1,
              "where the frozen slab sits inside its own validity window, a "
              "construction setting"),
    Exclusion("6.3", "the 95% Student-t interval", 1,
              "the confidence level of every stated interval, a frozen convention"),
    Exclusion("6.4", "Result (95% CI)", 1,
              "the confidence level in the result-table header, same convention"),
    Exclusion("6.7", "`sqrt(-g) = r^2 sin(theta)`", 1,
              "the measure identity the extension rests on, an exact formula"),
    Exclusion("6.7", "exactly as `det g = -1` served the plane wave", 1,
              "6.2's construction identity, quoted back as the analogy"),
    Exclusion("6.7", "AUC CI95 lower bound above 0.60", 1,
              "English rendering of a frozen Korean sentence; the original "
              "is bound to this wording by FROZEN_RENDERINGS"),
    Exclusion("6.7", "the joint 95% lower bound of out-of-sample balanced "
                     "accuracy exceeds 0.60", 1,
              "English rendering of a frozen Korean sentence; the original "
              "is bound to this wording by FROZEN_RENDERINGS"),
)


_CROSS_REFERENCE = re.compile(
    r"Sections? \d+(?:\.\d+)?(?:[-–]\d+(?:\.\d+)?)?"  # Section 4.6, Sections 4-6.8
    r"|Figures? \d+|Appendix [A-Z]|Table \d+"              # figure/appendix pointers
    r"|\bP\d+(?:[-–]P\d+)?\b|\bexp\d+\b"              # programme/producer ids
    r"|\b\d\+\dD\b|\b\dD\b"                                # 1+1D, 4D
    r"|\([a-d]\)"                                          # figure panel letters
)

_NUMBER = re.compile(r"(?<![\w.])[-+]?\d(?:[\d,]*\d)?(?:\.\d+)?(?:[eE][-+]?\d+)?")

_SCANNED = ("4", "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "5",
            "6", "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7")


class Rendering(NamedTuple):
    """One frozen sentence, and what the manuscript prints for it."""

    section: str      # where the rendering must appear
    original: str     # the frozen sentence, as its artifact holds it
    english: str      # the rendering printed in its place


#: The manuscript used to quote each frozen sentence in Korean, byte for
#: byte, and let an English gloss follow it. It no longer does: the bytes
#: were unverifiable to the paper's readers, the gloss the reader
#: actually relied on was the half the paper marked NON-frozen, and IOP
#: house style is roman-only. Appendix B now records the originals'
#: language and location once, and this table carries what removing them
#: would otherwise have cost -- each original, pinned here against what
#: `_frozen_sentences()` reads from the artifact, bound to the rendering
#: the manuscript prints in its place. A re-freeze fails the first half;
#: a translation that softens, strengthens or outruns its original fails
#: the second. Order matches `_frozen_sentences()`.
FROZEN_RENDERINGS = (
    Rendering(
        "6.4", "paired ensemble 평균 이동이 ε_Δ를 넘는다",
        "the paired ensemble mean shift exceeds epsilon_Delta"),
    Rendering(
        "6.4", "P3-C 분리를 독립적으로 재현했다.",
        "the probe chain's separation is independently reproduced"),
    Rendering(
        "6.7",
        "동결된 Schwarzschild 좌표·도메인(M=1, r∈[10,20], 극관각 캡 1.0, T=40)의 "
        "공통 측도 위에서, paired 앙상블 평균 이동의 identified CI95가 동결 방향으로 "
        "검출문턱 ε_det = 0.0036을 초과했다 — 유한밀도 인과 census가 type D 진공 "
        "곡률의 빛원뿔 변형을 C1-급으로 검출했다(프로그램 내부 진술).",
        "on the frozen Schwarzschild coordinates and domain, the paired "
        "ensemble mean shift exceeds the frozen detection threshold in the "
        "frozen direction — a finite-density causal census detects the "
        "light-cone deformation of type-D vacuum curvature at C1 grade; "
        "a program-internal statement"),
    Rendering(
        "6.7",
        "S4 블록과 S3 탐색 블록의 독립 두-표본 차이의 identified Welch CI95가 "
        "±ε_rep = ±0.0012 안에 들어, 탐색 효과가 정량적으로 재현됐다.",
        "the independent two-sample difference between the S4 and S3 blocks "
        "lies inside the replication band: the exploration effect is "
        "quantitatively reproduced"),
    Rendering(
        "6.7",
        "동결된 Schwarzschild 도메인·밀도에서, 단일 causal set의 global relation "
        "fraction은 flat/Schwarzschild 앙상블을 우연 수준보다 판별하는 정보를 "
        "운반한다 (AUC CI95 하한 > 0.60, 프로그램 내부 진술).",
        "on the frozen Schwarzschild domain and density, the global relation "
        "fraction of a single causal set carries information that "
        "discriminates flat from Schwarzschild ensembles above chance — AUC "
        "CI95 lower bound above 0.60; a program-internal statement"),
    Rendering(
        "6.7",
        "out-of-sample balanced accuracy의 결합 95% 하한이 0.60을 넘어, "
        "학습-외 판별이 확인됐다 (secondary).",
        "the joint 95% lower bound of out-of-sample balanced accuracy "
        "exceeds 0.60: out-of-training discrimination holds; a secondary "
        "verdict"),
)


def _frozen_sentences() -> tuple[str, ...]:
    """The frozen verdict sentences as the artifacts hold them, in the
    order `FROZEN_RENDERINGS` pairs them. Read from the artifacts, never
    copied, so a re-frozen sentence cannot leave a stale entry there."""

    plane = _json("p14_prereg_results.json")
    return (*(plane[k]["sentence"] for k in ("c1", "c2")),
            *_json("p14_s4_results.json")["sentences"],
            *_json("p14_s5_results.json")["sentences"])


def _residue(number: str) -> str:
    """What a section still says once everything already accounted for
    is struck out: each claim's matched sentence (removed ONCE, so a
    duplicated sentence survives), and pointers like "Section 4.6" that
    reference rather than measure.

    The frozen sentences used to be struck here too, because the
    manuscript quoted them in Korean and their numerals are the RULE's,
    not the paper's. It no longer does (`FROZEN_RENDERINGS`), so the
    numerals now reaching this scan are the ones the ENGLISH renderings
    print -- and those are the paper's own responsibility, accounted for
    by `ACCEPTED_EXCLUSIONS` like any other figure."""

    body = SECTIONS.get(number, "")
    for claim_ in CLAIMS:
        if claim_.section != number:
            continue
        matched = claim_.context.format(*claim_.derive())
        assert matched in body, (number, matched)
        body = body.replace(matched, " ", 1)
    return _CROSS_REFERENCE.sub(" ", body)


def _unexplained(number: str) -> list[str]:
    """Numerals left in a section once its claims AND its exclusions'
    exact wordings are struck out. Anything here is a figure nothing
    in this file has accounted for."""

    body = _residue(number)
    for excluded in ACCEPTED_EXCLUSIONS:
        if excluded.section != number:
            continue
        found = body.count(excluded.phrase)
        assert found == excluded.count, (
            number, excluded.phrase, found, excluded.count)
        body = body.replace(excluded.phrase, " ")
    return _NUMBER.findall(body)


def test_every_number_is_derived_or_explicitly_excluded():
    """The guard the [P2] review asked for, in the form two rounds of
    review drove it to. Every numeral printed in Sections 4-6.7 sits
    either inside a sentence some claim re-derived from an artifact, or
    inside a specific quoted phrase carrying its reason. A figure added
    anywhere in the range lands here, not in a silent gap."""

    leftover = {number: found for number in _SCANNED
                if (found := _unexplained(number))}
    assert not leftover, leftover


def test_no_exclusion_is_stale_or_unreasoned():
    """`_unexplained` already fails when a phrase stops occurring the
    stated number of times, which is staleness. What is left to check
    is that an exclusion is the kind of thing it claims to be: it must
    excuse an actual numeral, and it must say why."""

    for excluded in ACCEPTED_EXCLUSIONS:
        assert _NUMBER.search(excluded.phrase), excluded.phrase
        assert excluded.reason.strip(), excluded.phrase
        assert excluded.section in _SCANNED, excluded.section
    for number in _SCANNED:
        _unexplained(number)


def test_the_guard_binds_a_number_to_its_wording():
    """Guard the guard, against the hole R2 named. 4.3's `300` is
    excused as a row selector. Repoint that same occurrence at a new
    unverified result -- same section, same literal, same count -- and
    the excuse must stop covering it."""

    selector = next(e for e in ACCEPTED_EXCLUSIONS if e.section == "4.3")
    body = _residue("4.3").replace(selector.phrase, "the run took 300 ms", 1)
    for excluded in ACCEPTED_EXCLUSIONS:
        if excluded.section == "4.3":
            body = body.replace(excluded.phrase, " ")
    assert _NUMBER.findall(body) == ["300"]


def test_the_claim_table_covers_every_section_and_legacy_table():
    """Guard the guard: a claim table that quietly lost a section, or a
    committed legacy table that no claim reads, would leave exactly the
    hole this file exists to close."""

    covered = {c.section for c in CLAIMS}
    assert covered == {"4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7",
                       "5", "6.2", "6.3", "6.4", "6.7", "6.9"}, sorted(covered)

    read = " ".join(c.source for c in CLAIMS)
    unread = [p.name for p in sorted(_FIG_DATA.glob("*.csv"))
              if p.name not in read]
    # exp06's proxy table is cited in Section 5 as explicitly NOT a
    # validated estimator, and the manuscript quotes no figure from it.
    assert unread == ["spacelike_distance_proxy_summary.csv"], unread


@pytest.mark.parametrize(
    "claim_", BOUNDARY_CLAIMS, ids=lambda c: c.derive.__name__)
def test_the_claim_boundary_restates_the_artifact_values(claim_: Claim):
    """Whatever the claim boundary repeats from Sections 6.4 and 6.7 is
    re-derived from the same artifact, so the manuscript and the
    boundary cannot drift apart from each other or from the run."""

    expected = claim_.context.format(*claim_.derive())
    flat = " ".join(CLAIM_BOUNDARY.read_text(encoding="utf-8").split())
    assert expected in flat, (claim_.source, expected)
