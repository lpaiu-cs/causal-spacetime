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
selector, a construction setting) are NOT skipped: they are listed in
`ACCEPTED_EXCLUSIONS` with a one-line reason, and
`test_every_number_is_derived_or_explicitly_excluded` fails if a
number in Sections 4-6.7 is neither derived above nor listed there.
That guard is the actual promise -- the claim list alone could always
be out-run by a sentence added later.
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
       "(mean chain proper-time error {})")
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


def _weighted(key: str, places: int = 3) -> tuple[str, str]:
    rows = sorted((r for r in _csv("weighted_conformal_volume_summary.csv")
                   if r["profile"] == "constant_1.5"),
                  key=lambda r: float(r["N"]))
    assert len(rows) == 3, len(rows)
    return (f"{float(rows[0][key]):.{places}f}",
            f"{float(rows[-1][key]):.{places}f}")


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "decreasing relative RMSE ({} -> {} across the tested N)")
def _weighted_rmse() -> tuple[str, ...]:
    return _weighted("weighted_relative_rmse")


@claim("4.6", "exp19 weighted_conformal_volume_summary.csv",
       "no material RMSE improvement ({} -> {}, bias ~ {})")
def _unweighted() -> tuple[str, ...]:
    rows = [r for r in _csv("weighted_conformal_volume_summary.csv")
            if r["profile"] == "constant_1.5"]
    bias = statistics.fmean(_col(rows, "unweighted_volume_bias"))
    return (*_weighted("unweighted_relative_rmse"), f"{bias:.2f}")


@claim("4.6", "exp20 conformal_volume_exact_sanity.csv",
       "the analytic volume/proper-time formulas are verified to ~{}")
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


# ================== 4.7 Horizon analogue -- Rindler inaccessibility

@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "precision = recall = {}, zero false positives and false negatives")
def _rindler_exact() -> tuple[str, ...]:
    rows = _csv("rindler_horizon_reconstruction_summary.csv")
    values = set(_col(rows, "precision")) | set(_col(rows, "recall"))
    assert len(values) == 1, values
    return (f"{values.pop():.1f}",)


@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "(accessible fraction ~{})")
def _rindler_fraction() -> tuple[str, ...]:
    rows = _csv("rindler_horizon_reconstruction_summary.csv")
    return (f'{statistics.fmean(_col(rows, "wedge_accessible_fraction")):.2f}',)


@claim("4.7", "exp16 rindler_horizon_reconstruction_summary.csv",
       "radar-time RMSE falling with tick resolution (e.g. {} -> {})")
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
       "return the identical single-observer distance {}, while the two-chain "
       "oriented protocol recovers the signed positions {} and {}")
def _reflection() -> tuple[str, ...]:
    rows = [r for r in _csv("single_observer_reflection_degeneracy.csv")
            if abs(float(r["true_x"])) == 0.1]
    assert len(rows) == 2, len(rows)
    pair = sorted(rows, key=lambda r: -float(r["true_x"]))
    distances = {f'{float(r["single_observer_radar_distance"]):.1f}' for r in pair}
    assert len(distances) == 1, distances
    return (distances.pop(),
            *(f'{float(r["two_chain_signed_position"]):+.1f}' for r in pair))


@claim("5", "exp05 finite_speed_lattice_growth.csv",
       "while finite-t counts differ (t = 5: {} vs {}; t = 30: {} vs {})")
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

@claim("6.6", "p14_s1_cost.json",
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
       "the Schwarzschild discrimination is strong but imperfect (AUC ≈ {};")
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

    # 4.6: "negligible observed bias" against "a large low bias"
    weighted = [r for r in _csv("weighted_conformal_volume_summary.csv")
                if r["profile"] == "constant_1.5"]
    assert max(abs(float(r["weighted_volume_bias"])) for r in weighted) < 0.01
    assert all(float(r["unweighted_volume_bias"]) < -0.25 for r in weighted)

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

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    for sentence in (*(plane[k]["sentence"] for k in ("c1", "c2")),
                     *s4["sentences"], *s5["sentences"]):
        assert sentence in manuscript, sentence[:40]

    section = _section("6.7")
    assert "Stage verdict: CONFIRMED." in section
    assert "zero ambiguous pairs and one escalated pair" in section
    assert "Zero ambiguous and zero escalated pairs" in section


# ----------------------------------------------- the completeness guard

#: Numbers printed in Sections 4-6.7 that no committed artifact
#: produces as a RESULT. Each entry is the literal as printed, mapped
#: to why it is not derived. This dict is the honest edge of the
#: contract: everything not listed here is re-derived above, and the
#: guard below fails if a number is neither.
ACCEPTED_EXCLUSIONS: dict[str, str] = {
    "0": "the A = 0 control reading and 'zero' counts, construction settings",
    "1": "rung/claim labels (R1, C1, S1, M = 1) and unit-scale prose",
    "2": "rung labels and the dimension-2 series label, not measurements",
    "3": "rung labels and the dimension-3 series label, not measurements",
    "4": "rung labels and the dimension-4 series label, not measurements",
    "5": "rung label R5 and the t = 5 row selector of the exp05 claim",
    "8": "the S4 manifest file count, derived inside its own claim",
    "10": "the r-shell floor and 'about 10%', derived inside their claims",
    "16": "exp11 tick-count selector, quoted inside the radar-RMSE claim",
    "20": "the r-shell ceiling, derived inside the S4 domain claim",
    "25": "the 25% retention selector of the exp23 uncorrected-RMSE claim",
    "30": "the t = 30 row selector of the exp05 lattice-count claim",
    "32": "tick-count selector, quoted inside the exp13/exp14 claims",
    "40": "the coordinate-time extent, derived inside the S4 domain claim",
    "95": "the confidence level of every stated interval, a frozen convention",
    "128": "tick-count selector, quoted inside the exp13/exp14 claims",
    "200": "exp03 N selector, quoted inside the chain-convergence claim",
    "300": "N selector for exp10/exp11/exp13 and the P14 operating point",
    "600": "exp10 N selector, quoted inside the dimension-track claim",
    "1200": "exp10 N selector, quoted inside the dimension-track claim",
    "2000": "exp03 N selector, quoted inside the chain-convergence claim",
    "2400": "exp10/exp07 N selector, quoted inside the RMSE claims",
    "0.1": "the +-0.1 exp12 target positions, derived inside their claim",
    "0.3": "exp13 beta selector, quoted inside the Lorentz-map claim",
    "0.6": "exp13 beta selector, quoted inside the Lorentz-map claim",
    "1.0": "unit scale/weight/fraction labels and the wT = 1 selector",
    "1.5": "the conformal scale label of the exp18/exp19 profiles",
    "2.0": "the wT = 2 selector and the conformal scale label",
    "0.25": "'about a quarter' prose beside the derived ~0.25 fraction",
    "-1": "det g = -1, an exact property of the Brinkmann construction",
}


_CROSS_REFERENCE = re.compile(
    r"Sections? \d+(?:\.\d+)?(?:[-–]\d+(?:\.\d+)?)?"   # Section 4.6, Sections 4-6.8
    r"|Figures? \d+|Appendix [A-Z]|Table \d+"               # figure/appendix pointers
    r"|\bP\d+(?:[-–]P\d+)?\b|\bexp\d+\b"               # programme and producer ids
    r"|\b\d\+\dD\b|\b\dD\b"                                 # 1+1D, 4D
    r"|\([a-d]\)"                                           # figure panel letters
)

_NUMBER = re.compile(r"(?<![\w.])[-+]?\d(?:[\d,]*\d)?(?:\.\d+)?(?:[eE][-+]?\d+)?")

_SCANNED = ("4", "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "5",
            "6", "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7")


def _frozen_sentences() -> tuple[str, ...]:
    """The frozen verdict sentences, quoted in Korean byte-for-byte.
    Their numerals are the RULE's, not the paper's, and the quotes are
    already pinned against their artifacts verbatim by
    `test_the_preregistered_verdicts_are_the_artifacts_verdicts` --
    a stricter check than a numeral scan. Read from the artifacts, so
    a re-frozen sentence cannot leave a stale hole here."""

    plane = _json("p14_prereg_results.json")
    return (*(plane[k]["sentence"] for k in ("c1", "c2")),
            *_json("p14_s4_results.json")["sentences"],
            *_json("p14_s5_results.json")["sentences"])


def _scannable(number: str) -> str:
    body = SECTIONS.get(number, "")
    for sentence in _frozen_sentences():
        body = body.replace(sentence, " ")
    return _CROSS_REFERENCE.sub(" ", body)


def _claimed_literals() -> set[str]:
    """Every number that a claim above actually derived and matched."""

    found: set[str] = set()
    for claim_ in CLAIMS:
        found.update(_NUMBER.findall(claim_.context.format(*claim_.derive())))
    return found


def test_every_number_is_derived_or_explicitly_excluded():
    """The guard the [P2] review asked for. Every numeral printed in
    Sections 4-6.7 must be one a claim above re-derived from an
    artifact, or one named in ACCEPTED_EXCLUSIONS with its reason. A
    number added to the manuscript later lands here rather than in a
    silent gap."""

    claimed = _claimed_literals()
    uncovered: dict[str, str] = {}
    for number in _SCANNED:
        for token in _NUMBER.findall(_scannable(number)):
            if token not in claimed and token not in ACCEPTED_EXCLUSIONS:
                uncovered[token] = number
    assert not uncovered, uncovered


def test_no_exclusion_is_stale():
    """An exclusion matching nothing in Sections 4-6.7 is a standing
    licence for a number that is no longer there; drop it instead."""

    present: set[str] = set()
    for number in _SCANNED:
        present.update(_NUMBER.findall(_scannable(number)))
    assert not (set(ACCEPTED_EXCLUSIONS) - present), \
        sorted(set(ACCEPTED_EXCLUSIONS) - present)
    assert all(reason.strip() for reason in ACCEPTED_EXCLUSIONS.values())


def test_the_claim_table_covers_every_section_and_legacy_table():
    """Guard the guard: a claim table that quietly lost a section, or a
    committed legacy table that no claim reads, would leave exactly the
    hole this file exists to close."""

    covered = {c.section for c in CLAIMS}
    assert covered == {"4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7",
                       "5", "6.2", "6.3", "6.4", "6.6", "6.7"}, sorted(covered)

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
