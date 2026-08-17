"""S6 mass-ladder foundation contract tests.

The exact-path rung constants (KAPPA, dt, SCALE) pinned against the
PI's independent check values; the central rung reproduced EXACTLY
(dt == 8.5, SCALE bit-exact vs the untouched o4_sizing); every
mass-genericity lemma condition certified per rung; and M = 1
fallback DETECTORS: the geometry path must not read the module
default mass, and the wrapper-contract path must deliver the rung
mass to every solver call. No seeds, no long runs, no results."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_sizing as sz  # noqa: E402
import o4b_g3a as g3a  # noqa: E402
import s6_rungs as s6  # noqa: E402

# ------------------------------------------------ the exact path

def test_kappa_and_the_central_rung_reproduce_exactly():
    assert s6.KAPPA == 1.224782580680992
    assert s6.dt(1.0) == 8.5                      # EXACT, ratio form
    assert s6.t_min(1.0) == 6.0 + 2.0 * math.log(16.0 / 10.0)
    g = s6.rung_geometry(1.0)
    assert g["scale"] == sz.SCALE                 # bit-exact
    assert g["r_lo"] == sz.R_LO and g["r_hi"] == sz.R_HI
    assert g["psi_max"] == sz.PSI_MAX
    assert g["l_max_ub"] == sz.L_MAX_UB


def test_the_frozen_ladder_constants_are_the_derived_ones():
    """The PI's independent check values, at full precision, and the
    import-time drift gate that guards them."""

    assert s6.LADDER == (1.0, 1.4, 1.8)
    assert s6.RUNG_CONSTANTS[1.4]["dt"] == 9.070565190742672
    assert s6.RUNG_CONSTANTS[1.4]["scale"] == 18141.08004658374
    assert s6.RUNG_CONSTANTS[1.8]["dt"] == 9.725248174609407
    assert s6.RUNG_CONSTANTS[1.8]["scale"] == 13679.093767152488
    for m, want in s6.RUNG_CONSTANTS.items():
        g = s6.rung_geometry(m)
        assert g["dt"] == want["dt"] and g["scale"] == want["scale"]
        assert g["mu"] == want["mu"] == 2.0 * m / 15.0


def test_mu_is_the_pre_frozen_indicator_and_genuinely_varies():
    mus = [s6.mu(m) for m in s6.LADDER]
    assert mus == sorted(mus)
    assert mus[-1] / mus[0] == pytest.approx(1.8, abs=1e-12)
    # no two rungs are isometric copies: same shell, different mu
    assert len({round(x, 12) for x in mus}) == 3


def test_every_lemma_condition_certifies_on_every_rung():
    for m in s6.LADDER:
        tab = s6.lemma_table(m)
        assert tab["all_pass"], (m, tab["rows"])
        for name, row in tab["rows"].items():
            assert row["certified"], (m, name)
            assert row["margin"] > 0.0, (m, name)


def test_the_winding_margin_grows_with_mass():
    """The L5 cross-check gets SAFER deeper: w_glob grows with M
    while psi_max shrinks."""

    margins = [s6.lemma_table(m)["rows"]["l5_winding"]["margin"]
               for m in s6.LADDER]
    assert margins == sorted(margins)
    assert margins[0] > 3.0


def test_out_of_domain_masses_are_refused():
    with pytest.raises(ValueError, match="outside"):
        s6.t_min(0.0)
    with pytest.raises(ValueError, match="outside"):
        s6.t_min(6.0)


# ------------------------------------------------ fallback detectors

def test_geometry_never_reads_the_module_default_mass(monkeypatch):
    """Poison cft.M_DEFAULT: the parametric path must be unaffected.
    A silent M = 1 fallback anywhere in the geometry chain would
    change the derived constants and fail here."""

    import certified_flight_time as cft

    want = s6.rung_geometry(1.4)
    monkeypatch.setattr(cft, "M_DEFAULT", 99.0)
    got = s6.rung_geometry(1.4)
    assert got == want


def test_the_wrapper_contract_delivers_the_rung_mass(monkeypatch):
    """Every solver call inside the G3a table machinery must carry
    the EXPLICIT rung mass -- never s1.M (post-ruling: several
    runners were M = 1-fixed even though the solver takes m)."""

    import s1_schwarzschild_cost as s1

    seen = {"m": set()}
    real_ft = s1.flight_time
    real_cr = s1.causal_relation

    def spy_ft(r1, r2, dpsi, m, tol, details=None):
        seen["m"].add(m)
        return real_ft(r1, r2, dpsi, m, tol, details)

    def spy_cr(p, q, m, tol):
        seen["m"].add(m)
        return real_cr(p, q, m, tol)

    monkeypatch.setattr(g3a.s1, "flight_time", spy_ft)
    monkeypatch.setattr(g3a.s1, "causal_relation", spy_cr)
    case = g3a.check_case("probe/radial", 13.0, 17.0, 0.0,
                          s1.DEFAULT_TOL, m=1.4)
    assert seen["m"] == {1.4}
    assert case["rows"]                            # the case really ran


def test_geometry_radii_accepts_a_rung_box():
    """The shell-edge cases must resolve to the RUNG's own sampling
    box, not the frozen M = 1 edges (review PR #90 R2); the anchors
    stay absolute either way."""

    g = s6.rung_geometry(1.8)
    box = (g["r_lo"], g["r_hi"])
    assert g3a.geometry_radii("x", "R_LO", "R_HI", box) == box
    assert g3a.geometry_radii("x", "R_LO", "R_HI") == (
        sz.R_LO, sz.R_HI)
    assert g3a.geometry_radii("anchors", None, None, box) == (
        sz.R_IN, sz.R_OUT)


def test_the_full_wrapper_table_runs_on_the_rung_geometry(
        monkeypatch):
    """The REAL G3a table at M = 1.8 with that rung's own box and
    psi_max: the angle resolver must receive the rung cap, the
    shell-edge cases the rung edges, and the whole contract must
    PASS at the deep rung's actual boundaries -- the case the review
    named (a wrapper preflight passing on M = 1 geometry says nothing
    about the rung's own edges)."""

    g = s6.rung_geometry(1.8)
    import o4_g3_redesign as g3

    seen = {"caps": set(), "radii": set()}
    real_resolve = g3.resolve_theta

    def spy_resolve(spec, r1, r2, m, tol, cap):
        seen["caps"].add(cap)
        seen["radii"].update((r1, r2))
        assert m == 1.8
        return real_resolve(spec, r1, r2, m, tol, cap)

    monkeypatch.setattr(g3a.g3, "resolve_theta", spy_resolve)
    out = g3a.run_preflight(m=1.8, geometry=g)
    assert out["passed"], out["failed_conditions"]
    assert seen["caps"] == {g["psi_max"]}
    assert g["r_lo"] in seen["radii"] and g["r_hi"] in seen["radii"]
    assert sz.R_LO not in seen["radii"]        # no M = 1 edge leaked


def test_the_wrapper_contract_defaults_stay_the_frozen_m1(monkeypatch):
    """The default path is byte-compatible with the executed
    campaigns: no explicit m means s1.M, and ONLY s1.M."""

    import s1_schwarzschild_cost as s1

    seen = set()
    real_ft = s1.flight_time

    def spy_ft(r1, r2, dpsi, m, tol, details=None):
        seen.add(m)
        return real_ft(r1, r2, dpsi, m, tol, details)

    monkeypatch.setattr(g3a.s1, "flight_time", spy_ft)
    out = g3a.solver_determinism(s1.DEFAULT_TOL)
    assert out["bit_identical"] is True
    assert seen == {s1.M}


# ------------------------------------------------ the addendum

def test_the_addendum_freezes_the_ladder_discipline():
    doc = (_REPO / "docs" / "theory"
           / "schwarzschild_volume_oracle_certification.md"
           ).read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "S6 addendum: mass-general instantiation" in flat
    assert "KAPPA = 8.5 / T_min(1)" in flat
    assert "1.224782580680992" in flat
    assert "MU = 2M/r_c" in flat
    assert "never an isometric copy" in flat
    assert "never copied from another rung" in flat
    assert "M = 0 is out of scope" in flat
    # the frozen table rows quote the exact constants
    assert "18141.08004658374" in flat
    assert "13679.093767152488" in flat
