"""G2 instrumentation audit + density-scaling characterization (T1 v0.5).

Question. The T1 roadmap's ``error ~ rho^{-1/2}`` law presumes tick
statistics coupled to the sprinkling density as ``lam ~ rho * ell``.
With the density-coupled protocols now built
(``causal_spacetime_lab.density_coupled_clocks``), does that law hold --
and for WHICH protocol?

Three arms, same estimator (``d_hat = (W - 1) / (2 lam)``):

- ``thinned``: Poisson clock at exact x0 with lam = rho * ell (the
  coupling the document presumes). Here ``W - 1 ~ Poisson(2 lam d)`` by
  the Lemma 2 identity, so ``sd(d_hat) = sqrt(d / (2 lam))`` is a PROVED
  corollary of Theorem 2's setup, and RMSE ~ rho^{-1/2} follows. This arm
  is asserted against its exact prediction, not just eyeballed.
- ``harvest_fixed``: longest sprinkled chain in a fixed-width tube. The
  chain rate is a *measurement* (expected ~ sqrt(rho), the discreteness
  scale); the fixed tube width also sets a position-wiggle floor, so the
  error exponent is expected to flatten. Characterization only.
- ``harvest_scaled``: tube width ~ 3 / sqrt(rho) (a few discreteness
  lengths), removing the fixed floor. Characterization only: the chain's
  fluctuation class is NOT Poisson (a maximal path concentrates), so no
  exponent is asserted beyond a sanity band -- distinguishing the
  Poisson-rate guess (-1/4, from lam ~ sqrt(rho)) from a KPZ-like -1/3
  needs a dedicated study and is recorded as an open question.
- ``harvest_order_only``: longest chain between two DESIGNATED anchor
  events (nearest sprinkled events to the window endpoints on the
  observer line -- a coordinate-assisted DESIGNATION, made once, like
  placing observers). The SELECTION rule then reads order data alone,
  answering the order-only design question the tube protocol left open.
  The chain is pinned at the anchors but otherwise free to wander
  transversally; the arm therefore also measures the wandering, as the
  pooled transverse RMS of tick positions about the anchor line.

Pre-stated expectations for the order-only arm (written before the
grid was run; recorded as outcomes in the results table, and NOT gates
except where marked):

- E1 (gate, sanity band): the chain rate couples at the discreteness
  scale, lam-exponent in (0.35, 0.70), as for the tube arms.
- E2 (gate, same form as the scaled tube): the RMSE law is distinctly
  shallower than the thinned clock's rho^{-1/2}.
- E3 (recorded, directional): with the tube constraint removed, the
  chain wanders. If longest-chain wandering is KPZ-class (wandering
  exponent 2/3), the transverse RMS scales like rho^{-1/6} and, where
  it dominates, drags the RMSE exponent shallower than the scaled
  tube's; a diffusive alternative predicts rho^{-1/4}. The outcome
  reports the fitted transverse exponent and which candidate is
  nearer.

Deliberately NOT done here: no frozen gate, no confirmatory claim, no
pooling with any preregistered artifact. This is theory-track
instrumentation and its characterization.

Usage:
    python t1_g2_density_scaling.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.density_coupled_clocks import (
    bracket_width_against_worldline,
    chain_is_causal,
    harvest_chain_from_sprinkling_1p1,
    harvest_order_only_chain_1p1,
    make_poisson_clock_chain_1p1,
    nearest_event_index,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

ROOT = Path(__file__).resolve().parents[2]
# The full-grid table is a TRACKED artifact (committed evidence for the
# doc's execution-outcome numbers), so the writer targets the committed
# path directly: a rerun after any parameter or code change shows up as
# a git diff on the cited file instead of going stale beside a fresh
# copy in the ignored outputs/ tree.
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_g2_density_scaling_results.json"

DIAMOND_T = 2.0
DIAMOND_AREA = DIAMOND_T * DIAMOND_T / 2.0
OBSERVERS = (-0.3, 0.3)
TICK_WINDOW = (-0.6, 0.6)
ELL = 0.1
TUBE_SCALE = 3.0
RHO_GRID = (500, 1000, 2000, 4000, 8000, 16000, 32000)
SEEDS = 16
TARGETS_PER_SCENE = 30


def audit_harvested_chains(seed: int = 11, scenes: int = 5) -> dict:
    """Constructor audit: chain property, containment, determinism."""

    rng = np.random.default_rng(seed)
    all_causal = True
    all_contained = True
    all_simple = True
    deterministic = True
    for _ in range(scenes):
        rho = float(rng.choice([1000, 4000, 16000]))
        bulk = sprinkle_1p1_causal_diamond(
            rng.poisson(rho * DIAMOND_AREA), T=DIAMOND_T, seed=rng
        )
        x0 = float(rng.choice(OBSERVERS))
        width = TUBE_SCALE / np.sqrt(rho)
        idx = harvest_chain_from_sprinkling_1p1(
            bulk, x0, width, TICK_WINDOW[0], TICK_WINDOW[1]
        )
        idx_again = harvest_chain_from_sprinkling_1p1(
            bulk, x0, width, TICK_WINDOW[0], TICK_WINDOW[1]
        )
        deterministic = deterministic and np.array_equal(idx, idx_again)
        chain = bulk[idx]
        all_causal = all_causal and chain_is_causal(chain)
        all_contained = all_contained and bool(
            np.all(np.abs(chain[:, 1] - x0) <= width / 2.0)
            and np.all(chain[:, 0] >= TICK_WINDOW[0])
            and np.all(chain[:, 0] <= TICK_WINDOW[1])
        )
        all_simple = all_simple and bool(
            np.all(np.diff(chain[:, 0]) > 0)
        )
    return {
        "chains_causal": all_causal,
        "chains_contained": all_contained,
        "tick_times_strictly_increasing": all_simple,
        "harvest_deterministic": deterministic,
        "passed": bool(
            all_causal and all_contained and all_simple and deterministic
        ),
    }


def audit_order_only_chains(seed: int = 13, scenes: int = 5) -> dict:
    """Constructor audit for the order-only harvest: chain property,
    anchor endpoints, interval containment, determinism, REFLECTION
    invariance (a spatial reflection preserves the labelled order and
    must preserve the harvested chain exactly -- the tie-break reads
    labels, not coordinates), and -- at low density, where the matrix
    is affordable -- longest-chain LENGTH equality against an
    independent causal-matrix DP (the patience-sorting implementation
    must reproduce the order-theoretic optimum, not just some
    chain)."""

    rng = np.random.default_rng(seed)
    all_causal = True
    endpoints_ok = True
    contained = True
    all_simple = True
    deterministic = True
    reflection_invariant = True
    dp_lengths_match = True
    dp_checks = 0
    for _ in range(scenes):
        rho = float(rng.choice([500, 1000, 2000]))
        bulk = sprinkle_1p1_causal_diamond(
            rng.poisson(rho * DIAMOND_AREA), T=DIAMOND_T, seed=rng
        )
        x0 = float(rng.choice(OBSERVERS))
        bottom = nearest_event_index(bulk, TICK_WINDOW[0], x0)
        top = nearest_event_index(bulk, TICK_WINDOW[1], x0)
        idx = harvest_order_only_chain_1p1(bulk, bottom, top)
        idx_again = harvest_order_only_chain_1p1(bulk, bottom, top)
        deterministic = deterministic and np.array_equal(idx, idx_again)
        reflected = bulk.copy()
        reflected[:, 1] = -reflected[:, 1]
        idx_reflected = harvest_order_only_chain_1p1(reflected, bottom, top)
        reflection_invariant = reflection_invariant and np.array_equal(
            idx, idx_reflected
        )
        chain = bulk[idx]
        all_causal = all_causal and chain_is_causal(chain)
        endpoints_ok = endpoints_ok and bool(
            idx[0] == bottom and idx[-1] == top
        )
        u = bulk[:, 0] + bulk[:, 1]
        v = bulk[:, 0] - bulk[:, 1]
        interior = idx[1:-1]
        contained = contained and bool(
            np.all(u[interior] > u[bottom]) and np.all(v[interior] > v[bottom])
            and np.all(u[interior] < u[top]) and np.all(v[interior] < v[top])
        )
        all_simple = all_simple and bool(np.all(np.diff(chain[:, 0]) > 0))
        inside = (u > u[bottom]) & (v > v[bottom]) & (u < u[top]) & (v < v[top])
        sub = bulk[np.flatnonzero(inside)]
        sub = sub[np.argsort(sub[:, 0], kind="stable")]
        causal = causal_matrix_1p1(sub)
        best = np.ones(sub.shape[0], dtype=int)
        for j in range(sub.shape[0]):
            preds = np.flatnonzero(causal[:, j])
            if preds.size:
                best[j] = best[preds].max() + 1
        dp_len = (int(best.max()) if sub.shape[0] else 0) + 2
        dp_lengths_match = dp_lengths_match and dp_len == idx.size
        dp_checks += 1
    return {
        "chains_causal": all_causal,
        "anchor_endpoints": endpoints_ok,
        "interval_contained": contained,
        "tick_times_strictly_increasing": all_simple,
        "harvest_deterministic": deterministic,
        "reflection_invariant": reflection_invariant,
        "dp_length_checks": dp_checks,
        "dp_lengths_match": dp_lengths_match,
        "passed": bool(
            all_causal and endpoints_ok and contained and all_simple
            and deterministic and reflection_invariant and dp_lengths_match
        ),
    }


def _scene_targets(rng: np.random.Generator) -> np.ndarray:
    return np.column_stack([
        rng.uniform(-0.05, 0.05, size=TARGETS_PER_SCENE),
        rng.uniform(-0.15, 0.15, size=TARGETS_PER_SCENE),
    ])


def run_arm(arm: str, rho_grid=RHO_GRID, seeds=SEEDS) -> list[dict]:
    """Measure (lam, rmse, predicted rmse where proved) per density."""

    rows = []
    for rho in rho_grid:
        lam_values: list[float] = []
        errors: list[float] = []
        distances: list[float] = []
        transverse_sq: list[float] = []
        unreachable = 0
        short_clocks = 0
        for s in range(seeds):
            rng = np.random.default_rng(100_000 * s + int(rho))
            targets = _scene_targets(rng)
            bulk = sprinkle_1p1_causal_diamond(
                rng.poisson(rho * DIAMOND_AREA), T=DIAMOND_T, seed=rng
            )
            for x0 in OBSERVERS:
                if arm == "thinned":
                    lam_true = rho * ELL
                    ticks = make_poisson_clock_chain_1p1(
                        TICK_WINDOW[0], TICK_WINDOW[1], lam_true, x0, seed=rng
                    )
                    lam_decode = lam_true  # known by construction (proved arm)
                elif arm == "harvest_order_only":
                    bottom = nearest_event_index(bulk, TICK_WINDOW[0], x0)
                    top = nearest_event_index(bulk, TICK_WINDOW[1], x0)
                    try:
                        idx = harvest_order_only_chain_1p1(bulk, bottom, top)
                    except ValueError:
                        # anchors not causally related: a clock FAILURE,
                        # counted like a short clock, never silently skipped
                        unreachable += TARGETS_PER_SCENE
                        short_clocks += 1
                        continue
                    ticks = bulk[idx]
                else:
                    width = (
                        ELL if arm == "harvest_fixed"
                        else TUBE_SCALE / np.sqrt(rho)
                    )
                    idx = harvest_chain_from_sprinkling_1p1(
                        bulk, x0, width, TICK_WINDOW[0], TICK_WINDOW[1]
                    )
                    ticks = bulk[idx]
                if len(ticks) < 4:
                    # A clock too short to measure is a clock FAILURE:
                    # its targets count as unreachable so that sparse
                    # settings or a harvester regression cannot hide by
                    # shrinking the denominator of the surviving clocks.
                    unreachable += TARGETS_PER_SCENE
                    short_clocks += 1
                    continue
                if arm != "thinned":
                    # no true rate exists: the decoder estimates it
                    lam_decode = (len(ticks) - 1) / (
                        ticks[-1, 0] - ticks[0, 0]
                    )
                    transverse_sq.extend(
                        np.square(ticks[:, 1] - x0).tolist()
                    )
                lam_values.append(
                    (len(ticks) - 1) / (ticks[-1, 0] - ticks[0, 0])
                )
                for t, x in targets:
                    w = bracket_width_against_worldline(ticks, float(t), float(x))
                    if np.isnan(w):
                        unreachable += 1
                        continue
                    d_true = abs(x - x0)
                    errors.append(abs((w - 1) / (2 * lam_decode) - d_true))
                    distances.append(d_true)
        rmse = float(np.sqrt(np.mean(np.square(errors))))
        row = {
            "rho": rho,
            "lam_mean": float(np.nanmean(lam_values)),
            "rmse": rmse,
            "n_measurements": len(errors),
            "unreachable": unreachable,
            "short_clocks": short_clocks,
        }
        if arm == "thinned":
            row["rmse_predicted"] = float(
                np.sqrt(np.mean(np.asarray(distances) / (2 * rho * ELL)))
            )
        else:
            # pooled transverse RMS of tick positions about the nominal
            # worldline: the wandering observable of expectation E3
            row["transverse_rms"] = float(
                np.sqrt(np.mean(transverse_sq))
            )
        rows.append(row)
    return rows


def fit_exponent(rows: list[dict], key: str) -> float:
    log_rho = np.log([row["rho"] for row in rows])
    log_val = np.log([row[key] for row in rows])
    return float(np.polyfit(log_rho, log_val, 1)[0])


def main() -> None:

    audit = audit_harvested_chains()
    print("[AUDIT]", json.dumps(audit))
    if not audit["passed"]:
        raise SystemExit("harvested-chain audit failed; no scaling is run")
    audit_order_only = audit_order_only_chains()
    print("[AUDIT order-only]", json.dumps(audit_order_only))
    if not audit_order_only["passed"]:
        raise SystemExit("order-only-chain audit failed; no scaling is run")

    results: dict = {
        "config": {
            "rho_grid": list(RHO_GRID),
            "seeds_per_density": SEEDS,
            "targets_per_scene": TARGETS_PER_SCENE,
            "ell": ELL,
            "tube_scale": TUBE_SCALE,
            "observers": list(OBSERVERS),
            "tick_window": list(TICK_WINDOW),
            "diamond_T": DIAMOND_T,
        },
        "audit": audit,
        "audit_order_only": audit_order_only,
        "arms": {},
    }
    verdicts: dict[str, bool] = {}

    for arm in (
        "thinned", "harvest_fixed", "harvest_scaled", "harvest_order_only"
    ):
        rows = run_arm(arm)
        lam_slope = fit_exponent(rows, "lam_mean")
        rmse_slope = fit_exponent(rows, "rmse")
        print(f"\n[{arm}] lam(rho) exponent = {lam_slope:+.3f}, "
              f"rmse(rho) exponent = {rmse_slope:+.3f}")
        for row in rows:
            if "rmse_predicted" in row:
                extra = f"  predicted={row['rmse_predicted']:.5f}"
            elif "transverse_rms" in row:
                extra = f"  trans_rms={row['transverse_rms']:.5f}"
            else:
                extra = ""
            print(f"  rho={row['rho']:6d}  lam={row['lam_mean']:8.2f}  "
                  f"rmse={row['rmse']:.5f}{extra}  "
                  f"unreachable={row['unreachable']}")
        results["arms"][arm] = {
            "rows": rows,
            "lam_exponent": lam_slope,
            "rmse_exponent": rmse_slope,
        }

        # Clock failures are asserted, not just reported: a run with
        # short clocks or unreachable targets must not write
        # all_passed = true even if the broad exponent bands survive.
        verdicts[f"{arm}_no_unreachable"] = bool(
            all(row["unreachable"] == 0 for row in rows)
        )
        verdicts[f"{arm}_no_short_clocks"] = bool(
            all(row["short_clocks"] == 0 for row in rows)
        )

        if arm == "thinned":
            # PROVED corollary (Theorem 2 setup + lam = rho * ell):
            # RMSE must track sqrt(mean d / (2 lam)) and fall ~ rho^{-1/2}.
            # Tolerances account for intra-chain correlation: the ~30
            # targets on one chain share its tick realization (the shared
            # regions of Theorem 2 Step 1!), so the effective sample per
            # density is the number of independent chains (2 observers x
            # SEEDS), giving each per-density RMSE a relative sd of
            # roughly 1/sqrt(2 x 2 x SEEDS) ~ 12%. The grid-mean ratio
            # over 7 independent densities is ~5% -- asserted at 15%;
            # per-density ratios get a loose [0.6, 1.6] guard.
            ratios = [row["rmse"] / row["rmse_predicted"] for row in rows]
            verdicts["thinned_lam_exponent_near_1"] = bool(
                abs(lam_slope - 1.0) < 0.10
            )
            verdicts["thinned_mean_prediction_ratio_near_1"] = bool(
                abs(float(np.mean(ratios)) - 1.0) < 0.15
            )
            verdicts["thinned_per_density_ratios_sane"] = bool(
                all(0.6 < ratio < 1.6 for ratio in ratios)
            )
            verdicts["thinned_rmse_exponent_near_minus_half"] = bool(
                -0.60 < rmse_slope < -0.40
            )
        elif arm == "harvest_fixed":
            verdicts["harvest_lam_exponent_near_half"] = bool(
                0.35 < lam_slope < 0.70
            )
        elif arm == "harvest_scaled":
            verdicts["harvest_scaled_lam_exponent_near_half"] = bool(
                0.35 < fit_exponent(rows, "lam_mean") < 0.70
            )
            # Characterization sanity only: distinctly shallower than -1/2,
            # in the band spanned by the Poisson-rate guess (-1/4) and a
            # KPZ-like -1/3. The sharp exponent is an OPEN question.
            verdicts["harvest_scaled_rmse_shallower_than_thinned"] = bool(
                rmse_slope > results["arms"]["thinned"]["rmse_exponent"] + 0.10
            )
            verdicts["harvest_scaled_rmse_in_sanity_band"] = bool(
                -0.45 < rmse_slope < -0.10
            )
        else:
            # Order-only arm gates (E1, E2 of the header): discreteness-
            # scale rate coupling; RMSE distinctly shallower than the
            # thinned clock's proved rho^{-1/2}; plus a wide instrument-
            # sanity bound. The wandering class itself (E3) is recorded,
            # not gated -- see order_only_recorded_expectations.
            verdicts["order_only_lam_exponent_near_half"] = bool(
                0.35 < lam_slope < 0.70
            )
            verdicts["order_only_rmse_shallower_than_thinned"] = bool(
                rmse_slope > results["arms"]["thinned"]["rmse_exponent"] + 0.10
            )
            verdicts["order_only_rmse_in_wide_sanity_bound"] = bool(
                -0.50 < rmse_slope < 0.02
            )

    order_rows = results["arms"]["harvest_order_only"]["rows"]
    transverse_slope = fit_exponent(order_rows, "transverse_rms")
    scaled_rmse = results["arms"]["harvest_scaled"]["rmse_exponent"]
    order_rmse = results["arms"]["harvest_order_only"]["rmse_exponent"]
    candidates = {"kpz_wandering_-1/6": -1.0 / 6.0, "diffusive_-1/4": -0.25}
    nearest = min(
        candidates, key=lambda name: abs(candidates[name] - transverse_slope)
    )
    results["order_only_recorded_expectations"] = {
        "gating": False,
        "stated_before_run": {
            "E3_wandering": (
                "KPZ-class wandering predicts transverse RMS ~ rho^{-1/6} "
                "and, where it dominates, an RMSE exponent shallower than "
                "the scaled tube's; diffusive predicts rho^{-1/4}"
            ),
        },
        "outcomes": {
            "transverse_rms_exponent": transverse_slope,
            "nearest_wandering_candidate": nearest,
            "rmse_exponent_order_only": order_rmse,
            "rmse_exponent_scaled_tube": scaled_rmse,
            "order_only_shallower_than_scaled_tube": bool(
                order_rmse > scaled_rmse
            ),
        },
    }
    print(
        "\n[order-only E3] transverse exponent "
        f"{transverse_slope:+.3f} (nearest: {nearest}); "
        f"rmse {order_rmse:+.3f} vs scaled tube {scaled_rmse:+.3f}"
    )

    results["verdicts"] = verdicts
    all_passed = all(verdicts.values())
    results["all_passed"] = bool(all_passed)
    results["status"] = (
        "theory-track characterization; no gate; no frozen artifact; "
        "order-only harvest measured (wandering exponent recorded); "
        "harvested-chain count-fluctuation class left open"
    )

    RESULTS_PATH.write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    print(f"\nverdicts: {json.dumps(verdicts, indent=2)}")
    print(f"all_passed = {all_passed}")
    print(f"wrote {RESULTS_PATH}")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
