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

from causal_spacetime_lab.density_coupled_clocks import (
    bracket_width_against_worldline,
    chain_is_causal,
    harvest_chain_from_sprinkling_1p1,
    make_poisson_clock_chain_1p1,
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
        "arms": {},
    }
    verdicts: dict[str, bool] = {}

    for arm in ("thinned", "harvest_fixed", "harvest_scaled"):
        rows = run_arm(arm)
        lam_slope = fit_exponent(rows, "lam_mean")
        rmse_slope = fit_exponent(rows, "rmse")
        print(f"\n[{arm}] lam(rho) exponent = {lam_slope:+.3f}, "
              f"rmse(rho) exponent = {rmse_slope:+.3f}")
        for row in rows:
            extra = (
                f"  predicted={row['rmse_predicted']:.5f}"
                if "rmse_predicted" in row else ""
            )
            print(f"  rho={row['rho']:6d}  lam={row['lam_mean']:8.2f}  "
                  f"rmse={row['rmse']:.5f}{extra}  "
                  f"unreachable={row['unreachable']}")
        results["arms"][arm] = {
            "rows": rows,
            "lam_exponent": lam_slope,
            "rmse_exponent": rmse_slope,
        }

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
        else:
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

    results["verdicts"] = verdicts
    all_passed = all(verdicts.values())
    results["all_passed"] = bool(all_passed)
    results["status"] = (
        "theory-track characterization; no gate; no frozen artifact; "
        "harvested-chain fluctuation class left open"
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
