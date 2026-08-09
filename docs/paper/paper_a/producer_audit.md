# Paper A legacy-table producer audit (read-only)

Status: **audit complete, 2026-08-09.** Baseline `325df55`; audited against
current `origin/main` (`36caa79`). This document MODIFIES NO PRODUCER — it is
the pure audit the integration plan's C-1 step called for, and the branch
decision it feeds: with the findings below, the debt closes with a **single
regeneration PR** whose one prerequisite is recording the exact environment,
because no remediation of any producer is required.

## Headline findings

1. **All 18 producer scripts (19 manifest rows; exp14 produces two tables)
   are byte-identical to the baseline** — `git diff --stat 325df55..HEAD`
   is empty for every `experiments/exp*.py` producer.
2. **Exactly one shared library dependency changed since baseline:**
   `src/causal_spacetime_lab/causal.py` (+11/−2, commit `4e5b16e`), which
   extracts the literal `1e-12` into `DEFAULT_CAUSAL_ATOL` and uses it as the
   default `atol` — **numerically identical, behavior-neutral for
   regeneration.** Every other imported `src/` module is unchanged (the rest
   of the src diff since baseline is new files not imported by any producer).
3. **No producer is UNSEEDED.** Every stochastic producer seeds every RNG
   entry point explicitly (details per producer below); exp12 and exp20 use
   no RNG at all. No module-level `np.random.*`, no `os.urandom`, no
   time-based seeding, no parallelism anywhere in the producer chain.
4. **The one realistic regeneration risk is environmental:** `pyproject.toml`
   pins floors only (`numpy>=2.0`, `matplotlib>=3.8`, `requires-python
   >=3.11`; no lock file, CI floats). NumPy freezes bit-generator streams but
   reserves the right to change `Generator` distribution-method streams
   (`uniform`/`normal`/`choice`/`integers`) across feature releases, and the
   original run's versions are unrecorded — so bit-for-bit agreement with the
   cited numbers is guaranteed only under the same numpy feature version.
5. Minor cross-platform caveat for exp16/exp17 only: Rindler tick math runs
   through libm transcendentals (`sinh/cosh/log`); near-boundary
   accessibility classifications could in principle flip at the 1-ulp level
   across OS/CPU (the causal predicate itself is cushioned by
   `atol = 1e-12`).

## Consequence for the regeneration PR

- Regenerate all 14 missing tables (and re-verify the 5 committed copies)
  in ONE pinned, recorded environment; record interpreter, numpy, and
  matplotlib versions plus the regeneration commit in `artifact_manifest.md`
  (its closing paragraph already demands this).
- Run from the repository root (all output paths are relative
  `outputs/...`); use each script's defaults — the cited tables correspond
  to default arguments in every case.
- If any regenerated number differs from a cited one (the numpy-version risk
  above), the manuscript number is corrected to the regenerated value with a
  changelog entry — a rerun under a recorded environment supersedes an
  unrecorded one; no producer edit is implied.

## Per-producer summary

Seed-derivation formulas are quoted from the scripts; "default cmd" means
`python experiments/<script>.py` with argparse defaults (or no argparse).
All producers depend only on numpy (+ matplotlib for figures, except exp15
and exp20 which are numpy-only). None import scipy or numba.

| # | Producer → table(s) | RNG / seeds | Rows | Verdict |
| --- | --- | --- | --- | --- |
| exp03 | `timelike_reconstruction_summary.csv` | `BASE_SEED = 20260424`; sprinkle seeds `BASE_SEED + offset` over 4 event counts | 4 | DETERMINISTIC-SEEDED |
| exp05 | `finite_speed_lattice_growth.csv` | `SEED = 20260424`; lattice arm fully deterministic | 31 | DETERMINISTIC-SEEDED |
| exp06 | `spacelike_distance_proxy_summary.csv` | `SEED = 20260424`, single `default_rng` threaded into sprinkle + pair choice | ≤350 (None-interval pairs skipped; seed-determined) | DETERMINISTIC-SEEDED |
| exp07 | `timelike_pair_reconstruction_summary.csv` (committed copy) | `--seed 0`; `seed + 10_000*n_index + rep` (sprinkle), `+ 20_000*n_index + rep` (pairs) | 4 | DETERMINISTIC-SEEDED |
| exp08 | `probe_pair_statistical_calibration_summary.csv` | `--seed 0`; support/probe streams split by `10_000`/`20_000` offsets | 4 | DETERMINISTIC-SEEDED |
| exp09 | `longest_chain_calibration_summary.csv` | `--seed 0`; sprinkle `seed + rep`, probes `seed + 10_000 + rep` | 60 (pair-level) | DETERMINISTIC-SEEDED |
| exp10 | `dimension_reconstruction_summary.csv` (committed copy) | `--seed 0`; `seed + 100_000*dim + 1_000*n + rep`; inversion deterministic bisection | 12 | DETERMINISTIC-SEEDED |
| exp11 | `discrete_radar_reconstruction_summary.csv` (committed copy) | `--seed 0`; sprinkle `seed + 10_000*n + rep`; chains analytic | 16 | DETERMINISTIC-SEEDED |
| exp12 | `single_observer_reflection_degeneracy.csv` | none — 8 hardcoded events, closed-form | 8 | DETERMINISTIC (no RNG) |
| exp13 | `oriented_radar_lorentz_summary.csv` (committed copy) | `--seed 0`; `seed + 100_000*n + 1_000*rep` (beta reuses support by design); beta fit = deterministic grid + argmin | 18 | DETERMINISTIC-SEEDED |
| exp14 | `observer_atlas_transition_summary.csv`, `observer_atlas_loop_summary.csv` | `--seed 0`; sprinkle `seed + 100_000*n + 1_000*rep`; invariant pairs `seed + 10_000*n + 100*rep + transition` | 135 / 45 | DETERMINISTIC-SEEDED |
| exp15 | `exact_poincare_map_sanity.csv` | hardcoded `seed=0` (not CLI-exposed); fits on grid `num_grid=1901` | 4 | DETERMINISTIC-SEEDED |
| exp16 | `rindler_horizon_reconstruction_summary.csv` | `--seed 0`; `seed + 100_000*accel + 10_000*n + rep`; Rindler math analytic | 18 | DETERMINISTIC-SEEDED (libm caveat) |
| exp17 | `inertial_vs_rindler_accessibility.csv` | `--seed 0`; single sprinkle `seed=config.seed` | 800 (event-level) | DETERMINISTIC-SEEDED (libm caveat) |
| exp18 | `conformal_order_ambiguity_summary.csv` | `--seed 0`; single sprinkle; conformal integrals deterministic `np.trapezoid` (`num_t=4096`) | 4 | DETERMINISTIC-SEEDED |
| exp19 | `weighted_conformal_volume_summary.csv` (committed copy) | `--seed 0`; support `seed + 10_000*n + rep`, probes `seed + 100_000*n + 1_000*rep` | 9 | DETERMINISTIC-SEEDED |
| exp20 | `conformal_volume_exact_sanity.csv` | none — pure quadrature | 3 | DETERMINISTIC (no RNG) |
| exp23 | `thinning_coarse_graining_summary.csv` | `--seed 0`; support `seed + rep`, probes `seed + 10_000 + rep`, thinning `seed + 100_000*rep + int(1000*keep_p)`; `keep_p = 1.0` bypasses RNG | 4 | DETERMINISTIC-SEEDED |

Committed figure-source copies under `docs/paper/paper_a/figures/data/` map
to exp07, exp10, exp11, exp13, exp19; each committed header was verified to
match its producer's summarize function exactly.

## Cross-cutting determinism notes

- No multiprocessing/threading/joblib anywhere in the producer chain; no
  wall-clock or `os.urandom` inputs; no hash-order-dependent output
  (`PYTHONHASHSEED`-irrelevant: iteration is insertion-ordered, explicit, or
  `sorted`).
- CSV float text is shortest-roundtrip `repr(float)` via `csv.DictWriter` —
  stable across platforms and Python ≥ 3.11.
- Fits and inversions are deterministic (vectorized grid search with
  first-minimum `argmin` tie-break; fixed-tolerance bisection). No BLAS-path
  operations, so cross-machine float variation reduces to libm
  transcendentals (exp16/17 note above) and numpy's per-version pairwise
  summation.
- Cosmetic only, no action: some seed-offset schemes can collide across
  roles/blocks (e.g. exp07's `10_000*2 + r = 20_000*1 + r`); this is
  deterministic and does not affect reproducibility of the recorded tables.

## Environment record (as of this audit)

`pyproject.toml`: `requires-python >= 3.11`, `numpy >= 2.0` (the
`np.trapezoid` API floor), `matplotlib >= 3.8`; dev extras `pytest >= 8`,
`ruff >= 0.4`; experiments extra `numba >= 0.59` (used only by
`positive_control/accelerated_two_orders.py`, not by any producer). No lock
file; CI floats `python 3.11` with unpinned `pip install -e ".[dev,experiments]"`.
The regeneration PR must record the concrete versions it runs under.
