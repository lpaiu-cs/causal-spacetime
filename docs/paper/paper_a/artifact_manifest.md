# Paper A evidence-artifact manifest

Status: **draft evidence package audit, 2026-07-17**.

The manuscript cites 19 generated summary tables. A clean checkout currently
contains five figure-source copies; the other 14 expected tables live under the
gitignored `outputs/data/` path only after rerunning their producers. This
manifest records the boundary instead of treating regenerability as equivalent
to archival availability. Before submission, regenerate all rows from the
declared baseline, validate them, and commit a provenance-locked 19-table
package.

| Producer | Expected summary | Clean-checkout evidence |
| --- | --- | --- |
| `exp03_causalset_timelike_reconstruction.py` | `outputs/data/timelike_reconstruction_summary.csv` | regenerate |
| `exp05_finite_speed_lattice_counterexample.py` | `outputs/data/finite_speed_lattice_growth.csv` | regenerate |
| `exp06_spacelike_distance_reconstruction.py` | `outputs/data/spacelike_distance_proxy_summary.csv` | regenerate |
| `exp07_timelike_pair_reconstruction_convergence.py` | `outputs/data/timelike_pair_reconstruction_summary.csv` | committed figure-source copy |
| `exp08_probe_pair_statistical_calibration.py` | `outputs/data/probe_pair_statistical_calibration_summary.csv` | regenerate |
| `exp09_longest_chain_calibration.py` | `outputs/data/longest_chain_calibration_summary.csv` | regenerate |
| `exp10_dimension_reconstruction.py` | `outputs/data/dimension_reconstruction_summary.csv` | committed figure-source copy |
| `exp11_discrete_observer_radar_reconstruction.py` | `outputs/data/discrete_radar_reconstruction_summary.csv` | committed figure-source copy |
| `exp12_single_observer_reflection_degeneracy.py` | `outputs/data/single_observer_reflection_degeneracy.csv` | regenerate |
| `exp13_oriented_radar_lorentz_map_recovery.py` | `outputs/data/oriented_radar_lorentz_summary.csv` | committed figure-source copy |
| `exp14_observer_atlas_consistency.py` | `outputs/data/observer_atlas_transition_summary.csv` | regenerate |
| `exp14_observer_atlas_consistency.py` | `outputs/data/observer_atlas_loop_summary.csv` | regenerate |
| `exp15_exact_poincare_map_sanity.py` | `outputs/data/exact_poincare_map_sanity.csv` | regenerate |
| `exp16_rindler_horizon_reconstruction.py` | `outputs/data/rindler_horizon_reconstruction_summary.csv` | regenerate |
| `exp17_inertial_vs_rindler_accessibility.py` | `outputs/data/inertial_vs_rindler_accessibility.csv` | regenerate |
| `exp18_conformal_order_ambiguity.py` | `outputs/data/conformal_order_ambiguity_summary.csv` | regenerate |
| `exp19_weighted_conformal_volume_reconstruction.py` | `outputs/data/weighted_conformal_volume_summary.csv` | committed figure-source copy |
| `exp20_conformal_volume_exact_sanity.py` | `outputs/data/conformal_volume_exact_sanity.csv` | regenerate |
| `exp23_thinning_coarse_graining_stability.py` | `outputs/data/thinning_coarse_graining_summary.csv` | regenerate |

## Committed figure-source hashes

These files are under `docs/paper/paper_a/figures/data/`.

| File | SHA-256 |
| --- | --- |
| `dimension_reconstruction_summary.csv` | `b3f72954e2aee94b46b1e83f2f9a4c543c9bf9b11883c56205683eaa4c268ee3` |
| `discrete_radar_reconstruction_summary.csv` | `8c797a174c35ed56bfdeca8ae4a57c2fb33364932e46326e8d80445eef606f06` |
| `oriented_radar_lorentz_summary.csv` | `fcab9242d8d7f0c4a53b9e85c584622befb62304393014cddf9cb3fb3c484ae4` |
| `timelike_pair_reconstruction_summary.csv` | `2210339e9420d12bb6b04bdd62e2763eec0ef1176471e66ce8932e958c631c8d` |
| `weighted_conformal_volume_summary.csv` | `6923274c8febdb4736ddcebff583468bdd39ca2f4b7860f07fa2ea69232e0d88` |

The foundation-layer baseline named by the manuscript is `325df55`. A
submission artifact should additionally record the actual regeneration commit,
interpreter and dependency versions, arguments/default configuration, seeds,
row counts, and hashes for all 19 tables.

## Section 6 capstone evidence bundle (P14) — no regeneration debt

Status: **complete, provenance-locked, 2026-08-09.** Every file below is
committed; the repository's tests recompute the capstone's metrics, verdicts,
and sentences from the stored raw samples, reproduce a prefix of every frozen
seed stream, and assert the commit-ancestry contract (`P' = 51875a2 ≺
F = b858b08 ≺ R = 5126ddf`; manuscript Appendix B). The submission gate for
the paper is therefore: the 19 legacy tables above, regenerated and committed,
PLUS this bundle as-is.

| File | SHA-256 |
| --- | --- |
| `docs/prereg/p14_prereg.md` | `9acca966b5f6d10c8c4bea244940bb64eff9de9cec56b93cc84ad74221d9f6fc` |
| `docs/prereg/p14_prereg_preflight.json` | `173b7063c1a6d45204bc7f22155e8253a153e18c66d0c4867e756c828e4ae34e` |
| `docs/prereg/p14_prereg_freeze.json` | `45eb6859dc755fd7ea3b9ecfb9cd19329895a0a9dd45633c26bad72b306bedd8` |
| `docs/prereg/p14_prereg_results.json` | `fe7e700fce1527960210b7b54f16308c90b30c83abe45aaabd51028626a488f4` |
| `docs/prereg/p14_s1_cost.json` | `63b8a6c4b82a897a869c0aecabed96fd68c4e6cc2baa94d26ed5d888634f7e9f` |
| `docs/prereg/p14_s1_cost.md` | `e5980c7af1ef510981d55446146597512799240fe1d667f836bb69721963597c` |
| `docs/prereg/p14_probe_p3e_results.json` | `f6ae9b1ed3154ab5d518dfd452cf0d5d85d4c0fd96268054f48a66a02a11ff0a` |
| `docs/prereg/p14_probe_p3c_results.json` | `e0dfffac971c28f7b1db0fbc5bc3364096755ac3260e4593dfb56f6a72c8166f` |
| `docs/prereg/p14_weyl_curvature.md` | `483c6870869d792997888ded1fc0fcc7c869dba17db0b76d06aafa3bbde80187` |
| `docs/prereg/p14_checks/p14_brinkmann_vacuum_check.py` | `83081bd7ee610329c30d719961b09e63f4858426c31514d0fb6f394e7ef41036` |
| `docs/prereg/p14_checks/p14_interval_volume_constant_a.py` | `ffb5acd109a41ca16a2cd12cf7a43dc83302e26cc95d89b443fe8652f5f6ec1f` |
| `experiments/positive_control/p14_prereg.py` | `f2be65bbfa663aa211b6ac7af97e7a51ba2e70abde0b502f1df63058fefa28db` |
| `tests/test_p14_prereg.py` | `7678bf52343437fbd11964c70f70dd921691c9d208aeeb74a3b0f93d00e7dbc0` |

Hashes are as of results commit `R` plus this paper-integration change to
`p14_weyl_curvature.md` (two stale `frozen/p14/` pointers corrected to the
actual bundle location; correction recorded in that document's status notes).
If any bundle file changes, this table must be regenerated in the same commit.

### Citation-to-artifact inventory (manuscript Section 6 / Appendix B)

| Manuscript location | Cited content | Evidence file |
| --- | --- | --- |
| 6.1 | 1+1D conformal ceiling, measure/counting channel framing | frozen P12/P13 records (`docs/prereg/frozen/p12/`, `docs/prereg/frozen/p13v3/`); `outputs/data/conformal_order_ambiguity_summary.csv` (exp18, legacy table debt) |
| 6.2 | det g = -1, vacuum, Weyl != 0 (exact) | `docs/prereg/p14_checks/p14_brinkmann_vacuum_check.py` |
| 6.2 | diamond volume +0.4% (wT=1) / +7.3% (wT=2) | `docs/prereg/p14_checks/p14_interval_volume_constant_a.py`; `docs/prereg/p14_weyl_curvature.md` §4.4 |
| 6.3 | operating point and primary-statistic selection | `docs/prereg/p14_probe_p3e_results.json` |
| 6.3 | frozen claims, margins, sizes, seed plan, verdict grammar | `docs/prereg/p14_prereg.md` |
| 6.3 | joint/null certifications (4000/4000; CP-lower 0.958 / 0.906) | `docs/prereg/p14_prereg_preflight.json` |
| 6.4 | C1/C2 result table, sentences, ambiguity zeros, complete separation | `docs/prereg/p14_prereg_results.json` |
| 6.4 | replication target of C2 | `docs/prereg/p14_probe_p3c_results.json` |
| 6.6 | S1 price (0.77 ms/pair, ~360x, projections) | `docs/prereg/p14_s1_cost.json`, `docs/prereg/p14_s1_cost.md` |
| App. B | freeze manifest contents, code_version chain | `docs/prereg/p14_prereg_freeze.json`; ancestry/gate tests in `tests/test_p14_prereg.py` |
| §10 | runner modes and gates | `experiments/positive_control/p14_prereg.py` |
