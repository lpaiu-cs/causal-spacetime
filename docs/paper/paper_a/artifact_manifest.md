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
