# Paper A evidence-artifact manifest

Status: **complete 19-table evidence package, 2026-08-09** (supersedes the
2026-07-17 draft audit that recorded a 5-committed / 14-regenerate boundary).

Every summary table the manuscript cites is now committed under
`docs/paper/paper_a/figures/data/` with its SHA-256 and row count below, all
regenerated in ONE recorded environment. A contract test
(`tests/test_paper_a_manifest.py`) recomputes every hash in this manifest
against the committed files, so manifest and evidence cannot drift silently.

| Producer | Expected summary | Clean-checkout evidence |
| --- | --- | --- |
| `exp03_causalset_timelike_reconstruction.py` | `outputs/data/timelike_reconstruction_summary.csv` | committed |
| `exp05_finite_speed_lattice_counterexample.py` | `outputs/data/finite_speed_lattice_growth.csv` | committed |
| `exp06_spacelike_distance_reconstruction.py` | `outputs/data/spacelike_distance_proxy_summary.csv` | committed |
| `exp07_timelike_pair_reconstruction_convergence.py` | `outputs/data/timelike_pair_reconstruction_summary.csv` | committed |
| `exp08_probe_pair_statistical_calibration.py` | `outputs/data/probe_pair_statistical_calibration_summary.csv` | committed |
| `exp09_longest_chain_calibration.py` | `outputs/data/longest_chain_calibration_summary.csv` | committed |
| `exp10_dimension_reconstruction.py` | `outputs/data/dimension_reconstruction_summary.csv` | committed |
| `exp11_discrete_observer_radar_reconstruction.py` | `outputs/data/discrete_radar_reconstruction_summary.csv` | committed |
| `exp12_single_observer_reflection_degeneracy.py` | `outputs/data/single_observer_reflection_degeneracy.csv` | committed |
| `exp13_oriented_radar_lorentz_map_recovery.py` | `outputs/data/oriented_radar_lorentz_summary.csv` | committed |
| `exp14_observer_atlas_consistency.py` | `outputs/data/observer_atlas_transition_summary.csv` | committed |
| `exp14_observer_atlas_consistency.py` | `outputs/data/observer_atlas_loop_summary.csv` | committed |
| `exp15_exact_poincare_map_sanity.py` | `outputs/data/exact_poincare_map_sanity.csv` | committed |
| `exp16_rindler_horizon_reconstruction.py` | `outputs/data/rindler_horizon_reconstruction_summary.csv` | committed |
| `exp17_inertial_vs_rindler_accessibility.py` | `outputs/data/inertial_vs_rindler_accessibility.csv` | committed |
| `exp18_conformal_order_ambiguity.py` | `outputs/data/conformal_order_ambiguity_summary.csv` | committed |
| `exp19_weighted_conformal_volume_reconstruction.py` | `outputs/data/weighted_conformal_volume_summary.csv` | committed |
| `exp20_conformal_volume_exact_sanity.py` | `outputs/data/conformal_volume_exact_sanity.csv` | committed |
| `exp23_thinning_coarse_graining_stability.py` | `outputs/data/thinning_coarse_graining_summary.csv` | committed |

## Regeneration record (2026-08-09)

- **Code state:** every producer and its transitive dependencies as of merge
  `f0ac576` — byte-identical to the foundation baseline `325df55` for all
  producers per the read-only audit (`producer_audit.md`); the sole changed
  dependency (`causal.py` `DEFAULT_CAUSAL_ATOL` extraction) is numerically
  neutral.
- **Environment:** CPython 3.11.9, numpy 2.4.6, matplotlib 3.11.1
  (Windows, single process). All 18 producers run from the repository root
  with DEFAULT arguments; all seeds are the scripts' committed defaults.
- **Validation against the prior copies:** of the 5 previously committed
  figure-source copies, `discrete_radar_reconstruction_summary.csv` and
  `oriented_radar_lorentz_summary.csv` reproduced byte-for-byte; the other
  three differed only at relative magnitude <= 1.5e-13 (last-ulp reduction
  noise across numpy versions — same seeds, same sample streams). Every
  manuscript-cited value was re-checked against the regenerated tables and
  is unchanged at its quoted precision, so no manuscript number required
  correction. The prior 5 copies are superseded by the regenerated versions
  so the whole package carries one recorded environment.

## Committed evidence hashes (all 19 tables)

Files are under `docs/paper/paper_a/figures/data/`. All hashes in this
manifest are SHA-256 of the LF-NORMALIZED content -- the committed blob
bytes -- because working-tree line endings vary by platform (git's
autocrlf); the contract test normalizes the same way.

| File | Rows | SHA-256 |
| --- | --- | --- |
| `timelike_reconstruction_summary.csv` | 4 | `81e043ed885be660c671f3e8be7bb3ca29e16509ec0b177b57201293186756f9` |
| `finite_speed_lattice_growth.csv` | 31 | `cab0a170292bef4cd5f01423429d5e7b5a207895589564fb24528947d8e99e8d` |
| `spacelike_distance_proxy_summary.csv` | 350 | `1560d61fee1559f7e848714265b9676bcd7a47ffdc2ad9247db17411f47bc063` |
| `timelike_pair_reconstruction_summary.csv` | 4 | `ebdda0f7fa8cb4115bad1084b7d7c75fbf4e865a9a8f1583408dd3fbf7171781` |
| `probe_pair_statistical_calibration_summary.csv` | 4 | `b64ae0d3cbabe8d57c1c7fbb85f3232d4020472012a3e3ed3accb3a95957c8b1` |
| `longest_chain_calibration_summary.csv` | 60 | `39342d38369acb49a387718db642d9b99105b5f1839baaabb609b6d7fb880a13` |
| `dimension_reconstruction_summary.csv` | 12 | `6050451ee926fd34d2abf0ae65b4de13723c41d45f51e681813601ccd2efd5cb` |
| `discrete_radar_reconstruction_summary.csv` | 16 | `faf22817729980c222ba8059ff9c779f4abbf950c4ab5ee8286eda0d013d9ef6` |
| `single_observer_reflection_degeneracy.csv` | 8 | `6411d92156e8b18a020b4194ec775f4be4af15db9fae34c68ff8e45c54069b60` |
| `oriented_radar_lorentz_summary.csv` | 18 | `275768cef792949e3d834a71003b24f2c2766a9eb084aed7172ba6016919a768` |
| `observer_atlas_transition_summary.csv` | 135 | `6d11f44e155657539989b6b529410c8e00e19de47b13272e26a0c9beb5bb53d8` |
| `observer_atlas_loop_summary.csv` | 45 | `8ef32bb1fd5f4eb2d4921d4ab2ddcd5121c04c5de07badbc440a3e076e3d12bb` |
| `exact_poincare_map_sanity.csv` | 4 | `757a4cee8e15fe8909ac0af5ac3a7fc3c7e81ef1f9df0d7d145c8cf0f5e57b34` |
| `rindler_horizon_reconstruction_summary.csv` | 18 | `42185976019e0309bd1eb7c14f2e4daffa97f1444453efe55ed17491d6122085` |
| `inertial_vs_rindler_accessibility.csv` | 800 | `941dd5ef5745bcc044c7c28fec65509341a9f745e256b5037f4373e48d0d9f4a` |
| `conformal_order_ambiguity_summary.csv` | 4 | `25db0f57b9c90239f7b047b8bc6c76e8e093217826e9d037bd8575182ac53aca` |
| `weighted_conformal_volume_summary.csv` | 9 | `36ba52bc7b7c0c1dd0db5d768592c74ad99d10151c3c0a81b5d33806eb08ade1` |
| `conformal_volume_exact_sanity.csv` | 3 | `53efec7033999761cb25e9db725497986c7af7138d72e8883fb5a3708064098e` |
| `thinning_coarse_graining_summary.csv` | 4 | `b1c09b841dcfe54bee38264d1948ea687adf01a8324875b1f8e955db05c79716` |

The foundation-layer baseline named by the manuscript is `325df55`.

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
