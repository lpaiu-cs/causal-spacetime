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
manifest are SHA-256 of the raw committed bytes. Line endings of every
hashed artifact are pinned to LF at the storage boundary
(`.gitattributes`, `eol=lf`), so the same digest is reported by plain
`sha256sum` on any checkout of any platform, and any byte change --
including a line-ending change -- breaks the contract test.

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

### Section 6.7 extension bundle (S4 Schwarzschild + supporting evidence)

Same digest convention (raw SHA-256 of the committed LF blob bytes). Roles:
`p14_s4_executed_freeze_manifest.json` is the immutable snapshot of the
manifest that governed the executed campaign (byte-identical to the blob at
freeze commit `ceed85d`, verified against the historical blobs by
`tests/test_s4_prereg.py`); `p14_s4_freeze_manifest.json` is the CURRENT
manifest, which hashes the post-result replay surface and is what the runner
enforces today. The C2 feasibility audit is exploratory (stored-data only).

| File | SHA-256 |
| --- | --- |
| `docs/prereg/p14_s4_schwarzschild_c1.md` | `0820043dcde99c64cb73d49d519d12fb25723588b2a62413e0dd2f1e008371f4` |
| `docs/prereg/p14_s4_executed_freeze_manifest.json` | `9e604fbc3d4e1399d7cfdc410a2f97e4bc7645ff2d1ceab4638b4c3cd9661071` |
| `docs/prereg/p14_s4_freeze_manifest.json` | `6b1e4d6ab64a08dea389785afe605befb0c61383b6e9aba9b1d7c773c4a9b804` |
| `docs/prereg/p14_s4_results.json` | `3c728bbeaab6ee1b02a700c972d7d77d542da0726098cb46a3ec08ec0fc5fd35` |
| `docs/prereg/p14_s3_probe_results.json` | `1bd37acc1a01780c6078c895fd0fa77adb7a8acbd324b2c12a7d41d7afd7ebc8` |
| `docs/prereg/p14_c2_feasibility_audit.json` | `ab133b2030a82d17f6a70629b1fedc58be52de1e5fe03bd8f5ad127548efa886` |
| `experiments/positive_control/s4_schwarzschild_c1.py` | `71a13cbc052d1a71e59119d4dcf385b3d3fb830a2f142dd9253690c481823805` |
| `experiments/positive_control/probe_seed_ledger.py` | `16720b9cf6c8c53890539214d3442f05b6587d830444603e044cab8e0954dafc` |
| `experiments/positive_control/c2_feasibility_audit.py` | `21810f35d94ea854611195875c0ac61b4a144c862184e4baa50edb4b00781715` |
| `tests/test_s4_prereg.py` | `06642c6b0f491fd63fe63d8c4a76a83f274ab4a820e26f788c314cdf2af8fdaa` |
| `docs/theory/schwarzschild_volume_oracle_note.md` | `e7c10656cd1fb86b9978096374950a1c438cd7e49bbff647a1350027c62c5ad1` |

### Section 6.7 second-stage bundle (S5 C2-unpaired discrimination)

Same digest convention. Roles mirror S4:
`p14_s5_executed_freeze_manifest.json` is the immutable snapshot of the
manifest that governed the executed campaign (byte-identical to the blob at
freeze commit `86e3674`, verified against the historical blobs by
`tests/test_s5_prereg.py`); `p14_s5_freeze_manifest.json` is the CURRENT
manifest, hashing the post-result replay surface enforced by the runner
today. The seed ledger row above covers both stages.

| File | SHA-256 |
| --- | --- |
| `docs/prereg/p14_s5_schwarzschild_c2.md` | `c9f8997336b543400b10aee9548be66a175ac5eee16d1ee906398b132273c340` |
| `docs/prereg/p14_s5_executed_freeze_manifest.json` | `0b4d28976076ef0e59122326cb53dee976f87c23aeaff5bb517b672f3d1ded08` |
| `docs/prereg/p14_s5_freeze_manifest.json` | `36098a39993539d22bcfa94a27b111882c66be9c2c91d24ec1b3108d6e564a58` |
| `docs/prereg/p14_s5_results.json` | `d5f339fb4b33d3cd048c93587f630b45314350e7919ec848d560fe7e8b8a9ad6` |
| `experiments/positive_control/s5_schwarzschild_c2.py` | `4761e5f90e019b8823fd93c883ece1a2f55034799a45f65d40d2cbcb0f3adc19` |
| `tests/test_s5_prereg.py` | `29525ef24558c6142eaccf44eb58e82f65a689c65d88bcc5c691505434798a4a` |

### Auxiliary O4b instrument-audit bundle (oracle arc)

Same digest convention. This bundle carries the auxiliary audit of
manuscript Sections 9-10 and Appendix B — an instrument statement that
upgrades nothing in Sections 6-6.7. Roles:
`p14_o3_volume.json` is the certified volume the audit consumed (write-once,
O3 execution from freeze `785148e`); `p14_o4b_executed_freeze_manifest.json`
is the immutable snapshot of the 25-file surface that governed the executed
O4b campaign (byte-identical to the freeze blob at `715865a`, digest equal to
the `manifest_digest` recorded inside the incident, checkpoint, and result;
verified by `tests/test_o4b_incident.py`); `p14_o4b_freeze_manifest.json` is
the CURRENT manifest, re-pinned by protocol maintenance (PR #78);
`p14_o4b_incident.json` and `p14_o4b_checkpoint.json` are the campaign run's
preserved records, verbatim as emitted (their SHA-256 is additionally pinned
inside `o4b_recover.py`, which refuses recovery on any byte change);
`p14_o4b_results.json` is the recovered verdict
(`run_kind = recovered_completed_campaign`: no new seed, no solver call, no
resampling, no gate change — stage verdict CONCORDANT per the frozen table).
No campaign rerun stands behind any row of this table.

| File | SHA-256 |
| --- | --- |
| `docs/prereg/p14_o3_volume.json` | `aed2ce111901a53fc59b343af77ec3499afcfa48a54fcb68bf76134d25a44a5f` |
| `docs/prereg/p14_o3_executed_freeze_manifest.json` | `350cc1786a858f099bbe7e094a8682ae25915287fd091503297500f60fd00907` |
| `docs/prereg/p14_o4b_freeze_manifest.json` | `793a8358ff7a037de0ba016882fe1af359357c3056608ef89981d222642a2981` |
| `docs/prereg/p14_o4b_executed_freeze_manifest.json` | `cec650b9391af0fc11e4b6bb94455cdbc1a037c18ecec83149d0d4693e7d7be2` |
| `docs/prereg/p14_o4b_incident.json` | `8106b16f0d03efe1acc81941f7ca3149cb8ddee0fbaaf00cfccb5c8376d4cc75` |
| `docs/prereg/p14_o4b_checkpoint.json` | `5dd6de7eea25c8384329b26ebf9f61c9dae51ee451c6d5c5717dbd463dc6679a` |
| `docs/prereg/p14_o4b_results.json` | `fe47e43b9f5c6618bbf30dd201a5a8b8d9ae69a712dbd12ea959c1b122d681d8` |
| `experiments/oracle/o4b_recover.py` | `d11c288771cd6401e558563b3bac3a683e8cff9ae033dabb4a0dab74b1ff70cb` |
| `tests/test_o4b_recover.py` | `9312270e7dd2665c1f00e177f32c680e1a5ad879929f78d8f8cce16c344f4263` |
| `tests/test_o4b_wiring_fixes.py` | `cbb88a9f0ba76ac37afb730f4cf279f3a6339948326e509c7161d5d63b37aa5f` |

### Citation-to-artifact inventory (manuscript Section 6 / Appendix B)

| Manuscript location | Cited content | Evidence file |
| --- | --- | --- |
| 6.1 | 1+1D conformal ceiling, measure/counting channel framing | frozen P12/P13 records (`docs/prereg/frozen/p12/`, `docs/prereg/frozen/p13v3/`); `figures/data/conformal_order_ambiguity_summary.csv` (exp18, committed) |
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
| 6.7 | frozen rule, gates, margins, power certification | `docs/prereg/p14_s4_schwarzschild_c1.md` |
| 6.7 | result table, verdict, frozen sentences, ambiguity/escalation | `docs/prereg/p14_s4_results.json` |
| 6.7 | exploration block (comparison arm of the replication gate) | `docs/prereg/p14_s3_probe_results.json` |
| 6.7 | C2 feasibility audit (exploratory, stored-data) | `docs/prereg/p14_c2_feasibility_audit.json`; `experiments/positive_control/c2_feasibility_audit.py` |
| 9 | oracle analytic reductions, certification closure, audit relation | `docs/theory/schwarzschild_volume_oracle_note.md`; `docs/theory/schwarzschild_volume_oracle_certification.md` |
| App. B | executed-freeze snapshot vs current replay manifest; seed ledger | `docs/prereg/p14_s4_executed_freeze_manifest.json`; `docs/prereg/p14_s4_freeze_manifest.json`; `experiments/positive_control/probe_seed_ledger.py`; contract tests in `tests/test_s4_prereg.py` |
| §10 | S4 runner gates and replay ownership | `experiments/positive_control/s4_schwarzschild_c1.py` |
| 6.7 | S5 frozen rule, margins, four-way gate, certification | `docs/prereg/p14_s5_schwarzschild_c2.md` |
| 6.7 | S5 result table, DETECTED outcome, BA verdict, frozen sentences | `docs/prereg/p14_s5_results.json` |
| App. B | S5 executed-freeze snapshot vs current replay manifest | `docs/prereg/p14_s5_executed_freeze_manifest.json`; `docs/prereg/p14_s5_freeze_manifest.json`; contract tests in `tests/test_s5_prereg.py` |
| §10 | S5 runner gates, fail-closed CLI, replay ownership | `experiments/positive_control/s5_schwarzschild_c2.py` |
| 9 | certified volume V = [56.212737, 57.348019], target-met | `docs/prereg/p14_o3_volume.json` |
| 9 | O4b audit verdict CONCORDANT, G1/G2 numbers, recovery flags | `docs/prereg/p14_o4b_results.json` |
| §10 | O4b recovery gates: input authentication, bit-exact reproduction, tamper refusals | `experiments/oracle/o4b_recover.py`; contract tests in `tests/test_o4b_recover.py` |
| §10 | O4b real-`main()` end-to-end regression (the wiring the abort exposed) | `tests/test_o4b_wiring_fixes.py` |
| App. B | O4b executed-freeze snapshot vs current manifest; preserved run records | `docs/prereg/p14_o4b_executed_freeze_manifest.json`; `docs/prereg/p14_o4b_freeze_manifest.json`; `docs/prereg/p14_o4b_incident.json`; `docs/prereg/p14_o4b_checkpoint.json`; contract tests in `tests/test_o4b_incident.py` |
| App. B | O4b incident narrative and recovery link (historical record, not regraded) | `docs/prereg/p14_o4b_incident.md` |
