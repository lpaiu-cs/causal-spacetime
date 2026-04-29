# causal-spacetime-lab

`causal-spacetime-lab` is a small Python research simulation project for
testing operational reconstructions of spacetime quantities from causal order
and related information-accessibility structure.

The project is scientifically conservative: these simulations are sanity checks
for reconstruction procedures and known relativistic behavior. They do not prove
a new theory of spacetime.

Milestones 18-26 refactor the theory-facing framing around locally finite
state-change causal trigger order and add the first finite diagnostics at that
primitive layer:

```text
state-changing events + causal trigger order + local finiteness
  -> observer-chain time order
  -> observer-slice-relative distance order
  -> calibrated effective metric representation
```

Metric geometry, seconds, meters, ratios, metric tensors, and curvature values
are treated as representation-layer objects in this program, not primitive
inputs. This does not mean metric geometry is unnecessary; it means a metric is
treated as an effective representation when order structures are sufficiently
consistent, calibrated, and representable.

## Current Project State

The project has grown into a lightweight validation suite for 1+1D
special-relativistic, causal-set, observer-protocol, and measure-dependent
reconstruction experiments:

- decompose events into radar coordinates for stationary and inertial observers,
- test Lorentz length contraction as lab-simultaneous endpoint event selection,
- sample events uniformly inside a causal diamond,
- build the causal order matrix,
- compute longest causal chains,
- count Alexandrov interval elements,
- estimate timelike proper time for internal event pairs when event density is
  supplied as additional scale information,
- compare independent probe-pair reconstruction errors against Poisson and
  fixed-N binomial finite-sampling expectations,
- estimate spacetime dimension from causal-order statistics in controlled flat
  Alexandrov intervals,
- reconstruct observer-chain radar coordinates from causal order plus supplied
  observer clock labels,
- reconstruct signed 1+1D radar coordinates when a second synchronized beacon
  chain supplies orientation reference structure,
- test Lorentz-map recovery between two oriented inertial observer protocols,
- test observer-atlas consistency across overlapping reconstructed charts with
  affine Lorentz/Poincare transition maps,
- test Rindler horizon reconstruction-inaccessibility as a controlled
  flat-spacetime horizon analogue,
- demonstrate conformal ambiguity and measure-dependent reconstruction in
  controlled 1+1D conformal toy models,
- test physical-volume sprinkling, local measure-shape recovery, and
  coarse-graining stability under random thinning with density rescaling,
- test order-first diagnostics for radar-return distance order, monotone
  invariance, calibration-driven ratio stability, and finite metric
  representability conditions,
- test ordinal embedding as a finite diagnostic for when observer-relative
  distance order admits a low-dimensional effective metric representation,
- test held-out order validation, null-model baselines, and subset stability
  for effective metric representations,
- test simultaneity-sliced observer-relative distance order using radar-time
  bins derived from causal order and observer tick order,
- test protocol-dependent cross-slice identification, transport, anchors, and
  persistence as representation-layer structure,
- test transport-gauge relational spatial evolution from persistence plus
  slice-local pair-distance order histories,
- test persistence ambiguity, identity matching, partial-label constraints, and
  hypothesis-dependent relational histories,
- consolidate the theory-facing layer around locally finite state-changing
  events and causal trigger order,
- test a minimal finite state-change causal trigger network as a strict
  partial order diagnostic,
- rank reference-chain candidates in finite state-change trigger networks
  using coverage, two-sided bracketing, interval profiles, and
  protocol-reference choice diagnostics,
- compute order-level predecessor/successor brackets, radar-time ranks, and
  bracket-width ranks from selected reference chains,
- compute same-emission echo-return positions and echo-delay ranks from
  selected reference chains,
- plant controlled echo-response motifs and test recovery of planted
  echo-delay ranks,
- classify shortcut returns and causal interference using return spectra,
  targeted shortcut stress tests, generic background edge perturbations,
  motif robustness checks, and reference-protocol dependence,
- test return-spectrum stability under closure-preserving event thinning,
  immediate-edge thinning, and reference-chain subsampling,
- test echo-response order signatures over motif populations, stable
  response-order cores, and protocol-variant response-rank order,
- test finite-speed lattice counterexamples and exploratory spacelike-distance
  proxies.

Natural units are used throughout, with `c = 1`.

## Installation

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

## Run Tests

```bash
pytest
```

## Main Validation Experiment

The main non-tautological timelike validation experiment from Milestone 2
samples many internal timelike pairs inside a larger causal diamond. It uses
global event density as additional scale information, then reconstructs each
pair's proper time from its own Alexandrov interval cardinality.

```bash
python experiments/exp07_timelike_pair_reconstruction_convergence.py
```

It writes:

- `outputs/data/timelike_pair_reconstruction_pairs.csv`
- `outputs/data/timelike_pair_reconstruction_summary.csv`
- `outputs/figures/timelike_pair_reconstruction_scatter.png`
- `outputs/figures/timelike_pair_reconstruction_error_vs_N.png`

Milestone 3 adds an independent probe-pair statistical calibration experiment.
Probe endpoints are sampled independently from the support sprinkle and are not
inserted into the support event set.

```bash
python experiments/exp08_probe_pair_statistical_calibration.py
```

It writes:

- `outputs/data/probe_pair_statistical_calibration_pairs.csv`
- `outputs/data/probe_pair_statistical_calibration_summary.csv`
- `outputs/data/probe_pair_statistical_calibration_binned_by_tau.csv`
- `outputs/figures/probe_pair_tau_scatter.png`
- `outputs/figures/probe_pair_error_vs_tau.png`
- `outputs/figures/probe_pair_rmse_vs_N.png`
- `outputs/figures/probe_pair_relative_error_by_tau_bin.png`

An optional longest-chain calibration can be run with:

```bash
python experiments/exp09_longest_chain_calibration.py
```

It writes:

- `outputs/data/longest_chain_calibration_summary.csv`
- `outputs/figures/longest_chain_calibration.png`

Milestone 4 adds Myrheim-Meyer dimension reconstruction:

```bash
python experiments/exp10_dimension_reconstruction.py
```

It writes:

- `outputs/data/dimension_reconstruction_results.csv`
- `outputs/data/dimension_reconstruction_summary.csv`
- `outputs/figures/dimension_estimate_vs_N.png`
- `outputs/figures/relation_fraction_vs_dimension.png`
- `outputs/figures/dimension_error_vs_N.png`

Milestone 5 adds causal-order-based observer-chain radar reconstruction:

```bash
python experiments/exp11_discrete_observer_radar_reconstruction.py
```

It writes:

- `outputs/data/discrete_radar_reconstruction_events.csv`
- `outputs/data/discrete_radar_reconstruction_summary.csv`
- `outputs/figures/discrete_radar_time_scatter.png`
- `outputs/figures/discrete_radar_distance_scatter.png`
- `outputs/figures/discrete_radar_error_vs_ticks.png`
- `outputs/figures/discrete_radar_accessible_fraction.png`

Milestone 6 adds oriented two-chain radar reconstruction and Lorentz-map
recovery:

```bash
python experiments/exp12_single_observer_reflection_degeneracy.py
python experiments/exp13_oriented_radar_lorentz_map_recovery.py
```

It writes:

- `outputs/data/single_observer_reflection_degeneracy.csv`
- `outputs/figures/single_observer_reflection_degeneracy.png`
- `outputs/data/oriented_radar_lorentz_events.csv`
- `outputs/data/oriented_radar_lorentz_summary.csv`
- `outputs/figures/oriented_radar_lab_position_scatter.png`
- `outputs/figures/oriented_radar_moving_position_scatter.png`
- `outputs/figures/oriented_radar_lorentz_residual_vs_ticks.png`
- `outputs/figures/oriented_radar_beta_fit_vs_ticks.png`
- `outputs/figures/oriented_radar_accessible_fraction.png`

Milestone 7 adds observer-atlas consistency:

```bash
python experiments/exp14_observer_atlas_consistency.py
python experiments/exp15_exact_poincare_map_sanity.py
```

It writes:

- `outputs/data/observer_atlas_chart_events.csv`
- `outputs/data/observer_atlas_transition_summary.csv`
- `outputs/data/observer_atlas_loop_summary.csv`
- `outputs/data/exact_poincare_map_sanity.csv`
- `outputs/figures/observer_atlas_transition_beta_error_vs_ticks.png`
- `outputs/figures/observer_atlas_transition_rmse_vs_ticks.png`
- `outputs/figures/observer_atlas_invariant_disagreement_vs_ticks.png`
- `outputs/figures/observer_atlas_overlap_fraction_vs_ticks.png`
- `outputs/figures/observer_atlas_loop_consistency_vs_ticks.png`

Milestone 8 adds Rindler horizon reconstruction-inaccessibility:

```bash
python experiments/exp16_rindler_horizon_reconstruction.py
python experiments/exp17_inertial_vs_rindler_accessibility.py
```

It writes:

- `outputs/data/rindler_horizon_reconstruction_events.csv`
- `outputs/data/rindler_horizon_reconstruction_summary.csv`
- `outputs/data/inertial_vs_rindler_accessibility.csv`
- `outputs/figures/rindler_accessibility_map.png`
- `outputs/figures/rindler_accessible_fraction_vs_ticks.png`
- `outputs/figures/rindler_radar_time_scatter.png`
- `outputs/figures/rindler_radar_distance_scatter.png`
- `outputs/figures/rindler_error_vs_ticks.png`
- `outputs/figures/rindler_false_positive_negative_vs_ticks.png`
- `outputs/figures/inertial_vs_rindler_accessibility.png`

Milestone 9 adds conformal ambiguity and measure-dependent reconstruction:

```bash
python experiments/exp18_conformal_order_ambiguity.py
python experiments/exp19_weighted_conformal_volume_reconstruction.py
python experiments/exp20_conformal_volume_exact_sanity.py
```

It writes:

- `outputs/data/conformal_order_ambiguity_summary.csv`
- `outputs/data/weighted_conformal_volume_pairs.csv`
- `outputs/data/weighted_conformal_volume_summary.csv`
- `outputs/data/conformal_volume_exact_sanity.csv`
- `outputs/figures/conformal_order_ambiguity_scales.png`
- `outputs/figures/weighted_conformal_volume_scatter.png`
- `outputs/figures/weighted_conformal_volume_rmse_vs_N.png`
- `outputs/figures/weighted_conformal_volume_bias_by_profile.png`

Milestone 10 adds measure encoding and coarse-graining stability:

```bash
python experiments/exp21_physical_measure_sprinkling.py
python experiments/exp22_local_measure_profile_estimation.py
python experiments/exp23_thinning_coarse_graining_stability.py
python experiments/exp24_measure_sprinkling_exact_sanity.py
```

It writes:

- `outputs/data/physical_measure_sprinkling_pairs.csv`
- `outputs/data/physical_measure_sprinkling_summary.csv`
- `outputs/data/local_measure_profile_bins.csv`
- `outputs/data/local_measure_profile_summary.csv`
- `outputs/data/thinning_coarse_graining_pairs.csv`
- `outputs/data/thinning_coarse_graining_summary.csv`
- `outputs/data/measure_sprinkling_exact_sanity.csv`
- `outputs/figures/physical_measure_volume_scatter.png`
- `outputs/figures/physical_measure_rmse_vs_N.png`
- `outputs/figures/physical_measure_bias_by_profile.png`
- `outputs/figures/local_measure_profile_shape.png`
- `outputs/figures/local_measure_profile_rmse_vs_N.png`
- `outputs/figures/thinning_volume_error_vs_keep_probability.png`
- `outputs/figures/thinning_dimension_vs_keep_probability.png`
- `outputs/figures/thinning_count_ratio.png`

Milestone 11 adds order-first diagnostics and metric-representation conditions:

```bash
python experiments/exp25_radar_return_distance_order.py
python experiments/exp26_metric_representation_scale_invariance.py
python experiments/exp27_ratio_stability_from_calibration.py
python experiments/exp28_oriented_chart_distance_order_preservation.py
python experiments/exp29_metric_representability_diagnostics.py
python experiments/exp30_ordinal_exact_sanity.py
```

It writes:

- `outputs/data/radar_return_distance_order.csv`
- `outputs/data/metric_representation_scale_invariance.csv`
- `outputs/data/ratio_stability_from_calibration.csv`
- `outputs/data/oriented_chart_distance_order_preservation.csv`
- `outputs/data/metric_representability_diagnostics.csv`
- `outputs/data/ordinal_exact_sanity.csv`
- `outputs/figures/radar_return_order_inversion_vs_ticks.png`
- `outputs/figures/radar_return_order_scatter.png`
- `outputs/figures/metric_representation_order_preservation.png`
- `outputs/figures/ratio_stability_from_calibration.png`
- `outputs/figures/oriented_chart_distance_order_inversion_vs_ticks.png`

Milestone 12 adds ordinal embedding and effective metric representation
diagnostics:

```bash
python experiments/exp31_ordinal_embedding_recovery.py
python experiments/exp32_ordinal_dimension_selection.py
python experiments/exp33_noisy_incomplete_order_embedding.py
python experiments/exp34_observer_distance_order_embedding.py
python experiments/exp35_ordinal_embedding_exact_sanity.py
```

It writes:

- `outputs/data/ordinal_embedding_recovery.csv`
- `outputs/data/ordinal_dimension_selection.csv`
- `outputs/data/noisy_incomplete_order_embedding.csv`
- `outputs/data/observer_distance_order_embedding.csv`
- `outputs/data/ordinal_embedding_exact_sanity.csv`
- `outputs/figures/ordinal_embedding_violation_vs_constraints.png`
- `outputs/figures/ordinal_embedding_rmse_vs_constraints.png`
- `outputs/figures/ordinal_dimension_stress_curve.png`
- `outputs/figures/ordinal_dimension_violation_curve.png`
- `outputs/figures/noisy_order_embedding_violation.png`
- `outputs/figures/noisy_order_embedding_rmse.png`
- `outputs/figures/observer_distance_order_embedding_violation_vs_ticks.png`
- `outputs/figures/observer_distance_order_embedding_rmse_vs_ticks.png`

Milestone 13 adds held-out validation, null-model baselines, and representation
stability diagnostics:

```bash
python experiments/exp36_heldout_ordinal_embedding_validation.py
python experiments/exp37_embedding_stability_under_subsampling.py
python experiments/exp38_observer_order_null_baseline.py
python experiments/exp39_effective_metric_complexity_curve.py
python experiments/exp40_representation_stability_exact_sanity.py
```

It writes:

- `outputs/data/heldout_ordinal_embedding_validation.csv`
- `outputs/data/embedding_stability_under_subsampling.csv`
- `outputs/data/observer_order_null_baseline.csv`
- `outputs/data/effective_metric_complexity_curve.csv`
- `outputs/data/representation_stability_exact_sanity.csv`
- `outputs/figures/heldout_violation_by_model.png`
- `outputs/figures/heldout_generalization_gap.png`
- `outputs/figures/embedding_procrustes_stability_vs_constraints.png`
- `outputs/figures/embedding_order_stability_vs_constraints.png`
- `outputs/figures/observer_order_vs_null_test_violation.png`
- `outputs/figures/observer_order_vs_null_alignment_rmse.png`
- `outputs/figures/observer_order_vs_null_distance_order_error.png`
- `outputs/figures/effective_metric_complexity_curve.png`
- `outputs/figures/effective_metric_penalized_score.png`

Milestone 14 adds simultaneity-sliced observer-relative distance-order
diagnostics:

```bash
python experiments/exp41_radar_time_order_from_brackets.py
python experiments/exp42_same_slice_distance_order_preservation.py
python experiments/exp43_sliced_observer_order_null_baseline.py
python experiments/exp44_slice_width_sensitivity.py
python experiments/exp45_spatial_slice_exact_sanity.py
```

It writes:

- `outputs/data/radar_time_order_from_brackets.csv`
- `outputs/data/same_slice_distance_order_preservation.csv`
- `outputs/data/sliced_observer_order_null_baseline.csv`
- `outputs/data/slice_width_sensitivity.csv`
- `outputs/data/spatial_slice_exact_sanity.csv`
- `outputs/figures/radar_time_order_inversion_vs_ticks.png`
- `outputs/figures/same_slice_distance_order_inversion_vs_ticks.png`
- `outputs/figures/same_slice_vs_all_pairs_inversion.png`
- `outputs/figures/same_slice_pair_count_vs_bin_width.png`
- `outputs/figures/sliced_observer_vs_null_test_violation.png`
- `outputs/figures/sliced_observer_vs_null_alignment_rmse.png`
- `outputs/figures/sliced_observer_vs_null_distance_order_error.png`
- `outputs/figures/slice_width_pair_count.png`
- `outputs/figures/slice_width_distance_order_error.png`

Milestone 15 adds protocol-dependent cross-slice identification, transport,
anchor, and persistence diagnostics:

```bash
python experiments/exp46_cross_slice_predicate_undefined.py
python experiments/exp47_sliced_constraint_graph_decomposition.py
python experiments/exp48_slice_local_embedding_validation.py
python experiments/exp49_slice_gauge_dependence.py
python experiments/exp50_anchor_constrained_transport.py
python experiments/exp51_persistence_dependent_velocity.py
python experiments/exp52_noisy_transport_sensitivity.py
python experiments/exp53_cross_slice_transport_exact_sanity.py
```

It writes:

- `outputs/data/cross_slice_predicate_undefined.csv`
- `outputs/data/sliced_constraint_graph_decomposition.csv`
- `outputs/data/slice_local_embedding_validation.csv`
- `outputs/data/slice_gauge_dependence.csv`
- `outputs/data/anchor_constrained_transport.csv`
- `outputs/data/persistence_dependent_velocity.csv`
- `outputs/data/noisy_transport_sensitivity.csv`
- `outputs/data/cross_slice_transport_exact_sanity.csv`
- `outputs/figures/sliced_constraint_components_vs_bin_width.png`
- `outputs/figures/sliced_constraint_largest_component_vs_bin_width.png`
- `outputs/figures/slice_local_order_error_vs_ticks.png`
- `outputs/figures/slice_local_rmse_vs_ticks.png`
- `outputs/figures/slice_gauge_same_slice_vs_global_error.png`
- `outputs/figures/slice_gauge_cross_slice_judgments.png`
- `outputs/figures/anchor_transport_global_rmse.png`
- `outputs/figures/anchor_transport_distance_order_error.png`
- `outputs/figures/persistence_velocity_by_transport.png`
- `outputs/figures/noisy_transport_global_rmse.png`
- `outputs/figures/noisy_transport_velocity_instability.png`

Milestone 16 adds transport-gauge relational spatial-evolution diagnostics:

```bash
python experiments/exp54_predicate_definability_table.py
python experiments/exp55_relational_shape_history_without_transport.py
python experiments/exp56_relational_history_gauge_invariance.py
python experiments/exp57_observer_slice_relational_evolution.py
python experiments/exp58_relational_invariants_vs_velocity.py
python experiments/exp59_relational_evolution_exact_sanity.py
```

It writes:

- `outputs/data/predicate_definability_table.csv`
- `outputs/data/relational_shape_history_without_transport.csv`
- `outputs/data/relational_history_gauge_invariance.csv`
- `outputs/data/observer_slice_relational_evolution.csv`
- `outputs/data/relational_invariants_vs_velocity.csv`
- `outputs/data/relational_evolution_exact_sanity.csv`
- `outputs/figures/relational_shape_change_rates.png`
- `outputs/figures/relational_history_gauge_invariance.png`
- `outputs/figures/observer_slice_relational_evolution_error.png`
- `outputs/figures/relational_change_vs_velocity_transport.png`

Milestone 17 adds persistence ambiguity, identity matching, and
relational-history hypothesis diagnostics:

```bash
python experiments/exp60_persistence_predicate_undefined.py
python experiments/exp61_symmetric_persistence_ambiguity.py
python experiments/exp62_relational_persistence_matching_recovery.py
python experiments/exp63_partial_label_constrained_persistence.py
python experiments/exp64_crossing_persistence_failure.py
python experiments/exp65_persistence_hypothesis_dependence.py
python experiments/exp66_persistence_matching_exact_sanity.py
```

It writes:

- `outputs/data/persistence_predicate_undefined.csv`
- `outputs/data/symmetric_persistence_ambiguity.csv`
- `outputs/data/relational_persistence_matching_recovery.csv`
- `outputs/data/partial_label_constrained_persistence.csv`
- `outputs/data/crossing_persistence_failure.csv`
- `outputs/data/persistence_hypothesis_dependence.csv`
- `outputs/data/persistence_matching_exact_sanity.csv`
- `outputs/figures/symmetric_persistence_ambiguity_gap.png`
- `outputs/figures/persistence_matching_accuracy_vs_motion.png`
- `outputs/figures/persistence_matching_ambiguity_gap_vs_motion.png`
- `outputs/figures/partial_label_matching_accuracy.png`
- `outputs/figures/partial_label_ambiguity_gap.png`
- `outputs/figures/crossing_persistence_track_error.png`
- `outputs/figures/persistence_hypothesis_change_rates.png`

Milestone 18 is a theory refactor rather than a new simulation milestone. It
adds:

- `docs/theory/state_change_causal_order.md`
- `docs/theory/order_vs_metric_representation.md`
- `docs/theory/relativity_compatibility_conditions.md`
- `docs/theory/quantum_extension_requirements.md`
- `docs/theory/milestone_reinterpretation_under_state_change_order.md`
- `docs/theory/paper_outline_state_change_order.md`

It also adds:

```bash
python scripts/check_theory_language.py
```

to guard against risky overclaim language in the theory docs.

Milestone 19 adds the first minimal state-change causal trigger toy model:

```bash
python experiments/exp67_state_change_exact_sanity.py
python experiments/exp68_state_change_toy_model.py
python experiments/exp69_state_change_observer_chain_diagnostic.py
```

It writes:

- `outputs/data/state_change_exact_sanity.csv`
- `outputs/data/state_change_toy_model_summary.csv`
- `outputs/data/state_change_observer_chain_diagnostic.csv`
- `outputs/figures/state_change_relation_density.png`
- `outputs/figures/state_change_max_interval_size.png`
- `outputs/figures/state_change_events_per_system.png`
- `outputs/figures/state_change_observer_chain_coverage.png`

This is a finite DAG diagnostic for the primitive state-change order layer. It
does not reconstruct metric geometry, extract observers automatically, add
finite-speed spatial geometry, derive quantum mechanics, or model curved
spacetime.

Milestone 20 adds reference-chain utility diagnostics:

```bash
python experiments/exp70_observer_chain_exact_sanity.py
python experiments/exp71_observer_chain_candidate_ranking.py
python experiments/exp72_observer_chain_coverage_vs_trigger_probability.py
python experiments/exp73_observer_chain_interval_profile.py
```

It writes:

- `outputs/data/observer_chain_exact_sanity.csv`
- `outputs/data/observer_chain_candidate_ranking.csv`
- `outputs/data/observer_chain_candidate_ambiguity.csv`
- `outputs/data/observer_chain_coverage_vs_trigger_probability.csv`
- `outputs/data/observer_chain_interval_profile.csv`
- `outputs/figures/observer_chain_score_by_source.png`
- `outputs/figures/observer_chain_bracketing_by_source.png`
- `outputs/figures/observer_chain_ambiguity_gap.png`
- `outputs/figures/observer_chain_coverage_vs_trigger_probability.png`
- `outputs/figures/observer_chain_bracketing_vs_trigger_probability.png`
- `outputs/figures/observer_chain_interval_cv_by_source.png`
- `outputs/figures/observer_chain_interval_profile_example.png`

These experiments compare local-system chains, order-only heuristic chains,
longest chains, and random baselines as candidate reference protocols. The
scores are finite reference-chain utility diagnostics, not clock calibration
and not an observer-reality claim.

Milestone 21 adds order-level bracket diagnostics from selected reference
chains:

```bash
python experiments/exp74_state_change_reference_bracket_diagnostics.py
python experiments/exp75_state_change_bracket_rank_reference_dependence.py
python experiments/exp76_state_change_reference_bracket_coverage_vs_trigger_density.py
python experiments/exp77_state_change_reference_bracket_exact_sanity.py
```

It writes:

- `outputs/data/state_change_reference_bracket_diagnostics.csv`
- `outputs/data/state_change_bracket_rank_reference_dependence.csv`
- `outputs/data/state_change_reference_bracket_coverage_vs_trigger_density.csv`
- `outputs/data/state_change_reference_bracket_exact_sanity.csv`
- bracket accessibility, rank-dependence, and trigger-density figures under
  `outputs/figures/`.

These rank-level diagnostics use only causal trigger order and selected
reference-chain positions. They do not reconstruct metric radar distance,
calibrate clocks, or define spatial geometry.

Milestone 22 adds same-emission echo-order diagnostics:

```bash
python experiments/exp78_state_change_echo_exact_sanity.py
python experiments/exp79_state_change_echo_order_diagnostics.py
python experiments/exp80_state_change_echo_reference_dependence.py
python experiments/exp81_state_change_echo_emission_sensitivity.py
python experiments/exp82_state_change_echo_coverage_vs_trigger_density.py
```

It writes:

- `outputs/data/state_change_echo_exact_sanity.csv`
- `outputs/data/state_change_echo_order_diagnostics.csv`
- `outputs/data/state_change_echo_reference_dependence.csv`
- `outputs/data/state_change_echo_emission_sensitivity.csv`
- `outputs/data/state_change_echo_coverage_vs_trigger_density.csv`
- echo reachability, reference-dependence, emission-sensitivity, and
  trigger-density figures under `outputs/figures/`.

These diagnostics fix a reference emission position and ask which target
events return later to the same selected reference chain. Echo-delay rank is a
rank/order-level diagnostic. It does not define physical distance, reconstruct
metric radar distance, implement finite-speed spatial geometry, or calibrate
time.

Milestone 23 adds controlled echo-response motif diagnostics:

```bash
python experiments/exp83_echo_motif_exact_sanity.py
python experiments/exp84_planted_echo_motif_recovery.py
python experiments/exp85_echo_motif_background_interference.py
python experiments/exp86_echo_motif_density_resolution.py
python experiments/exp87_echo_motif_reference_choice_visibility.py
```

It writes:

- `outputs/data/echo_motif_exact_sanity.csv`
- `outputs/data/planted_echo_motif_recovery.csv`
- `outputs/data/planted_echo_motif_recovery_summary.csv`
- `outputs/data/echo_motif_background_interference.csv`
- `outputs/data/echo_motif_density_resolution.csv`
- `outputs/data/echo_motif_reference_choice_visibility.csv`
- motif recovery, shortcut, tie-resolution, and reference-visibility figures
  under `outputs/figures/`.

These experiments insert controlled causal trigger structures with planted
echo-delay ranks. The planted rank is a validation label, not physical
distance or calibrated time. Shortcut returns are recorded as background
causal interference.

Milestone 24 adds echo shortcut and interference classification:

```bash
python experiments/exp88_echo_shortcut_exact_sanity.py
python experiments/exp89_echo_shortcut_injection_sweep.py
python experiments/exp90_echo_background_edge_perturbation.py
python experiments/exp91_echo_motif_path_length_robustness.py
python experiments/exp92_echo_shortcut_reference_dependence.py
python experiments/exp93_echo_interference_exact_sanity.py
```

It writes:

- `outputs/data/echo_shortcut_exact_sanity.csv`
- `outputs/data/echo_shortcut_injection_sweep.csv`
- `outputs/data/echo_background_edge_perturbation.csv`
- `outputs/data/echo_motif_path_length_robustness.csv`
- `outputs/data/echo_shortcut_reference_dependence.csv`
- `outputs/data/echo_interference_exact_sanity.csv`
- shortcut fraction, shortcut-depth, background-perturbation, path-robustness,
  and reference-dependence figures under `outputs/figures/`.

The key diagnostic is the return spectrum for a motif target relative to a
chosen reference chain and emission position. Targeted shortcut injection and
generic acyclic background edge perturbations are distinct stress tests; both
classify causal-order interference before any stronger interpretation is
attempted.

Milestone 25 adds return-spectrum stability diagnostics under coarse-graining:

```bash
python experiments/exp94_echo_coarse_graining_exact_sanity.py
python experiments/exp95_echo_event_thinning_stability.py
python experiments/exp96_echo_reference_subsampling_resolution.py
python experiments/exp97_echo_edge_thinning_fragility.py
python experiments/exp98_echo_shortcut_classification_under_coarse_graining.py
python experiments/exp99_echo_return_spectrum_stability_exact_sanity.py
```

It writes:

- `outputs/data/echo_coarse_graining_exact_sanity.csv`
- `outputs/data/echo_event_thinning_stability.csv`
- `outputs/data/echo_reference_subsampling_resolution.csv`
- `outputs/data/echo_edge_thinning_fragility.csv`
- `outputs/data/echo_shortcut_classification_under_coarse_graining.csv`
- `outputs/data/echo_return_spectrum_stability_exact_sanity.csv`
- event-thinning, reference-subsampling, edge-thinning, and classification
  stability figures under `outputs/figures/`.

Closure-preserving event thinning integrates out intermediate events while
preserving retained-event reachability. Immediate-edge thinning changes the
trigger graph before closure. Reference-chain subsampling reduces rank
resolution and can merge return ranks into ties.

Milestone 26 adds stable echo-response order signatures:

```bash
python experiments/exp100_echo_response_signature_exact_sanity.py
python experiments/exp101_layered_echo_response_order_recovery.py
python experiments/exp102_echo_response_signature_coarse_protocol_stability.py
python experiments/exp103_echo_response_shortcut_robust_core.py
python experiments/exp104_echo_response_reference_protocol_dependence.py
python experiments/exp105_echo_response_order_precondition_diagnostics.py
python experiments/exp106_echo_response_signature_stability_exact_sanity.py
```

It writes:

- `outputs/data/echo_response_signature_exact_sanity.csv`
- `outputs/data/layered_echo_response_order_recovery.csv`
- `outputs/data/echo_response_signature_coarse_protocol_stability.csv`
- `outputs/data/echo_response_shortcut_robust_core.csv`
- `outputs/data/echo_response_reference_protocol_dependence.csv`
- `outputs/data/echo_response_order_precondition_diagnostics.csv`
- `outputs/data/echo_response_signature_stability_exact_sanity.csv`
- response-signature recovery, stable-core, shortcut-robustness, reference
  dependence, and precondition figures under `outputs/figures/`.

Response-order signatures are ordinal diagnostics over motif populations. A
stable response-order core is candidate input for later representability tests,
not a spatial-distance claim or metric reconstruction.

## Other Experiments

The original full-diamond timelike reconstruction sanity check can be run with:

```bash
python experiments/exp03_causalset_timelike_reconstruction.py
```

The script writes:

- `outputs/data/timelike_reconstruction_summary.csv`
- `outputs/figures/timelike_reconstruction_error.png`

The Lorentz length-contraction experiment can be run with:

```bash
python experiments/exp02_lorentz_length_contraction.py
```

It writes:

- `outputs/data/lorentz_length_contraction_summary.csv`
- `outputs/figures/lorentz_length_contraction.png`

The finite-speed lattice counterexample can be run with:

```bash
python experiments/exp05_finite_speed_lattice_counterexample.py
```

It writes:

- `outputs/data/finite_speed_lattice_growth.csv`
- `outputs/figures/finite_speed_lattice_cones.png`
- `outputs/figures/finite_speed_lattice_count_growth.png`

The exploratory spacelike-distance proxy experiment can be run with:

```bash
python experiments/exp06_spacelike_distance_reconstruction.py
```

It writes:

- `outputs/data/spacelike_distance_proxy_summary.csv`
- `outputs/figures/spacelike_distance_proxy_scatter.png`

Run the lightweight suite with:

```bash
python experiments/run_all.py
```

## What The Result Means

The validation experiments check whether timelike separation in known 1+1D
Minkowski spacetime can be reconstructed using causal interval cardinality once
an event density is specified. Milestone 3 asks whether observed reconstruction
errors are consistent with finite-sampling expectations, or whether they suggest
bias, boundary effects, estimator mistakes, or incorrect conventions.

Milestone 4 asks whether dimension is recoverable as an order-statistical
observable in controlled flat Alexandrov intervals. This is part of a
mathematical reconstruction program:

```text
primitive causal/information-accessibility structure
  + counting measure / event density
  + observer protocol
  + orientation/reference protocol
  -> operational time, distance, dimension, coordinate transformations,
     atlas consistency diagnostics, horizon-limited reconstruction,
     measure-dependent metric-scale reconstruction, and coarse-graining
     stability checks
```

The longest chain is also reported as a causal-order observable with a simple
1+1D asymptotic normalization. Its normalization is finite-size sensitive.

This result does not show that spacetime is made of information, and it does not
derive metric scale from causal order alone. The interval-count estimate uses
event density as additional structure, which is exactly the point being tested.
The radar-coordinate functions implement a standard operational coordinate
protocol for specified observers; they are not a new theory of spacetime.
The length-contraction script illustrates the standard special-relativistic
event-selection issue: a lab-frame length uses endpoint events simultaneous in
the lab frame.
The finite-speed lattice script shows a counterexample to the weaker claim that
finite signal speed alone implies Lorentz-invariant spacetime structure.
The spacelike-distance proxy experiment is exploratory. Common-past,
common-future, and enclosing-interval counts are boundary-dependent diagnostics,
not validated estimators of spacelike distance.
Agreement with Poisson or binomial sampling-noise estimates is a consistency
check in a known spacetime model, not a proof of a new theory.
Dimension reconstruction is likewise controlled validation inside known causal
intervals; it does not show that dimension is purely information or that
spacetime has been derived.
Observer-chain radar reconstruction uses a supplied observer protocol and clock
labels. It tests operational spatial decomposition from causal accessibility,
but it does not show that causal order alone gives radar distance.
Milestone 6 makes the reflection degeneracy explicit: a single observer gives
unsigned distance only. Signed coordinates require a supplied orientation
reference such as a synchronized beacon chain with known separation. Lorentz-map
recovery is therefore a controlled validation between observer protocols, not a
derivation of Lorentz transformations from causal order alone.
Milestone 7 extends this to a small observer atlas. It fits affine
Lorentz/Poincare transition maps on chart overlaps and checks invariant
interval agreement and loop consistency. Observer origins, clocks, beacon
separations, and orientation references remain supplied protocol structure.
Milestone 8 uses a Rindler observer in flat 1+1D Minkowski spacetime to test
reconstruction-inaccessibility. It distinguishes ideal Rindler wedge
accessibility from finite-chain coverage. This is a controlled horizon analogue,
not a black hole simulation or a derivation of horizons from causal order alone.
Milestone 9 makes the conformal ambiguity explicit. Positive conformal
rescalings preserve causal order while changing physical volume and clock
scale. Weighted conformal volume reconstruction uses supplied measure weights;
the conformal factor is not derived from causal order alone.
Milestone 10 tests two ways that measure information can enter a reconstruction
program: as physical-volume event distribution and as a density scale that must
be rescaled under thinning. Local relative measure shape can be recovered
statistically from nonuniform counts in controlled conformal toy models, but
global constant conformal scale remains underdetermined without an absolute
density scale. Random thinning is a coarse-graining check; corrected density
should remain stable while uncorrected density should fail.
Milestone 11 reframes these reconstructions as representation-layer tests.
Primitive temporal structure is causal order; primitive spatial structure is
observer-relative distance order. Positive monotone transformations can preserve
order while changing ratios, so ratio stability requires calibration,
concatenation, repeated processes, or dynamics. Not every distance order admits
a useful low-dimensional metric representation.
Milestone 12 tests this representability question directly with ordinal
embedding. The experiments ask whether low-dimensional coordinates can act as a
low-complexity compression of distance-order constraints, how candidate
dimension affects ordinal stress, how noise and incomplete comparisons degrade
the representation, and whether observer-derived distance order supports an
effective 1D spatial embedding. These are finite diagnostics, not
representation theorems or proofs of spacetime emergence.
Milestone 13 adds held-out validation and null-model baselines. It asks whether
structured geometric or observer-derived order constraints generalize better
than shuffled or random constraints, whether independent subsets yield stable
embeddings, and whether low-dimensional complexity curves distinguish
structured order from null order. These checks support the interpretation of
metric geometry as a stable effective representation only when order data has
consistent low-dimensional regularity.
Milestone 14 refines the spatial layer by requiring same-slice comparison.
Radar-time ranks are reconstructed from tick brackets using causal order and
observer tick order. Spatial distance-order comparisons are then restricted to
observer-relative radar-time bins. This avoids treating all accessible events
as if they belonged to one spatial slice and keeps spatial distance explicitly
observer-relative and slice-protocol dependent.
Milestone 15 adds the cross-slice identification layer. Without a transport,
anchor, persistence, or calibration rule, same-position, same-direction,
velocity, constant-velocity, and spatial-evolution questions are undefined
rather than false. When a transport rule is supplied, those become
transport-relative statements. Anchor and persistence experiments show how
additional protocol structure can constrain the per-slice translation,
reflection, orientation, and scale freedoms, while noisy transport degrades the
derived cross-slice quantities.
Milestone 16 identifies weaker transport-gauge relational content. With
supplied persistence labels, pair-distance order histories can record ordinal
shape changes across slices without identifying absolute same positions. These
relational invariants remain weaker than velocity, metric dynamics, or
quantitative spatial evolution, all of which still require transport and
calibration.
Milestone 17 makes persistence itself explicit. Without supplied object labels
or a persistence hypothesis, cross-slice object identity and pair-distance
order histories are undefined. Relational-continuity matching, partial labels,
and anchors can restrict identity hypotheses, but they do not derive object
identity from causal order alone. Different compatible persistence hypotheses
can produce different relational-evolution claims.
Milestone 18 consolidates these layers under a state-change causal trigger
order thesis. Primitive objects are state-changing events, trigger order, and
local finiteness. Observer chains, radar-time ranks, and observer slices are
derived protocol structures. Metric geometry, velocity, curvature values,
seconds, meters, ratios, and quantitative dynamics remain calibrated effective
representations. Quantum compatibility is future work and would require an
additional amplitude or Hilbert-space layer.
Milestone 19 implements a minimal finite trigger graph. Generated trigger
networks are checked as locally finite strict partial orders, and local system
chains are tested as reference ordered protocols without introducing global
physical time.
Milestone 20 ranks reference-chain candidates by coverage, two-sided
bracketing, interval profiles, and protocol-reference choice dependence. The
reference-chain utility score is a finite diagnostic, not a claim of observer
reality.
Milestone 21 computes passive order-level brackets from selected reference
chains. Predecessor and successor reference positions induce radar-time ranks,
bracket-width ranks, and rank slices, all without metric calibration.
Milestone 22 fixes an emission position on a selected reference chain and
computes same-emission echo-return positions and echo-delay ranks. These
rank/order diagnostics study causal response structure relative to a chosen
reference protocol; they do not define physical distance or reconstruct metric
radar distance.
Milestone 23 plants controlled echo-response motifs. Clean motif networks
test whether planted echo-delay ranks are recovered by the echo protocol,
while interference and reference-choice experiments record shortcut returns,
ties, and visibility changes as order-level diagnostics.
Milestone 24 classifies shortcut returns using return spectra. It separates
targeted shortcut-return stress tests from generic background edge
perturbations, reports shortcut count and shortcut depth, and treats the
results as causal-order interference diagnostics.
Milestone 25 studies return-spectrum stability under coarse-graining. It
separates closure-preserving event thinning, immediate-edge thinning, and
reference-chain subsampling, showing which echo classifications survive and
which reflect resolution choices in the finite protocol.
Milestone 26 studies echo-response order signatures over motif populations.
It compares planted, recovered, shortcut-perturbed, coarse-protocol, and
reference-dependent response-rank orders. Stable response-order cores are
treated as inputs for future representability tests, not as metric geometry.

Milestone 27 cleans echo-order semantics before representability checks. It
separates passive `W_rank` from fixed-emission `D_echo`, distinguishes
`S_full`, `S_retained`, and `S_immediate`, compares predeclared gated echo
protocols against the earliest-return rule, and checks whether stable
response-order cores are scalar-rank representable. These are ordinal
precondition diagnostics, not metric reconstruction, clock calibration, or
spatial geometry.

Milestone 28 studies underdetermination of pairwise distance order by scalar
response order. A single reference response signature is not a distance-order
structure. Multiple reference-protocol columns form richer response profiles,
but those profiles are still pre-metric. The representability ladder records
which additional assumptions are needed before pairwise distance-order,
ordinal-embedding, or calibrated metric-representation tests are meaningful.

Milestone 29 defines admissible pairwise response-profile comparison
protocols. These protocols produce response-profile dissimilarities and
response-comparison constraints under declared missing-data policies. They are
pre-metric and do not produce spatial distance. Null baselines are added before
any ordinal embedding attempt.

Milestone 30 validates response-comparison constraint pools. It adds
held-out protocol agreement, bootstrap stability, null-baseline separation,
margin filtering, and coverage diagnostics. It does not perform ordinal
embedding. Passing validation gates makes a pool eligible for future
representability experiments only under explicitly stated assumptions, and
passing gates does not make constraints spatial distances.

Milestone 31 exports predeclared handoff manifests. These manifests are frozen
input specifications for future experiments. They record protocol choices,
validation thresholds, constraint splits, null settings, stop rules, and
forbidden interpretations. They do not contain fitted embeddings and do not
imply distance or geometry. Future representability experiments must report
both eligible and failed manifests.

Milestone 32 performs latent ordinal representation diagnostics from frozen
manifests. The fitted variables are mathematical representation variables, not
spatial coordinates. The manifest train/held-out split is consumed exactly as
exported and is not changed. Held-out performance is a diagnostic of
response-comparison constraint generalization, not physical geometry. Null
baselines test artifacts, not geometry, and failed manifests remain explicit
no-fit rows.

Milestone 33 compares frozen manifest families under preregistered settings.
The comparison does not retune thresholds after seeing fit results.
Target-label permutation is a symmetry control, not a destructive null.
Destructive nulls and marginal-preserving nulls are reported separately.
Failed and ineligible manifests are included in accounting. No physical
geometry is inferred.

Milestone 34 defines carry-forward eligibility for future stress tests.
Families are classified as carry-forward, provisional, blocked, report-only,
or failed-control using fixed diagnostic thresholds. Carry-forward eligibility
is not evidence of metric geometry. Failed, ineligible, blocked, and
provisional families must remain visible in reports. No thresholds are retuned
after seeing fit outcomes.
