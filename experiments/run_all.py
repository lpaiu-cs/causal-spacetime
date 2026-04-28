"""Run the lightweight validation and demonstration experiment suite."""

from __future__ import annotations

from pathlib import Path

import exp02_lorentz_length_contraction as exp02
import exp03_causalset_timelike_reconstruction as exp03
import exp05_finite_speed_lattice_counterexample as exp05
import exp06_spacelike_distance_reconstruction as exp06
import exp12_single_observer_reflection_degeneracy as exp12
from exp07_timelike_pair_reconstruction_convergence import (
    ExperimentConfig,
)
from exp07_timelike_pair_reconstruction_convergence import (
    run_experiment as run_exp07,
)
from exp07_timelike_pair_reconstruction_convergence import (
    save_figures as save_exp07_figures,
)
from exp07_timelike_pair_reconstruction_convergence import (
    write_outputs as write_exp07_outputs,
)
from exp08_probe_pair_statistical_calibration import (
    ExperimentConfig as Exp08Config,
)
from exp08_probe_pair_statistical_calibration import (
    run_experiment as run_exp08,
)
from exp08_probe_pair_statistical_calibration import (
    save_figures as save_exp08_figures,
)
from exp08_probe_pair_statistical_calibration import (
    write_outputs as write_exp08_outputs,
)
from exp10_dimension_reconstruction import (
    ExperimentConfig as Exp10Config,
)
from exp10_dimension_reconstruction import (
    run_experiment as run_exp10,
)
from exp10_dimension_reconstruction import (
    save_figures as save_exp10_figures,
)
from exp10_dimension_reconstruction import (
    write_outputs as write_exp10_outputs,
)
from exp11_discrete_observer_radar_reconstruction import (
    ExperimentConfig as Exp11Config,
)
from exp11_discrete_observer_radar_reconstruction import (
    run_experiment as run_exp11,
)
from exp11_discrete_observer_radar_reconstruction import (
    save_figures as save_exp11_figures,
)
from exp11_discrete_observer_radar_reconstruction import (
    write_outputs as write_exp11_outputs,
)
from exp13_oriented_radar_lorentz_map_recovery import (
    ExperimentConfig as Exp13Config,
)
from exp13_oriented_radar_lorentz_map_recovery import (
    run_experiment as run_exp13,
)
from exp13_oriented_radar_lorentz_map_recovery import (
    save_figures as save_exp13_figures,
)
from exp13_oriented_radar_lorentz_map_recovery import (
    write_outputs as write_exp13_outputs,
)
from exp14_observer_atlas_consistency import (
    ExperimentConfig as Exp14Config,
)
from exp14_observer_atlas_consistency import (
    run_experiment as run_exp14,
)
from exp14_observer_atlas_consistency import (
    save_figures as save_exp14_figures,
)
from exp14_observer_atlas_consistency import (
    write_outputs as write_exp14_outputs,
)
from exp15_exact_poincare_map_sanity import (
    run_experiment as run_exp15,
)
from exp15_exact_poincare_map_sanity import (
    write_outputs as write_exp15_outputs,
)
from exp16_rindler_horizon_reconstruction import (
    ExperimentConfig as Exp16Config,
)
from exp16_rindler_horizon_reconstruction import (
    run_experiment as run_exp16,
)
from exp16_rindler_horizon_reconstruction import (
    save_figures as save_exp16_figures,
)
from exp16_rindler_horizon_reconstruction import (
    write_outputs as write_exp16_outputs,
)
from exp17_inertial_vs_rindler_accessibility import (
    ExperimentConfig as Exp17Config,
)
from exp17_inertial_vs_rindler_accessibility import (
    run_experiment as run_exp17,
)
from exp17_inertial_vs_rindler_accessibility import (
    save_plot as save_exp17_plot,
)
from exp17_inertial_vs_rindler_accessibility import (
    write_outputs as write_exp17_outputs,
)
from exp18_conformal_order_ambiguity import (
    ExperimentConfig as Exp18Config,
)
from exp18_conformal_order_ambiguity import (
    run_experiment as run_exp18,
)
from exp18_conformal_order_ambiguity import (
    save_plot as save_exp18_plot,
)
from exp18_conformal_order_ambiguity import (
    write_outputs as write_exp18_outputs,
)
from exp19_weighted_conformal_volume_reconstruction import (
    ExperimentConfig as Exp19Config,
)
from exp19_weighted_conformal_volume_reconstruction import (
    run_experiment as run_exp19,
)
from exp19_weighted_conformal_volume_reconstruction import (
    save_figures as save_exp19_figures,
)
from exp19_weighted_conformal_volume_reconstruction import (
    write_outputs as write_exp19_outputs,
)
from exp20_conformal_volume_exact_sanity import (
    run_experiment as run_exp20,
)
from exp20_conformal_volume_exact_sanity import (
    write_outputs as write_exp20_outputs,
)
from exp21_physical_measure_sprinkling import (
    ExperimentConfig as Exp21Config,
)
from exp21_physical_measure_sprinkling import (
    run_experiment as run_exp21,
)
from exp21_physical_measure_sprinkling import (
    save_figures as save_exp21_figures,
)
from exp21_physical_measure_sprinkling import (
    write_outputs as write_exp21_outputs,
)
from exp22_local_measure_profile_estimation import (
    ExperimentConfig as Exp22Config,
)
from exp22_local_measure_profile_estimation import (
    run_experiment as run_exp22,
)
from exp22_local_measure_profile_estimation import (
    save_figures as save_exp22_figures,
)
from exp22_local_measure_profile_estimation import (
    write_outputs as write_exp22_outputs,
)
from exp23_thinning_coarse_graining_stability import (
    ExperimentConfig as Exp23Config,
)
from exp23_thinning_coarse_graining_stability import (
    run_experiment as run_exp23,
)
from exp23_thinning_coarse_graining_stability import (
    save_figures as save_exp23_figures,
)
from exp23_thinning_coarse_graining_stability import (
    write_outputs as write_exp23_outputs,
)
from exp24_measure_sprinkling_exact_sanity import (
    run_experiment as run_exp24,
)
from exp24_measure_sprinkling_exact_sanity import (
    write_outputs as write_exp24_outputs,
)
from exp25_radar_return_distance_order import (
    ExperimentConfig as Exp25Config,
)
from exp25_radar_return_distance_order import (
    run_experiment as run_exp25,
)
from exp25_radar_return_distance_order import (
    save_figures as save_exp25_figures,
)
from exp25_radar_return_distance_order import (
    write_outputs as write_exp25_outputs,
)
from exp26_metric_representation_scale_invariance import (
    run_experiment as run_exp26,
)
from exp26_metric_representation_scale_invariance import (
    save_plot as save_exp26_plot,
)
from exp26_metric_representation_scale_invariance import (
    write_outputs as write_exp26_outputs,
)
from exp27_ratio_stability_from_calibration import (
    run_experiment as run_exp27,
)
from exp27_ratio_stability_from_calibration import (
    save_plot as save_exp27_plot,
)
from exp27_ratio_stability_from_calibration import (
    write_outputs as write_exp27_outputs,
)
from exp28_oriented_chart_distance_order_preservation import (
    ExperimentConfig as Exp28Config,
)
from exp28_oriented_chart_distance_order_preservation import (
    run_experiment as run_exp28,
)
from exp28_oriented_chart_distance_order_preservation import (
    save_plot as save_exp28_plot,
)
from exp28_oriented_chart_distance_order_preservation import (
    write_outputs as write_exp28_outputs,
)
from exp29_metric_representability_diagnostics import (
    run_experiment as run_exp29,
)
from exp29_metric_representability_diagnostics import (
    write_outputs as write_exp29_outputs,
)
from exp30_ordinal_exact_sanity import (
    run_experiment as run_exp30,
)
from exp30_ordinal_exact_sanity import (
    write_outputs as write_exp30_outputs,
)
from exp31_ordinal_embedding_recovery import (
    ExperimentConfig as Exp31Config,
)
from exp31_ordinal_embedding_recovery import (
    run_experiment as run_exp31,
)
from exp31_ordinal_embedding_recovery import (
    save_figures as save_exp31_figures,
)
from exp31_ordinal_embedding_recovery import (
    write_outputs as write_exp31_outputs,
)
from exp32_ordinal_dimension_selection import (
    ExperimentConfig as Exp32Config,
)
from exp32_ordinal_dimension_selection import (
    run_experiment as run_exp32,
)
from exp32_ordinal_dimension_selection import (
    save_figures as save_exp32_figures,
)
from exp32_ordinal_dimension_selection import (
    write_outputs as write_exp32_outputs,
)
from exp33_noisy_incomplete_order_embedding import (
    ExperimentConfig as Exp33Config,
)
from exp33_noisy_incomplete_order_embedding import (
    run_experiment as run_exp33,
)
from exp33_noisy_incomplete_order_embedding import (
    save_figures as save_exp33_figures,
)
from exp33_noisy_incomplete_order_embedding import (
    write_outputs as write_exp33_outputs,
)
from exp34_observer_distance_order_embedding import (
    ExperimentConfig as Exp34Config,
)
from exp34_observer_distance_order_embedding import (
    run_experiment as run_exp34,
)
from exp34_observer_distance_order_embedding import (
    save_figures as save_exp34_figures,
)
from exp34_observer_distance_order_embedding import (
    write_outputs as write_exp34_outputs,
)
from exp35_ordinal_embedding_exact_sanity import (
    run_experiment as run_exp35,
)
from exp35_ordinal_embedding_exact_sanity import (
    write_outputs as write_exp35_outputs,
)
from exp36_heldout_ordinal_embedding_validation import (
    ExperimentConfig as Exp36Config,
)
from exp36_heldout_ordinal_embedding_validation import (
    run_experiment as run_exp36,
)
from exp36_heldout_ordinal_embedding_validation import (
    save_figures as save_exp36_figures,
)
from exp36_heldout_ordinal_embedding_validation import (
    write_outputs as write_exp36_outputs,
)
from exp37_embedding_stability_under_subsampling import (
    ExperimentConfig as Exp37Config,
)
from exp37_embedding_stability_under_subsampling import (
    run_experiment as run_exp37,
)
from exp37_embedding_stability_under_subsampling import (
    save_figures as save_exp37_figures,
)
from exp37_embedding_stability_under_subsampling import (
    write_outputs as write_exp37_outputs,
)
from exp38_observer_order_null_baseline import (
    ExperimentConfig as Exp38Config,
)
from exp38_observer_order_null_baseline import (
    run_experiment as run_exp38,
)
from exp38_observer_order_null_baseline import (
    save_figures as save_exp38_figures,
)
from exp38_observer_order_null_baseline import (
    write_outputs as write_exp38_outputs,
)
from exp39_effective_metric_complexity_curve import (
    ExperimentConfig as Exp39Config,
)
from exp39_effective_metric_complexity_curve import (
    run_experiment as run_exp39,
)
from exp39_effective_metric_complexity_curve import (
    save_figures as save_exp39_figures,
)
from exp39_effective_metric_complexity_curve import (
    write_outputs as write_exp39_outputs,
)
from exp40_representation_stability_exact_sanity import (
    run_experiment as run_exp40,
)
from exp40_representation_stability_exact_sanity import (
    write_outputs as write_exp40_outputs,
)
from exp41_radar_time_order_from_brackets import (
    ExperimentConfig as Exp41Config,
)
from exp41_radar_time_order_from_brackets import (
    run_experiment as run_exp41,
)
from exp41_radar_time_order_from_brackets import (
    save_figure as save_exp41_figure,
)
from exp41_radar_time_order_from_brackets import (
    write_outputs as write_exp41_outputs,
)
from exp42_same_slice_distance_order_preservation import (
    ExperimentConfig as Exp42Config,
)
from exp42_same_slice_distance_order_preservation import (
    run_experiment as run_exp42,
)
from exp42_same_slice_distance_order_preservation import (
    save_figures as save_exp42_figures,
)
from exp42_same_slice_distance_order_preservation import (
    write_outputs as write_exp42_outputs,
)
from exp43_sliced_observer_order_null_baseline import (
    ExperimentConfig as Exp43Config,
)
from exp43_sliced_observer_order_null_baseline import (
    run_experiment as run_exp43,
)
from exp43_sliced_observer_order_null_baseline import (
    save_figures as save_exp43_figures,
)
from exp43_sliced_observer_order_null_baseline import (
    write_outputs as write_exp43_outputs,
)
from exp44_slice_width_sensitivity import (
    ExperimentConfig as Exp44Config,
)
from exp44_slice_width_sensitivity import (
    run_experiment as run_exp44,
)
from exp44_slice_width_sensitivity import (
    save_figures as save_exp44_figures,
)
from exp44_slice_width_sensitivity import (
    write_outputs as write_exp44_outputs,
)
from exp45_spatial_slice_exact_sanity import (
    run_experiment as run_exp45,
)
from exp45_spatial_slice_exact_sanity import (
    write_outputs as write_exp45_outputs,
)
from exp46_cross_slice_predicate_undefined import (
    run_experiment as run_exp46,
)
from exp46_cross_slice_predicate_undefined import (
    write_outputs as write_exp46_outputs,
)
from exp47_sliced_constraint_graph_decomposition import (
    ExperimentConfig as Exp47Config,
)
from exp47_sliced_constraint_graph_decomposition import (
    run_experiment as run_exp47,
)
from exp47_sliced_constraint_graph_decomposition import (
    save_figures as save_exp47_figures,
)
from exp47_sliced_constraint_graph_decomposition import (
    write_outputs as write_exp47_outputs,
)
from exp48_slice_local_embedding_validation import (
    ExperimentConfig as Exp48Config,
)
from exp48_slice_local_embedding_validation import (
    run_experiment as run_exp48,
)
from exp48_slice_local_embedding_validation import (
    save_figures as save_exp48_figures,
)
from exp48_slice_local_embedding_validation import (
    write_outputs as write_exp48_outputs,
)
from exp49_slice_gauge_dependence import (
    ExperimentConfig as Exp49Config,
)
from exp49_slice_gauge_dependence import (
    run_experiment as run_exp49,
)
from exp49_slice_gauge_dependence import (
    save_figures as save_exp49_figures,
)
from exp49_slice_gauge_dependence import (
    write_outputs as write_exp49_outputs,
)
from exp50_anchor_constrained_transport import (
    ExperimentConfig as Exp50Config,
)
from exp50_anchor_constrained_transport import (
    run_experiment as run_exp50,
)
from exp50_anchor_constrained_transport import (
    save_figures as save_exp50_figures,
)
from exp50_anchor_constrained_transport import (
    write_outputs as write_exp50_outputs,
)
from exp51_persistence_dependent_velocity import (
    run_experiment as run_exp51,
)
from exp51_persistence_dependent_velocity import (
    save_figure as save_exp51_figure,
)
from exp51_persistence_dependent_velocity import (
    write_outputs as write_exp51_outputs,
)
from exp52_noisy_transport_sensitivity import (
    ExperimentConfig as Exp52Config,
)
from exp52_noisy_transport_sensitivity import (
    run_experiment as run_exp52,
)
from exp52_noisy_transport_sensitivity import (
    save_figures as save_exp52_figures,
)
from exp52_noisy_transport_sensitivity import (
    write_outputs as write_exp52_outputs,
)
from exp53_cross_slice_transport_exact_sanity import (
    run_experiment as run_exp53,
)
from exp53_cross_slice_transport_exact_sanity import (
    write_outputs as write_exp53_outputs,
)
from exp54_predicate_definability_table import (
    run_experiment as run_exp54,
)
from exp54_predicate_definability_table import (
    write_outputs as write_exp54_outputs,
)
from exp55_relational_shape_history_without_transport import (
    ExperimentConfig as Exp55Config,
)
from exp55_relational_shape_history_without_transport import (
    run_experiment as run_exp55,
)
from exp55_relational_shape_history_without_transport import (
    save_figure as save_exp55_figure,
)
from exp55_relational_shape_history_without_transport import (
    write_outputs as write_exp55_outputs,
)
from exp56_relational_history_gauge_invariance import (
    ExperimentConfig as Exp56Config,
)
from exp56_relational_history_gauge_invariance import (
    run_experiment as run_exp56,
)
from exp56_relational_history_gauge_invariance import (
    save_figure as save_exp56_figure,
)
from exp56_relational_history_gauge_invariance import (
    write_outputs as write_exp56_outputs,
)
from exp57_observer_slice_relational_evolution import (
    ExperimentConfig as Exp57Config,
)
from exp57_observer_slice_relational_evolution import (
    run_experiment as run_exp57,
)
from exp57_observer_slice_relational_evolution import (
    save_figure as save_exp57_figure,
)
from exp57_observer_slice_relational_evolution import (
    write_outputs as write_exp57_outputs,
)
from exp58_relational_invariants_vs_velocity import (
    ExperimentConfig as Exp58Config,
)
from exp58_relational_invariants_vs_velocity import (
    run_experiment as run_exp58,
)
from exp58_relational_invariants_vs_velocity import (
    save_figure as save_exp58_figure,
)
from exp58_relational_invariants_vs_velocity import (
    write_outputs as write_exp58_outputs,
)
from exp59_relational_evolution_exact_sanity import (
    run_experiment as run_exp59,
)
from exp59_relational_evolution_exact_sanity import (
    write_outputs as write_exp59_outputs,
)
from exp60_persistence_predicate_undefined import (
    run_experiment as run_exp60,
)
from exp60_persistence_predicate_undefined import (
    write_outputs as write_exp60_outputs,
)
from exp61_symmetric_persistence_ambiguity import (
    ExperimentConfig as Exp61Config,
)
from exp61_symmetric_persistence_ambiguity import (
    run_experiment as run_exp61,
)
from exp61_symmetric_persistence_ambiguity import (
    save_figure as save_exp61_figure,
)
from exp61_symmetric_persistence_ambiguity import (
    write_outputs as write_exp61_outputs,
)
from exp62_relational_persistence_matching_recovery import (
    ExperimentConfig as Exp62Config,
)
from exp62_relational_persistence_matching_recovery import (
    run_experiment as run_exp62,
)
from exp62_relational_persistence_matching_recovery import (
    save_figures as save_exp62_figures,
)
from exp62_relational_persistence_matching_recovery import (
    write_outputs as write_exp62_outputs,
)
from exp63_partial_label_constrained_persistence import (
    ExperimentConfig as Exp63Config,
)
from exp63_partial_label_constrained_persistence import (
    run_experiment as run_exp63,
)
from exp63_partial_label_constrained_persistence import (
    save_figures as save_exp63_figures,
)
from exp63_partial_label_constrained_persistence import (
    write_outputs as write_exp63_outputs,
)
from exp64_crossing_persistence_failure import (
    ExperimentConfig as Exp64Config,
)
from exp64_crossing_persistence_failure import (
    run_experiment as run_exp64,
)
from exp64_crossing_persistence_failure import (
    save_figure as save_exp64_figure,
)
from exp64_crossing_persistence_failure import (
    write_outputs as write_exp64_outputs,
)
from exp65_persistence_hypothesis_dependence import (
    ExperimentConfig as Exp65Config,
)
from exp65_persistence_hypothesis_dependence import (
    run_experiment as run_exp65,
)
from exp65_persistence_hypothesis_dependence import (
    save_figure as save_exp65_figure,
)
from exp65_persistence_hypothesis_dependence import (
    write_outputs as write_exp65_outputs,
)
from exp66_persistence_matching_exact_sanity import (
    run_experiment as run_exp66,
)
from exp66_persistence_matching_exact_sanity import (
    write_outputs as write_exp66_outputs,
)
from exp67_state_change_exact_sanity import (
    run_experiment as run_exp67,
)
from exp67_state_change_exact_sanity import (
    write_outputs as write_exp67_outputs,
)
from exp68_state_change_toy_model import (
    ExperimentConfig as Exp68Config,
)
from exp68_state_change_toy_model import (
    run_experiment as run_exp68,
)
from exp68_state_change_toy_model import (
    save_figures as save_exp68_figures,
)
from exp68_state_change_toy_model import (
    write_outputs as write_exp68_outputs,
)
from exp69_state_change_observer_chain_diagnostic import (
    ExperimentConfig as Exp69Config,
)
from exp69_state_change_observer_chain_diagnostic import (
    run_experiment as run_exp69,
)
from exp69_state_change_observer_chain_diagnostic import (
    save_figure as save_exp69_figure,
)
from exp69_state_change_observer_chain_diagnostic import (
    write_outputs as write_exp69_outputs,
)
from exp70_observer_chain_exact_sanity import (
    run_experiment as run_exp70,
)
from exp70_observer_chain_exact_sanity import (
    write_outputs as write_exp70_outputs,
)
from exp71_observer_chain_candidate_ranking import (
    ExperimentConfig as Exp71Config,
)
from exp71_observer_chain_candidate_ranking import (
    run_experiment as run_exp71,
)
from exp71_observer_chain_candidate_ranking import (
    save_figures as save_exp71_figures,
)
from exp71_observer_chain_candidate_ranking import (
    write_outputs as write_exp71_outputs,
)
from exp72_observer_chain_coverage_vs_trigger_probability import (
    ExperimentConfig as Exp72Config,
)
from exp72_observer_chain_coverage_vs_trigger_probability import (
    run_experiment as run_exp72,
)
from exp72_observer_chain_coverage_vs_trigger_probability import (
    save_figures as save_exp72_figures,
)
from exp72_observer_chain_coverage_vs_trigger_probability import (
    write_outputs as write_exp72_outputs,
)
from exp73_observer_chain_interval_profile import (
    ExperimentConfig as Exp73Config,
)
from exp73_observer_chain_interval_profile import (
    run_experiment as run_exp73,
)
from exp73_observer_chain_interval_profile import (
    save_figures as save_exp73_figures,
)
from exp73_observer_chain_interval_profile import (
    write_outputs as write_exp73_outputs,
)
from exp74_state_change_reference_bracket_diagnostics import (
    ExperimentConfig as Exp74Config,
)
from exp74_state_change_reference_bracket_diagnostics import (
    run_experiment as run_exp74,
)
from exp74_state_change_reference_bracket_diagnostics import (
    save_figures as save_exp74_figures,
)
from exp74_state_change_reference_bracket_diagnostics import (
    write_outputs as write_exp74_outputs,
)
from exp75_state_change_bracket_rank_reference_dependence import (
    ExperimentConfig as Exp75Config,
)
from exp75_state_change_bracket_rank_reference_dependence import (
    run_experiment as run_exp75,
)
from exp75_state_change_bracket_rank_reference_dependence import (
    save_figures as save_exp75_figures,
)
from exp75_state_change_bracket_rank_reference_dependence import (
    write_outputs as write_exp75_outputs,
)
from exp76_state_change_reference_bracket_coverage_vs_trigger_density import (
    ExperimentConfig as Exp76Config,
)
from exp76_state_change_reference_bracket_coverage_vs_trigger_density import (
    run_experiment as run_exp76,
)
from exp76_state_change_reference_bracket_coverage_vs_trigger_density import (
    save_figures as save_exp76_figures,
)
from exp76_state_change_reference_bracket_coverage_vs_trigger_density import (
    write_outputs as write_exp76_outputs,
)
from exp77_state_change_reference_bracket_exact_sanity import (
    run_experiment as run_exp77,
)
from exp77_state_change_reference_bracket_exact_sanity import (
    write_outputs as write_exp77_outputs,
)
from exp78_state_change_echo_exact_sanity import (
    run_experiment as run_exp78,
)
from exp78_state_change_echo_exact_sanity import (
    write_outputs as write_exp78_outputs,
)
from exp79_state_change_echo_order_diagnostics import (
    ExperimentConfig as Exp79Config,
)
from exp79_state_change_echo_order_diagnostics import (
    run_experiment as run_exp79,
)
from exp79_state_change_echo_order_diagnostics import (
    save_figures as save_exp79_figures,
)
from exp79_state_change_echo_order_diagnostics import (
    write_outputs as write_exp79_outputs,
)
from exp80_state_change_echo_reference_dependence import (
    ExperimentConfig as Exp80Config,
)
from exp80_state_change_echo_reference_dependence import (
    run_experiment as run_exp80,
)
from exp80_state_change_echo_reference_dependence import (
    save_figures as save_exp80_figures,
)
from exp80_state_change_echo_reference_dependence import (
    write_outputs as write_exp80_outputs,
)
from exp81_state_change_echo_emission_sensitivity import (
    ExperimentConfig as Exp81Config,
)
from exp81_state_change_echo_emission_sensitivity import (
    run_experiment as run_exp81,
)
from exp81_state_change_echo_emission_sensitivity import (
    save_figures as save_exp81_figures,
)
from exp81_state_change_echo_emission_sensitivity import (
    write_outputs as write_exp81_outputs,
)
from exp82_state_change_echo_coverage_vs_trigger_density import (
    ExperimentConfig as Exp82Config,
)
from exp82_state_change_echo_coverage_vs_trigger_density import (
    run_experiment as run_exp82,
)
from exp82_state_change_echo_coverage_vs_trigger_density import (
    save_figures as save_exp82_figures,
)
from exp82_state_change_echo_coverage_vs_trigger_density import (
    write_outputs as write_exp82_outputs,
)
from exp83_echo_motif_exact_sanity import (
    run_experiment as run_exp83,
)
from exp83_echo_motif_exact_sanity import (
    write_outputs as write_exp83_outputs,
)
from exp84_planted_echo_motif_recovery import (
    ExperimentConfig as Exp84Config,
)
from exp84_planted_echo_motif_recovery import (
    run_experiment as run_exp84,
)
from exp84_planted_echo_motif_recovery import (
    save_figures as save_exp84_figures,
)
from exp84_planted_echo_motif_recovery import (
    write_outputs as write_exp84_outputs,
)
from exp85_echo_motif_background_interference import (
    ExperimentConfig as Exp85Config,
)
from exp85_echo_motif_background_interference import (
    run_experiment as run_exp85,
)
from exp85_echo_motif_background_interference import (
    save_figures as save_exp85_figures,
)
from exp85_echo_motif_background_interference import (
    write_outputs as write_exp85_outputs,
)
from exp86_echo_motif_density_resolution import (
    ExperimentConfig as Exp86Config,
)
from exp86_echo_motif_density_resolution import (
    run_experiment as run_exp86,
)
from exp86_echo_motif_density_resolution import (
    save_figures as save_exp86_figures,
)
from exp86_echo_motif_density_resolution import (
    write_outputs as write_exp86_outputs,
)
from exp87_echo_motif_reference_choice_visibility import (
    ExperimentConfig as Exp87Config,
)
from exp87_echo_motif_reference_choice_visibility import (
    run_experiment as run_exp87,
)
from exp87_echo_motif_reference_choice_visibility import (
    save_figures as save_exp87_figures,
)
from exp87_echo_motif_reference_choice_visibility import (
    write_outputs as write_exp87_outputs,
)
from exp88_echo_shortcut_exact_sanity import (
    run_experiment as run_exp88,
)
from exp88_echo_shortcut_exact_sanity import (
    write_outputs as write_exp88_outputs,
)
from exp89_echo_shortcut_injection_sweep import (
    ExperimentConfig as Exp89Config,
)
from exp89_echo_shortcut_injection_sweep import (
    run_experiment as run_exp89,
)
from exp89_echo_shortcut_injection_sweep import (
    save_figures as save_exp89_figures,
)
from exp89_echo_shortcut_injection_sweep import (
    write_outputs as write_exp89_outputs,
)
from exp90_echo_background_edge_perturbation import (
    ExperimentConfig as Exp90Config,
)
from exp90_echo_background_edge_perturbation import (
    run_experiment as run_exp90,
)
from exp90_echo_background_edge_perturbation import (
    save_figures as save_exp90_figures,
)
from exp90_echo_background_edge_perturbation import (
    write_outputs as write_exp90_outputs,
)
from exp91_echo_motif_path_length_robustness import (
    ExperimentConfig as Exp91Config,
)
from exp91_echo_motif_path_length_robustness import (
    run_experiment as run_exp91,
)
from exp91_echo_motif_path_length_robustness import (
    save_figures as save_exp91_figures,
)
from exp91_echo_motif_path_length_robustness import (
    write_outputs as write_exp91_outputs,
)
from exp92_echo_shortcut_reference_dependence import (
    ExperimentConfig as Exp92Config,
)
from exp92_echo_shortcut_reference_dependence import (
    run_experiment as run_exp92,
)
from exp92_echo_shortcut_reference_dependence import (
    save_figures as save_exp92_figures,
)
from exp92_echo_shortcut_reference_dependence import (
    write_outputs as write_exp92_outputs,
)
from exp93_echo_interference_exact_sanity import (
    run_experiment as run_exp93,
)
from exp93_echo_interference_exact_sanity import (
    write_outputs as write_exp93_outputs,
)
from exp94_echo_coarse_graining_exact_sanity import (
    run_experiment as run_exp94,
)
from exp94_echo_coarse_graining_exact_sanity import (
    write_outputs as write_exp94_outputs,
)
from exp95_echo_event_thinning_stability import (
    ExperimentConfig as Exp95Config,
)
from exp95_echo_event_thinning_stability import (
    run_experiment as run_exp95,
)
from exp95_echo_event_thinning_stability import (
    save_figures as save_exp95_figures,
)
from exp95_echo_event_thinning_stability import (
    write_outputs as write_exp95_outputs,
)
from exp96_echo_reference_subsampling_resolution import (
    ExperimentConfig as Exp96Config,
)
from exp96_echo_reference_subsampling_resolution import (
    run_experiment as run_exp96,
)
from exp96_echo_reference_subsampling_resolution import (
    save_figures as save_exp96_figures,
)
from exp96_echo_reference_subsampling_resolution import (
    write_outputs as write_exp96_outputs,
)
from exp97_echo_edge_thinning_fragility import (
    ExperimentConfig as Exp97Config,
)
from exp97_echo_edge_thinning_fragility import (
    run_experiment as run_exp97,
)
from exp97_echo_edge_thinning_fragility import (
    save_figures as save_exp97_figures,
)
from exp97_echo_edge_thinning_fragility import (
    write_outputs as write_exp97_outputs,
)
from exp98_echo_shortcut_classification_under_coarse_graining import (
    ExperimentConfig as Exp98Config,
)
from exp98_echo_shortcut_classification_under_coarse_graining import (
    run_experiment as run_exp98,
)
from exp98_echo_shortcut_classification_under_coarse_graining import (
    save_figures as save_exp98_figures,
)
from exp98_echo_shortcut_classification_under_coarse_graining import (
    write_outputs as write_exp98_outputs,
)
from exp99_echo_return_spectrum_stability_exact_sanity import (
    run_experiment as run_exp99,
)
from exp99_echo_return_spectrum_stability_exact_sanity import (
    write_outputs as write_exp99_outputs,
)


def run_lorentz_length_contraction() -> None:
    rows = exp02.run_experiment()
    summary_path = exp02.write_summary(rows)
    figure_path = exp02.save_plot(rows)
    print(f"exp02 wrote {summary_path} and {figure_path}")


def run_legacy_timelike_reconstruction() -> None:
    rows = exp03.run_experiment()
    summary_path = exp03.write_summary(rows)
    figure_path = exp03.save_timelike_reconstruction_plot(rows, exp03.OUTPUT_FIGURE)
    print(f"exp03 wrote {summary_path} and {figure_path}")


def run_finite_speed_lattice_counterexample() -> None:
    rows, lattice_events, sprinkled_events = exp05.run_experiment()
    summary_path = exp05.write_summary(rows)
    cone_path = exp05.save_cone_plot(lattice_events, sprinkled_events)
    growth_path = exp05.save_growth_plot(rows)
    exp05.print_direction_summary()
    print(f"exp05 wrote {summary_path}, {cone_path}, and {growth_path}")


def run_spacelike_proxy_experiment() -> None:
    rows = exp06.run_experiment()
    summary_path = exp06.write_summary(rows)
    figure_path = exp06.save_proxy_plot(rows)
    print(f"exp06 wrote {summary_path} and {figure_path}")


def run_timelike_pair_convergence() -> None:
    config = ExperimentConfig(
        T=2.0,
        n_values=(300, 600),
        repetitions=2,
        pairs_per_repetition=40,
        seed=0,
        max_pairs_for_chain=10,
        skip_chain=False,
        output_dir=Path("outputs"),
    )
    pair_rows, summary_rows = run_exp07(config)
    pair_path, summary_path = write_exp07_outputs(
        pair_rows,
        summary_rows,
        config.output_dir,
    )
    scatter_path, error_path = save_exp07_figures(
        pair_rows,
        summary_rows,
        config.output_dir,
    )
    print(f"exp07 wrote {pair_path}, {summary_path}, {scatter_path}, and {error_path}")


def run_probe_pair_statistical_calibration() -> None:
    config = Exp08Config(
        T=2.0,
        n_values=(300, 600),
        repetitions=2,
        pairs_per_repetition=80,
        seed=0,
        min_tau=0.10,
        max_tau=None,
        output_dir=Path("outputs"),
    )
    pair_rows, summary_rows, binned_rows = run_exp08(config)
    pair_path, summary_path, binned_path = write_exp08_outputs(
        pair_rows,
        summary_rows,
        binned_rows,
        config.output_dir,
    )
    figure_paths = save_exp08_figures(
        pair_rows,
        summary_rows,
        binned_rows,
        config.output_dir,
    )
    print(
        "exp08 wrote "
        f"{pair_path}, {summary_path}, {binned_path}, and {len(figure_paths)} figures"
    )


def run_dimension_reconstruction() -> None:
    config = Exp10Config(
        spacetime_dims=(2, 3),
        n_values=(300, 600),
        repetitions=2,
        T=2.0,
        seed=0,
        output_dir=Path("outputs"),
    )
    result_rows, summary_rows = run_exp10(config)
    result_path, summary_path = write_exp10_outputs(
        result_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp10_figures(result_rows, summary_rows, config)
    print(
        "exp10 wrote "
        f"{result_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_discrete_observer_radar_reconstruction() -> None:
    config = Exp11Config(
        T=2.0,
        n_values=(300,),
        tick_values=(16, 32),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
        use_extracted_chain=False,
    )
    event_rows, summary_rows = run_exp11(config)
    event_path, summary_path = write_exp11_outputs(
        event_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp11_figures(event_rows, summary_rows, config.output_dir)
    print(
        "exp11 wrote "
        f"{event_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_single_observer_reflection_degeneracy() -> None:
    rows = exp12.run_experiment()
    summary_path = exp12.write_summary(rows)
    figure_path = exp12.save_plot(rows)
    print(f"exp12 wrote {summary_path} and {figure_path}")


def run_oriented_radar_lorentz_map_recovery() -> None:
    config = Exp13Config(
        T=2.0,
        n_values=(300,),
        tick_values=(32, 64),
        beta_values=(0.3,),
        beacon_separation=0.15,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    event_rows, summary_rows = run_exp13(config)
    event_path, summary_path = write_exp13_outputs(
        event_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp13_figures(event_rows, summary_rows, config.output_dir)
    print(
        "exp13 wrote "
        f"{event_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_observer_atlas_consistency() -> None:
    config = Exp14Config(
        T=2.0,
        n_values=(300,),
        tick_values=(32, 64),
        repetitions=2,
        seed=0,
        beacon_separation=0.15,
        output_dir=Path("outputs"),
        num_invariant_pairs=200,
    )
    chart_rows, transition_rows, loop_rows = run_exp14(config)
    chart_path, transition_path, loop_path = write_exp14_outputs(
        chart_rows,
        transition_rows,
        loop_rows,
        config.output_dir,
    )
    figure_paths = save_exp14_figures(
        transition_rows,
        loop_rows,
        config.output_dir,
    )
    print(
        "exp14 wrote "
        f"{chart_path}, {transition_path}, {loop_path}, and "
        f"{len(figure_paths)} figures"
    )


def run_exact_poincare_map_sanity() -> None:
    rows = run_exp15()
    output_path = write_exp15_outputs(rows)
    print(f"exp15 wrote {output_path}")


def run_rindler_horizon_reconstruction() -> None:
    config = Exp16Config(
        T=4.0,
        acceleration_values=(2.0,),
        n_values=(600,),
        tick_values=(32, 64),
        repetitions=2,
        seed=0,
        direction=1,
        output_dir=Path("outputs"),
    )
    event_rows, summary_rows = run_exp16(config)
    event_path, summary_path = write_exp16_outputs(
        event_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp16_figures(event_rows, summary_rows, config)
    print(
        "exp16 wrote "
        f"{event_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_inertial_vs_rindler_accessibility() -> None:
    config = Exp17Config(
        T=4.0,
        n_events=400,
        tick_count=64,
        acceleration=2.0,
        seed=0,
        direction=1,
        output_dir=Path("outputs"),
    )
    rows = run_exp17(config)
    data_path = write_exp17_outputs(
        rows,
        config.output_dir / "data" / "inertial_vs_rindler_accessibility.csv",
    )
    figure_path = save_exp17_plot(
        rows,
        config.output_dir / "figures" / "inertial_vs_rindler_accessibility.png",
    )
    print(f"exp17 wrote {data_path} and {figure_path}")


def run_conformal_order_ambiguity() -> None:
    config = Exp18Config(output_dir=Path("outputs"))
    rows = run_exp18(config)
    summary_path = write_exp18_outputs(rows, config.output_dir)
    figure_path = save_exp18_plot(rows, config.output_dir)
    print(f"exp18 wrote {summary_path} and {figure_path}")


def run_weighted_conformal_volume_reconstruction() -> None:
    config = Exp19Config(
        T=2.0,
        n_values=(600,),
        repetitions=2,
        pairs_per_repetition=100,
        seed=0,
        output_dir=Path("outputs"),
    )
    pair_rows, summary_rows = run_exp19(config)
    pair_path, summary_path = write_exp19_outputs(
        pair_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp19_figures(pair_rows, summary_rows, config.output_dir)
    print(
        "exp19 wrote "
        f"{pair_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_conformal_volume_exact_sanity() -> None:
    rows = run_exp20()
    output_path = write_exp20_outputs(rows)
    print(f"exp20 wrote {output_path}")


def run_physical_measure_sprinkling() -> None:
    config = Exp21Config(
        T=2.0,
        n_values=(600,),
        repetitions=2,
        pairs_per_repetition=100,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows, summary = run_exp21(config)
    pair_path, summary_path = write_exp21_outputs(rows, summary, config.output_dir)
    figure_paths = save_exp21_figures(rows, summary, config.output_dir)
    print(
        "exp21 wrote "
        f"{pair_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_local_measure_profile_estimation() -> None:
    config = Exp22Config(
        T=2.0,
        n_values=(600, 1200),
        repetitions=2,
        num_bins=12,
        seed=0,
        output_dir=Path("outputs"),
    )
    bin_rows, summary_rows = run_exp22(config)
    bin_path, summary_path = write_exp22_outputs(
        bin_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp22_figures(bin_rows, summary_rows, config.output_dir)
    print(
        "exp22 wrote "
        f"{bin_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_thinning_coarse_graining_stability() -> None:
    config = Exp23Config(
        T=2.0,
        n_events=1200,
        repetitions=2,
        keep_probabilities=(1.0, 0.5),
        pairs_per_repetition=100,
        seed=0,
        output_dir=Path("outputs"),
    )
    pair_rows, summary_rows = run_exp23(config)
    pair_path, summary_path = write_exp23_outputs(
        pair_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp23_figures(summary_rows, config)
    print(
        "exp23 wrote "
        f"{pair_path}, {summary_path}, and {len(figure_paths)} figures"
    )


def run_measure_sprinkling_exact_sanity() -> None:
    rows = run_exp24()
    output_path = write_exp24_outputs(rows)
    print(f"exp24 wrote {output_path}")


def run_radar_return_distance_order() -> None:
    config = Exp25Config(
        T=2.0,
        tick_values=(16, 32),
        target_counts=(50,),
        repetitions=2,
        seed=0,
        emission_tick_fraction=0.25,
        output_dir=Path("outputs"),
    )
    rows = run_exp25(config)
    data_path = write_exp25_outputs(rows, config.output_dir)
    figure_paths = save_exp25_figures(rows, config.output_dir)
    print(f"exp25 wrote {data_path} and {len(figure_paths)} figures")


def run_metric_representation_scale_invariance() -> None:
    rows = run_exp26()
    data_path = write_exp26_outputs(rows)
    figure_path = save_exp26_plot(rows)
    print(f"exp26 wrote {data_path} and {figure_path}")


def run_ratio_stability_from_calibration() -> None:
    rows = run_exp27()
    data_path = write_exp27_outputs(rows)
    figure_path = save_exp27_plot(rows)
    print(f"exp27 wrote {data_path} and {figure_path}")


def run_oriented_chart_distance_order_preservation() -> None:
    config = Exp28Config(
        T=2.0,
        n_values=(300,),
        tick_values=(32, 64),
        repetitions=2,
        seed=0,
        beacon_separation=0.15,
        pair_count=200,
        output_dir=Path("outputs"),
    )
    rows = run_exp28(config)
    data_path = write_exp28_outputs(rows, config.output_dir)
    figure_path = save_exp28_plot(rows, config.output_dir)
    print(f"exp28 wrote {data_path} and {figure_path}")


def run_metric_representability_diagnostics() -> None:
    rows = run_exp29()
    output_path = write_exp29_outputs(rows)
    print(f"exp29 wrote {output_path}")


def run_ordinal_exact_sanity() -> None:
    rows = run_exp30()
    output_path = write_exp30_outputs(rows)
    print(f"exp30 wrote {output_path}")


def run_ordinal_embedding_recovery() -> None:
    config = Exp31Config(
        true_dims=(1,),
        n_points_values=(20,),
        constraint_counts=(1000,),
        repetitions=2,
        seed=0,
        steps=600,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp31(config)
    data_path = write_exp31_outputs(rows, config.output_dir)
    figure_paths = save_exp31_figures(rows, config.output_dir)
    print(f"exp31 wrote {data_path} and {len(figure_paths)} figures")


def run_ordinal_dimension_selection() -> None:
    config = Exp32Config(
        true_dims=(1, 2),
        candidate_dims=(1, 2, 3),
        n_points=25,
        num_constraints=1500,
        repetitions=2,
        seed=0,
        steps=600,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp32(config)
    data_path = write_exp32_outputs(rows, config.output_dir)
    figure_paths = save_exp32_figures(rows, config.output_dir)
    print(f"exp32 wrote {data_path} and {len(figure_paths)} figures")


def run_noisy_incomplete_order_embedding() -> None:
    config = Exp33Config(
        true_dim=2,
        n_points=40,
        constraint_counts=(1000,),
        flip_probabilities=(0.0, 0.05),
        repetitions=2,
        seed=0,
        steps=600,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp33(config)
    data_path = write_exp33_outputs(rows, config.output_dir)
    figure_paths = save_exp33_figures(rows, config.output_dir)
    print(f"exp33 wrote {data_path} and {len(figure_paths)} figures")


def run_observer_distance_order_embedding() -> None:
    config = Exp34Config(
        T=2.0,
        n_events=250,
        tick_values=(32, 64),
        constraint_counts=(1000,),
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        steps=600,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp34(config)
    data_path = write_exp34_outputs(rows, config.output_dir)
    figure_paths = save_exp34_figures(rows, config.output_dir)
    print(f"exp34 wrote {data_path} and {len(figure_paths)} figures")


def run_ordinal_embedding_exact_sanity() -> None:
    rows = run_exp35()
    output_path = write_exp35_outputs(rows)
    print(f"exp35 wrote {output_path}")


def run_heldout_ordinal_embedding_validation() -> None:
    config = Exp36Config(
        true_dim=2,
        n_points_values=(30,),
        constraint_counts=(1000,),
        repetitions=2,
        seed=0,
        steps=400,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp36(config)
    data_path = write_exp36_outputs(rows, config.output_dir)
    figure_paths = save_exp36_figures(rows, config.output_dir)
    print(f"exp36 wrote {data_path} and {len(figure_paths)} figures")


def run_embedding_stability_under_subsampling() -> None:
    config = Exp37Config(
        true_dim=2,
        n_points=40,
        total_constraints=3000,
        subset_sizes=(500, 1000),
        num_subsets=3,
        repetitions=2,
        seed=0,
        steps=400,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp37(config)
    data_path = write_exp37_outputs(rows, config.output_dir)
    figure_paths = save_exp37_figures(rows, config.output_dir)
    print(f"exp37 wrote {data_path} and {len(figure_paths)} figures")


def run_observer_order_null_baseline() -> None:
    config = Exp38Config(
        T=2.0,
        n_events=250,
        tick_values=(32, 64),
        constraint_counts=(1000,),
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        steps=400,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp38(config)
    data_path = write_exp38_outputs(rows, config.output_dir)
    figure_paths = save_exp38_figures(rows, config.output_dir)
    print(f"exp38 wrote {data_path} and {len(figure_paths)} figures")


def run_effective_metric_complexity_curve() -> None:
    config = Exp39Config(
        true_dim=2,
        candidate_dims=(1, 2, 3),
        n_points=30,
        num_constraints=1500,
        repetitions=2,
        seed=0,
        steps=400,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp39(config)
    data_path = write_exp39_outputs(rows, config.output_dir)
    figure_paths = save_exp39_figures(rows, config.output_dir)
    print(f"exp39 wrote {data_path} and {len(figure_paths)} figures")


def run_representation_stability_exact_sanity() -> None:
    rows = run_exp40()
    output_path = write_exp40_outputs(rows)
    print(f"exp40 wrote {output_path}")


def run_radar_time_order_from_brackets() -> None:
    config = Exp41Config(
        T=2.0,
        n_values=(300,),
        tick_values=(16, 32),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp41(config)
    data_path = write_exp41_outputs(rows, config.output_dir)
    figure_path = save_exp41_figure(rows, config.output_dir)
    print(f"exp41 wrote {data_path} and {figure_path}")


def run_same_slice_distance_order_preservation() -> None:
    config = Exp42Config(
        T=2.0,
        n_values=(300,),
        tick_values=(32, 64),
        bin_width_values=(2,),
        repetitions=2,
        seed=0,
        beacon_separation=0.15,
        max_pairs_per_slice=100,
        output_dir=Path("outputs"),
    )
    rows = run_exp42(config)
    data_path = write_exp42_outputs(rows, config.output_dir)
    figure_paths = save_exp42_figures(rows, config.output_dir)
    print(f"exp42 wrote {data_path} and {len(figure_paths)} figures")


def run_sliced_observer_order_null_baseline() -> None:
    config = Exp43Config(
        T=2.0,
        n_events=300,
        tick_values=(32, 64),
        bin_width_values=(2,),
        constraint_counts=(500,),
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        steps=400,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp43(config)
    data_path = write_exp43_outputs(rows, config.output_dir)
    figure_paths = save_exp43_figures(rows, config.output_dir)
    print(f"exp43 wrote {data_path} and {len(figure_paths)} figures")


def run_slice_width_sensitivity() -> None:
    config = Exp44Config(
        T=2.0,
        n_events=400,
        tick_count=64,
        bin_width_values=(1, 2, 4),
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp44(config)
    data_path = write_exp44_outputs(rows, config.output_dir)
    figure_paths = save_exp44_figures(rows, config.output_dir)
    print(f"exp44 wrote {data_path} and {len(figure_paths)} figures")


def run_spatial_slice_exact_sanity() -> None:
    rows = run_exp45()
    output_path = write_exp45_outputs(rows)
    print(f"exp45 wrote {output_path}")


def run_cross_slice_predicate_undefined() -> None:
    rows = run_exp46()
    output_path = write_exp46_outputs(
        rows,
        Path("outputs/data/cross_slice_predicate_undefined.csv"),
    )
    print(f"exp46 wrote {output_path}")


def run_sliced_constraint_graph_decomposition() -> None:
    config = Exp47Config(
        T=2.0,
        n_events=300,
        tick_count=64,
        bin_width_values=(2, 4),
        constraint_count=1000,
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp47(config)
    data_path = write_exp47_outputs(rows, config.output_dir)
    figure_paths = save_exp47_figures(rows, config.output_dir)
    print(f"exp47 wrote {data_path} and {len(figure_paths)} figures")


def run_slice_local_embedding_validation() -> None:
    config = Exp48Config(
        T=2.0,
        n_values=(300,),
        tick_values=(64,),
        bin_width_values=(2,),
        constraint_count=1000,
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        steps=300,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp48(config)
    data_path = write_exp48_outputs(rows, config.output_dir)
    figure_paths = save_exp48_figures(rows, config.output_dir)
    print(f"exp48 wrote {data_path} and {len(figure_paths)} figures")


def run_slice_gauge_dependence() -> None:
    config = Exp49Config(
        T=2.0,
        n_events=400,
        tick_count=64,
        bin_width=4,
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp49(config)
    data_path = write_exp49_outputs(rows, config.output_dir)
    figure_paths = save_exp49_figures(rows, config.output_dir)
    print(f"exp49 wrote {data_path} and {len(figure_paths)} figures")


def run_anchor_constrained_transport() -> None:
    config = Exp50Config(
        T=2.0,
        n_events=400,
        tick_count=64,
        bin_width=4,
        constraint_count=1000,
        repetitions=2,
        beacon_separation=0.15,
        seed=0,
        steps=300,
        restarts=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp50(config)
    data_path = write_exp50_outputs(rows, config.output_dir)
    figure_paths = save_exp50_figures(rows, config.output_dir)
    print(f"exp50 wrote {data_path} and {len(figure_paths)} figures")


def run_persistence_dependent_velocity() -> None:
    rows = run_exp51()
    data_path = write_exp51_outputs(rows, Path("outputs"))
    figure_path = save_exp51_figure(rows, Path("outputs"))
    print(f"exp51 wrote {data_path} and {figure_path}")


def run_noisy_transport_sensitivity() -> None:
    config = Exp52Config(
        T=2.0,
        n_events=400,
        tick_count=64,
        bin_width=4,
        anchor_positions=(-0.35, 0.0, 0.35),
        noise_levels=(0.0, 0.05),
        repetitions=2,
        seed=0,
        beacon_separation=0.15,
        output_dir=Path("outputs"),
    )
    rows = run_exp52(config)
    data_path = write_exp52_outputs(rows, config.output_dir)
    figure_paths = save_exp52_figures(rows, config.output_dir)
    print(f"exp52 wrote {data_path} and {len(figure_paths)} figures")


def run_cross_slice_transport_exact_sanity() -> None:
    rows = run_exp53()
    output_path = write_exp53_outputs(rows)
    print(f"exp53 wrote {output_path}")


def run_predicate_definability_table() -> None:
    rows = run_exp54()
    output_path = write_exp54_outputs(rows, Path("outputs"))
    print(f"exp54 wrote {output_path}")


def run_relational_shape_history_without_transport() -> None:
    config = Exp55Config(
        slice_count=12,
        object_count=5,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp55(config)
    data_path = write_exp55_outputs(rows, config.output_dir)
    figure_path = save_exp55_figure(rows, config.output_dir)
    print(f"exp55 wrote {data_path} and {figure_path}")


def run_relational_history_gauge_invariance() -> None:
    config = Exp56Config(
        slice_count=8,
        object_count=5,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp56(config)
    data_path = write_exp56_outputs(rows, config.output_dir)
    figure_path = save_exp56_figure(rows, config.output_dir)
    print(f"exp56 wrote {data_path} and {figure_path}")


def run_observer_slice_relational_evolution() -> None:
    config = Exp57Config(
        T=2.0,
        slice_count=8,
        object_count=5,
        tick_count=64,
        bin_width=4,
        beacon_separation=0.15,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp57(config)
    data_path = write_exp57_outputs(rows, config.output_dir)
    figure_path = save_exp57_figure(rows, config.output_dir)
    print(f"exp57 wrote {data_path} and {figure_path}")


def run_relational_invariants_vs_velocity() -> None:
    config = Exp58Config(
        slice_count=12,
        object_count=5,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp58(config)
    data_path = write_exp58_outputs(rows, config.output_dir)
    figure_path = save_exp58_figure(rows, config.output_dir)
    print(f"exp58 wrote {data_path} and {figure_path}")


def run_relational_evolution_exact_sanity() -> None:
    rows = run_exp59()
    output_path = write_exp59_outputs(rows)
    print(f"exp59 wrote {output_path}")


def run_persistence_predicate_undefined() -> None:
    rows = run_exp60()
    output_path = write_exp60_outputs(rows, Path("outputs"))
    print(f"exp60 wrote {output_path}")


def run_symmetric_persistence_ambiguity() -> None:
    config = Exp61Config(
        slice_count=3,
        object_counts=(4,),
        spacing=1.0,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp61(config)
    data_path = write_exp61_outputs(rows, config.output_dir)
    figure_path = save_exp61_figure(rows, config.output_dir)
    print(f"exp61 wrote {data_path} and {figure_path}")


def run_relational_persistence_matching_recovery() -> None:
    config = Exp62Config(
        slice_count=5,
        object_counts=(4,),
        motion_scales=(0.05, 0.10),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp62(config)
    data_path = write_exp62_outputs(rows, config.output_dir)
    figure_paths = save_exp62_figures(rows, config.output_dir)
    print(f"exp62 wrote {data_path} and {len(figure_paths)} figures")


def run_partial_label_constrained_persistence() -> None:
    config = Exp63Config(
        slice_count=5,
        object_count=5,
        known_fractions=(0.0, 0.5),
        motion_scale=0.10,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp63(config)
    data_path = write_exp63_outputs(rows, config.output_dir)
    figure_paths = save_exp63_figures(rows, config.output_dir)
    print(f"exp63 wrote {data_path} and {len(figure_paths)} figures")


def run_crossing_persistence_failure() -> None:
    config = Exp64Config(
        slice_count=6,
        object_count=5,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp64(config)
    data_path = write_exp64_outputs(rows, config.output_dir)
    figure_path = save_exp64_figure(rows, config.output_dir)
    print(f"exp64 wrote {data_path} and {figure_path}")


def run_persistence_hypothesis_dependence() -> None:
    config = Exp65Config(
        slice_count=6,
        object_count=5,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp65(config)
    data_path = write_exp65_outputs(rows, config.output_dir)
    figure_path = save_exp65_figure(rows, config.output_dir)
    print(f"exp65 wrote {data_path} and {figure_path}")


def run_persistence_matching_exact_sanity() -> None:
    rows = run_exp66()
    output_path = write_exp66_outputs(rows)
    print(f"exp66 wrote {output_path}")


def run_state_change_exact_sanity() -> None:
    rows = run_exp67()
    output_path = write_exp67_outputs(rows)
    print(f"exp67 wrote {output_path}")


def run_state_change_toy_model() -> None:
    config = Exp68Config(
        num_systems_values=(5,),
        max_events_values=(100,),
        trigger_probability_values=(0.20,),
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp68(config)
    data_path = write_exp68_outputs(rows, config.output_dir)
    figure_paths = save_exp68_figures(rows, config.output_dir)
    print(f"exp68 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_observer_chain_diagnostic() -> None:
    config = Exp69Config(
        num_systems=5,
        max_events=100,
        trigger_probability=0.20,
        max_triggers_per_event=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp69(config)
    data_path = write_exp69_outputs(rows, config.output_dir)
    figure_path = save_exp69_figure(rows, config.output_dir)
    print(f"exp69 wrote {data_path} and {figure_path}")


def run_observer_chain_exact_sanity() -> None:
    rows = run_exp70()
    output_path = write_exp70_outputs(rows)
    print(f"exp70 wrote {output_path}")


def run_observer_chain_candidate_ranking() -> None:
    config = Exp71Config(
        num_systems_values=(5,),
        max_events_values=(100,),
        trigger_probability_values=(0.20,),
        max_triggers_per_event=2,
        repetitions=2,
        random_candidate_count=3,
        seed=0,
        output_dir=Path("outputs"),
    )
    ranking_rows, ambiguity_rows = run_exp71(config)
    ranking_path, ambiguity_path = write_exp71_outputs(
        ranking_rows,
        ambiguity_rows,
        config.output_dir,
    )
    figure_paths = save_exp71_figures(
        ranking_rows,
        ambiguity_rows,
        config.output_dir,
    )
    print(
        f"exp71 wrote {ranking_path}, {ambiguity_path}, "
        f"and {len(figure_paths)} figures"
    )


def run_observer_chain_coverage_vs_trigger_probability() -> None:
    config = Exp72Config(
        num_systems=5,
        max_events=150,
        trigger_probability_values=(0.10, 0.30),
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp72(config)
    data_path = write_exp72_outputs(rows, config.output_dir)
    figure_paths = save_exp72_figures(rows, config.output_dir)
    print(f"exp72 wrote {data_path} and {len(figure_paths)} figures")


def run_observer_chain_interval_profile() -> None:
    config = Exp73Config(
        num_systems=5,
        max_events=150,
        trigger_probability=0.20,
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows, example_profile = run_exp73(config)
    data_path = write_exp73_outputs(rows, config.output_dir)
    figure_paths = save_exp73_figures(rows, example_profile, config.output_dir)
    print(f"exp73 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_reference_bracket_diagnostics() -> None:
    config = Exp74Config(
        num_systems_values=(5,),
        max_events_values=(100,),
        trigger_probability_values=(0.20,),
        max_triggers_per_event=2,
        repetitions=2,
        random_candidate_count=3,
        seed=0,
        bin_width=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp74(config)
    data_path = write_exp74_outputs(rows, config.output_dir)
    figure_paths = save_exp74_figures(rows, config.output_dir)
    print(f"exp74 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_bracket_rank_reference_dependence() -> None:
    config = Exp75Config(
        num_systems=5,
        max_events=150,
        trigger_probability_values=(0.20,),
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        top_k=3,
        bin_width=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp75(config)
    data_path = write_exp75_outputs(rows, config.output_dir)
    figure_paths = save_exp75_figures(rows, config.output_dir)
    print(f"exp75 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_reference_bracket_coverage_vs_trigger_density() -> None:
    config = Exp76Config(
        num_systems=5,
        max_events=150,
        trigger_probability_values=(0.10, 0.30),
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        bin_width=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp76(config)
    data_path = write_exp76_outputs(rows, config.output_dir)
    figure_paths = save_exp76_figures(rows, config.output_dir)
    print(f"exp76 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_reference_bracket_exact_sanity() -> None:
    rows = run_exp77()
    output_path = write_exp77_outputs(rows)
    print(f"exp77 wrote {output_path}")


def run_state_change_echo_exact_sanity() -> None:
    rows = run_exp78()
    output_path = write_exp78_outputs(rows)
    print(f"exp78 wrote {output_path}")


def run_state_change_echo_order_diagnostics() -> None:
    config = Exp79Config(
        num_systems_values=(5,),
        max_events_values=(100,),
        trigger_probability_values=(0.20,),
        max_triggers_per_event=2,
        repetitions=2,
        random_candidate_count=3,
        emission_count=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp79(config)
    data_path = write_exp79_outputs(rows, config.output_dir)
    figure_paths = save_exp79_figures(rows, config.output_dir)
    print(f"exp79 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_echo_reference_dependence() -> None:
    config = Exp80Config(
        num_systems=5,
        max_events=150,
        trigger_probability_values=(0.20,),
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        top_k=3,
        output_dir=Path("outputs"),
    )
    rows = run_exp80(config)
    data_path = write_exp80_outputs(rows, config.output_dir)
    figure_paths = save_exp80_figures(rows, config.output_dir)
    print(f"exp80 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_echo_emission_sensitivity() -> None:
    config = Exp81Config(
        num_systems=5,
        max_events=150,
        trigger_probability=0.20,
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        emission_count=3,
        output_dir=Path("outputs"),
    )
    rows = run_exp81(config)
    data_path = write_exp81_outputs(rows, config.output_dir)
    figure_paths = save_exp81_figures(rows, config.output_dir)
    print(f"exp81 wrote {data_path} and {len(figure_paths)} figures")


def run_state_change_echo_coverage_vs_trigger_density() -> None:
    config = Exp82Config(
        num_systems=5,
        max_events=150,
        trigger_probability_values=(0.10, 0.30),
        max_triggers_per_event=2,
        repetitions=2,
        seed=0,
        emission_count=2,
        output_dir=Path("outputs"),
    )
    rows = run_exp82(config)
    data_path = write_exp82_outputs(rows, config.output_dir)
    figure_paths = save_exp82_figures(rows, config.output_dir)
    print(f"exp82 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_motif_exact_sanity() -> None:
    rows = run_exp83()
    output_path = write_exp83_outputs(rows)
    print(f"exp83 wrote {output_path}")


def run_planted_echo_motif_recovery() -> None:
    config = Exp84Config(
        reference_lengths=(16,),
        motif_counts=(10,),
        delay_ranks=(2, 3, 5),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    motif_rows, summary_rows = run_exp84(config)
    motif_path, summary_path = write_exp84_outputs(
        motif_rows,
        summary_rows,
        config.output_dir,
    )
    figure_paths = save_exp84_figures(summary_rows, config.output_dir)
    print(
        f"exp84 wrote {motif_path}, {summary_path}, "
        f"and {len(figure_paths)} figures"
    )


def run_echo_motif_background_interference() -> None:
    config = Exp85Config(
        num_systems=5,
        max_events=150,
        trigger_probability_values=(0.05, 0.20),
        motif_count=10,
        delay_ranks=(2, 3, 5),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp85(config)
    data_path = write_exp85_outputs(rows, config.output_dir)
    figure_paths = save_exp85_figures(rows, config.output_dir)
    print(f"exp85 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_motif_density_resolution() -> None:
    config = Exp86Config(
        reference_length=32,
        motif_counts=(10, 30),
        delay_rank_sets=("small", "medium"),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp86(config)
    data_path = write_exp86_outputs(rows, config.output_dir)
    figure_paths = save_exp86_figures(rows, config.output_dir)
    print(f"exp86 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_motif_reference_choice_visibility() -> None:
    config = Exp87Config(
        num_systems=5,
        max_events=150,
        trigger_probability=0.20,
        motif_count=10,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp87(config)
    data_path = write_exp87_outputs(rows, config.output_dir)
    figure_paths = save_exp87_figures(rows, config.output_dir)
    print(f"exp87 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_shortcut_exact_sanity() -> None:
    rows = run_exp88()
    output_path = write_exp88_outputs(rows)
    print(f"exp88 wrote {output_path}")


def run_echo_shortcut_injection_sweep() -> None:
    config = Exp89Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(3, 5, 8),
        shortcut_probabilities=(0.0, 0.3),
        shortcut_modes=("target_to_early_reference",),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp89(config)
    data_path = write_exp89_outputs(rows, config.output_dir)
    figure_paths = save_exp89_figures(rows, config.output_dir)
    print(f"exp89 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_background_edge_perturbation() -> None:
    config = Exp90Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(3, 5, 8),
        edge_probabilities=(0.0, 0.01),
        repetitions=2,
        seed=0,
        max_edges=100,
        output_dir=Path("outputs"),
    )
    rows = run_exp90(config)
    data_path = write_exp90_outputs(rows, config.output_dir)
    figure_paths = save_exp90_figures(rows, config.output_dir)
    print(f"exp90 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_motif_path_length_robustness() -> None:
    config = Exp91Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(5, 8),
        outward_steps_values=(0, 1),
        return_steps_values=(0, 1),
        shortcut_probability=0.3,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp91(config)
    data_path = write_exp91_outputs(rows, config.output_dir)
    figure_paths = save_exp91_figures(rows, config.output_dir)
    print(f"exp91 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_shortcut_reference_dependence() -> None:
    config = Exp92Config(
        num_systems=5,
        max_events=150,
        trigger_probability=0.20,
        motif_count=10,
        shortcut_probability=0.3,
        repetitions=2,
        seed=0,
        top_k=3,
        output_dir=Path("outputs"),
    )
    rows = run_exp92(config)
    data_path = write_exp92_outputs(rows, config.output_dir)
    figure_paths = save_exp92_figures(rows, config.output_dir)
    print(f"exp92 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_interference_exact_sanity() -> None:
    rows = run_exp93()
    output_path = write_exp93_outputs(rows)
    print(f"exp93 wrote {output_path}")


def run_echo_coarse_graining_exact_sanity() -> None:
    rows = run_exp94()
    output_path = write_exp94_outputs(rows)
    print(f"exp94 wrote {output_path}")


def run_echo_event_thinning_stability() -> None:
    config = Exp95Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(3, 5, 8),
        keep_probabilities=(1.0, 0.4),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp95(config)
    data_path = write_exp95_outputs(rows, config.output_dir)
    figure_paths = save_exp95_figures(rows, config.output_dir)
    print(f"exp95 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_reference_subsampling_resolution() -> None:
    config = Exp96Config(
        reference_length=64,
        motif_count=30,
        delay_ranks=(2, 3, 5, 8),
        strides=(1, 4),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp96(config)
    data_path = write_exp96_outputs(rows, config.output_dir)
    figure_paths = save_exp96_figures(rows, config.output_dir)
    print(f"exp96 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_edge_thinning_fragility() -> None:
    config = Exp97Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(3, 5, 8),
        removal_probabilities=(0.0, 0.3),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp97(config)
    data_path = write_exp97_outputs(rows, config.output_dir)
    figure_paths = save_exp97_figures(rows, config.output_dir)
    print(f"exp97 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_shortcut_classification_under_coarse_graining() -> None:
    config = Exp98Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(3, 5, 8),
        shortcut_probability=0.3,
        event_keep_probability=0.4,
        edge_removal_probability=0.15,
        reference_strides=(1, 4),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp98(config)
    data_path = write_exp98_outputs(rows, config.output_dir)
    figure_paths = save_exp98_figures(rows, config.output_dir)
    print(f"exp98 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_return_spectrum_stability_exact_sanity() -> None:
    rows = run_exp99()
    output_path = write_exp99_outputs(rows)
    print(f"exp99 wrote {output_path}")


def main() -> None:
    run_lorentz_length_contraction()
    run_legacy_timelike_reconstruction()
    run_finite_speed_lattice_counterexample()
    run_spacelike_proxy_experiment()
    run_timelike_pair_convergence()
    run_probe_pair_statistical_calibration()
    run_dimension_reconstruction()
    run_discrete_observer_radar_reconstruction()
    run_single_observer_reflection_degeneracy()
    run_oriented_radar_lorentz_map_recovery()
    run_observer_atlas_consistency()
    run_exact_poincare_map_sanity()
    run_rindler_horizon_reconstruction()
    run_inertial_vs_rindler_accessibility()
    run_conformal_order_ambiguity()
    run_weighted_conformal_volume_reconstruction()
    run_conformal_volume_exact_sanity()
    run_physical_measure_sprinkling()
    run_local_measure_profile_estimation()
    run_thinning_coarse_graining_stability()
    run_measure_sprinkling_exact_sanity()
    run_radar_return_distance_order()
    run_metric_representation_scale_invariance()
    run_ratio_stability_from_calibration()
    run_oriented_chart_distance_order_preservation()
    run_metric_representability_diagnostics()
    run_ordinal_exact_sanity()
    run_ordinal_embedding_recovery()
    run_ordinal_dimension_selection()
    run_noisy_incomplete_order_embedding()
    run_observer_distance_order_embedding()
    run_ordinal_embedding_exact_sanity()
    run_heldout_ordinal_embedding_validation()
    run_embedding_stability_under_subsampling()
    run_observer_order_null_baseline()
    run_effective_metric_complexity_curve()
    run_representation_stability_exact_sanity()
    run_radar_time_order_from_brackets()
    run_same_slice_distance_order_preservation()
    run_sliced_observer_order_null_baseline()
    run_slice_width_sensitivity()
    run_spatial_slice_exact_sanity()
    run_cross_slice_predicate_undefined()
    run_sliced_constraint_graph_decomposition()
    run_slice_local_embedding_validation()
    run_slice_gauge_dependence()
    run_anchor_constrained_transport()
    run_persistence_dependent_velocity()
    run_noisy_transport_sensitivity()
    run_cross_slice_transport_exact_sanity()
    run_predicate_definability_table()
    run_relational_shape_history_without_transport()
    run_relational_history_gauge_invariance()
    run_observer_slice_relational_evolution()
    run_relational_invariants_vs_velocity()
    run_relational_evolution_exact_sanity()
    run_persistence_predicate_undefined()
    run_symmetric_persistence_ambiguity()
    run_relational_persistence_matching_recovery()
    run_partial_label_constrained_persistence()
    run_crossing_persistence_failure()
    run_persistence_hypothesis_dependence()
    run_persistence_matching_exact_sanity()
    run_state_change_exact_sanity()
    run_state_change_toy_model()
    run_state_change_observer_chain_diagnostic()
    run_observer_chain_exact_sanity()
    run_observer_chain_candidate_ranking()
    run_observer_chain_coverage_vs_trigger_probability()
    run_observer_chain_interval_profile()
    run_state_change_reference_bracket_diagnostics()
    run_state_change_bracket_rank_reference_dependence()
    run_state_change_reference_bracket_coverage_vs_trigger_density()
    run_state_change_reference_bracket_exact_sanity()
    run_state_change_echo_exact_sanity()
    run_state_change_echo_order_diagnostics()
    run_state_change_echo_reference_dependence()
    run_state_change_echo_emission_sensitivity()
    run_state_change_echo_coverage_vs_trigger_density()
    run_echo_motif_exact_sanity()
    run_planted_echo_motif_recovery()
    run_echo_motif_background_interference()
    run_echo_motif_density_resolution()
    run_echo_motif_reference_choice_visibility()
    run_echo_shortcut_exact_sanity()
    run_echo_shortcut_injection_sweep()
    run_echo_background_edge_perturbation()
    run_echo_motif_path_length_robustness()
    run_echo_shortcut_reference_dependence()
    run_echo_interference_exact_sanity()
    run_echo_coarse_graining_exact_sanity()
    run_echo_event_thinning_stability()
    run_echo_reference_subsampling_resolution()
    run_echo_edge_thinning_fragility()
    run_echo_shortcut_classification_under_coarse_graining()
    run_echo_return_spectrum_stability_exact_sanity()


if __name__ == "__main__":
    main()
