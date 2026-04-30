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
from exp100_echo_response_signature_exact_sanity import (
    run_experiment as run_exp100,
)
from exp100_echo_response_signature_exact_sanity import (
    write_outputs as write_exp100_outputs,
)
from exp101_layered_echo_response_order_recovery import (
    ExperimentConfig as Exp101Config,
)
from exp101_layered_echo_response_order_recovery import (
    run_experiment as run_exp101,
)
from exp101_layered_echo_response_order_recovery import (
    save_figures as save_exp101_figures,
)
from exp101_layered_echo_response_order_recovery import (
    write_outputs as write_exp101_outputs,
)
from exp102_echo_response_signature_coarse_protocol_stability import (
    ExperimentConfig as Exp102Config,
)
from exp102_echo_response_signature_coarse_protocol_stability import (
    run_experiment as run_exp102,
)
from exp102_echo_response_signature_coarse_protocol_stability import (
    save_figures as save_exp102_figures,
)
from exp102_echo_response_signature_coarse_protocol_stability import (
    write_outputs as write_exp102_outputs,
)
from exp103_echo_response_shortcut_robust_core import (
    ExperimentConfig as Exp103Config,
)
from exp103_echo_response_shortcut_robust_core import (
    run_experiment as run_exp103,
)
from exp103_echo_response_shortcut_robust_core import (
    save_figures as save_exp103_figures,
)
from exp103_echo_response_shortcut_robust_core import (
    write_outputs as write_exp103_outputs,
)
from exp104_echo_response_reference_protocol_dependence import (
    ExperimentConfig as Exp104Config,
)
from exp104_echo_response_reference_protocol_dependence import (
    run_experiment as run_exp104,
)
from exp104_echo_response_reference_protocol_dependence import (
    save_figures as save_exp104_figures,
)
from exp104_echo_response_reference_protocol_dependence import (
    write_outputs as write_exp104_outputs,
)
from exp105_echo_response_order_precondition_diagnostics import (
    ExperimentConfig as Exp105Config,
)
from exp105_echo_response_order_precondition_diagnostics import (
    run_experiment as run_exp105,
)
from exp105_echo_response_order_precondition_diagnostics import (
    save_figures as save_exp105_figures,
)
from exp105_echo_response_order_precondition_diagnostics import (
    write_outputs as write_exp105_outputs,
)
from exp106_echo_response_signature_stability_exact_sanity import (
    run_experiment as run_exp106,
)
from exp106_echo_response_signature_stability_exact_sanity import (
    write_outputs as write_exp106_outputs,
)
from exp107_echo_order_semantics_exact_sanity import (
    run_experiment as run_exp107,
)
from exp107_echo_order_semantics_exact_sanity import (
    write_outputs as write_exp107_outputs,
)
from exp108_gated_echo_protocol_comparison import (
    ExperimentConfig as Exp108Config,
)
from exp108_gated_echo_protocol_comparison import (
    run_experiment as run_exp108,
)
from exp108_gated_echo_protocol_comparison import (
    save_figures as save_exp108_figures,
)
from exp108_gated_echo_protocol_comparison import (
    write_outputs as write_exp108_outputs,
)
from exp109_echo_delay_spacing_tie_resolution import (
    ExperimentConfig as Exp109Config,
)
from exp109_echo_delay_spacing_tie_resolution import (
    run_experiment as run_exp109,
)
from exp109_echo_delay_spacing_tie_resolution import (
    save_figures as save_exp109_figures,
)
from exp109_echo_delay_spacing_tie_resolution import (
    write_outputs as write_exp109_outputs,
)
from exp110_response_order_scalar_representability import (
    ExperimentConfig as Exp110Config,
)
from exp110_response_order_scalar_representability import (
    run_experiment as run_exp110,
)
from exp110_response_order_scalar_representability import (
    save_figures as save_exp110_figures,
)
from exp110_response_order_scalar_representability import (
    write_outputs as write_exp110_outputs,
)
from exp111_echo_terminology_audit import (
    run_experiment as run_exp111,
)
from exp111_echo_terminology_audit import (
    write_outputs as write_exp111_outputs,
)
from exp112_response_representability_exact_sanity import (
    run_experiment as run_exp112,
)
from exp112_response_representability_exact_sanity import (
    write_outputs as write_exp112_outputs,
)
from exp113_response_order_underdetermination_exact_sanity import (
    run_experiment as run_exp113,
)
from exp113_response_order_underdetermination_exact_sanity import (
    write_outputs as write_exp113_outputs,
)
from exp114_single_reference_response_order_underdetermination import (
    ExperimentConfig as Exp114Config,
)
from exp114_single_reference_response_order_underdetermination import (
    run_experiment as run_exp114,
)
from exp114_single_reference_response_order_underdetermination import (
    save_figures as save_exp114_figures,
)
from exp114_single_reference_response_order_underdetermination import (
    write_outputs as write_exp114_outputs,
)
from exp115_multi_reference_response_profile_diagnostics import (
    ExperimentConfig as Exp115Config,
)
from exp115_multi_reference_response_profile_diagnostics import (
    run_experiment as run_exp115,
)
from exp115_multi_reference_response_profile_diagnostics import (
    save_figures as save_exp115_figures,
)
from exp115_multi_reference_response_profile_diagnostics import (
    write_outputs as write_exp115_outputs,
)
from exp116_response_representability_requirement_table import (
    run_experiment as run_exp116,
)
from exp116_response_representability_requirement_table import (
    write_outputs as write_exp116_outputs,
)
from exp117_response_profile_stability_under_protocol_variants import (
    ExperimentConfig as Exp117Config,
)
from exp117_response_profile_stability_under_protocol_variants import (
    run_experiment as run_exp117,
)
from exp117_response_profile_stability_under_protocol_variants import (
    save_figures as save_exp117_figures,
)
from exp117_response_profile_stability_under_protocol_variants import (
    write_outputs as write_exp117_outputs,
)
from exp118_response_profile_requirement_exact_sanity import (
    run_experiment as run_exp118,
)
from exp118_response_profile_requirement_exact_sanity import (
    write_outputs as write_exp118_outputs,
)
from exp119_pairwise_response_profile_exact_sanity import (
    run_experiment as run_exp119,
)
from exp119_pairwise_response_profile_exact_sanity import (
    write_outputs as write_exp119_outputs,
)
from exp120_pairwise_response_protocol_comparison import (
    ExperimentConfig as Exp120Config,
)
from exp120_pairwise_response_protocol_comparison import (
    run_experiment as run_exp120,
)
from exp120_pairwise_response_protocol_comparison import (
    save_figures as save_exp120_figures,
)
from exp120_pairwise_response_protocol_comparison import (
    write_outputs as write_exp120_outputs,
)
from exp121_pairwise_response_null_baselines import (
    ExperimentConfig as Exp121Config,
)
from exp121_pairwise_response_null_baselines import (
    run_experiment as run_exp121,
)
from exp121_pairwise_response_null_baselines import (
    save_figures as save_exp121_figures,
)
from exp121_pairwise_response_null_baselines import (
    write_outputs as write_exp121_outputs,
)
from exp122_pairwise_response_protocol_variant_stability import (
    ExperimentConfig as Exp122Config,
)
from exp122_pairwise_response_protocol_variant_stability import (
    run_experiment as run_exp122,
)
from exp122_pairwise_response_protocol_variant_stability import (
    save_figures as save_exp122_figures,
)
from exp122_pairwise_response_protocol_variant_stability import (
    write_outputs as write_exp122_outputs,
)
from exp123_pairwise_response_missing_data_sensitivity import (
    ExperimentConfig as Exp123Config,
)
from exp123_pairwise_response_missing_data_sensitivity import (
    run_experiment as run_exp123,
)
from exp123_pairwise_response_missing_data_sensitivity import (
    save_figures as save_exp123_figures,
)
from exp123_pairwise_response_missing_data_sensitivity import (
    write_outputs as write_exp123_outputs,
)
from exp124_pairwise_response_protocol_choice_dependence import (
    ExperimentConfig as Exp124Config,
)
from exp124_pairwise_response_protocol_choice_dependence import (
    run_experiment as run_exp124,
)
from exp124_pairwise_response_protocol_choice_dependence import (
    save_figures as save_exp124_figures,
)
from exp124_pairwise_response_protocol_choice_dependence import (
    write_outputs as write_exp124_outputs,
)
from exp125_pairwise_response_null_admissibility_exact_sanity import (
    run_experiment as run_exp125,
)
from exp125_pairwise_response_null_admissibility_exact_sanity import (
    write_outputs as write_exp125_outputs,
)
from exp126_response_constraint_pool_exact_sanity import (
    run_experiment as run_exp126,
)
from exp126_response_constraint_pool_exact_sanity import (
    write_outputs as write_exp126_outputs,
)
from exp127_response_constraint_heldout_protocol_validation import (
    ExperimentConfig as Exp127Config,
)
from exp127_response_constraint_heldout_protocol_validation import (
    run_experiment as run_exp127,
)
from exp127_response_constraint_heldout_protocol_validation import (
    save_figures as save_exp127_figures,
)
from exp127_response_constraint_heldout_protocol_validation import (
    write_outputs as write_exp127_outputs,
)
from exp128_response_constraint_bootstrap_stability import (
    ExperimentConfig as Exp128Config,
)
from exp128_response_constraint_bootstrap_stability import (
    run_experiment as run_exp128,
)
from exp128_response_constraint_bootstrap_stability import (
    save_figures as save_exp128_figures,
)
from exp128_response_constraint_bootstrap_stability import (
    write_outputs as write_exp128_outputs,
)
from exp129_response_constraint_null_separation import (
    ExperimentConfig as Exp129Config,
)
from exp129_response_constraint_null_separation import (
    run_experiment as run_exp129,
)
from exp129_response_constraint_null_separation import (
    save_figures as save_exp129_figures,
)
from exp129_response_constraint_null_separation import (
    write_outputs as write_exp129_outputs,
)
from exp130_response_constraint_pool_coverage import (
    ExperimentConfig as Exp130Config,
)
from exp130_response_constraint_pool_coverage import (
    run_experiment as run_exp130,
)
from exp130_response_constraint_pool_coverage import (
    save_figures as save_exp130_figures,
)
from exp130_response_constraint_pool_coverage import (
    write_outputs as write_exp130_outputs,
)
from exp131_response_constraint_validation_gate_summary import (
    ExperimentConfig as Exp131Config,
)
from exp131_response_constraint_validation_gate_summary import (
    run_experiment as run_exp131,
)
from exp131_response_constraint_validation_gate_summary import (
    save_figures as save_exp131_figures,
)
from exp131_response_constraint_validation_gate_summary import (
    write_outputs as write_exp131_outputs,
)
from exp132_response_constraint_validation_exact_sanity import (
    run_experiment as run_exp132,
)
from exp132_response_constraint_validation_exact_sanity import (
    write_outputs as write_exp132_outputs,
)
from exp133_response_handoff_exact_sanity import (
    run_experiment as run_exp133,
)
from exp133_response_handoff_exact_sanity import (
    write_outputs as write_exp133_outputs,
)
from exp134_response_handoff_threshold_sensitivity import (
    ExperimentConfig as Exp134Config,
)
from exp134_response_handoff_threshold_sensitivity import (
    run_experiment as run_exp134,
)
from exp134_response_handoff_threshold_sensitivity import (
    save_figures as save_exp134_figures,
)
from exp134_response_handoff_threshold_sensitivity import (
    write_outputs as write_exp134_outputs,
)
from exp135_response_handoff_protocol_selection import (
    ExperimentConfig as Exp135Config,
)
from exp135_response_handoff_protocol_selection import (
    run_experiment as run_exp135,
)
from exp135_response_handoff_protocol_selection import (
    save_figures as save_exp135_figures,
)
from exp135_response_handoff_protocol_selection import (
    write_outputs as write_exp135_outputs,
)
from exp136_response_handoff_manifest_export import (
    ExperimentConfig as Exp136Config,
)
from exp136_response_handoff_manifest_export import (
    run_experiment as run_exp136,
)
from exp136_response_handoff_manifest_export import (
    write_outputs as write_exp136_outputs,
)
from exp137_response_handoff_failure_catalog import (
    ExperimentConfig as Exp137Config,
)
from exp137_response_handoff_failure_catalog import (
    run_experiment as run_exp137,
)
from exp137_response_handoff_failure_catalog import (
    save_figures as save_exp137_figures,
)
from exp137_response_handoff_failure_catalog import (
    write_outputs as write_exp137_outputs,
)
from exp138_response_handoff_preregistration_rules import (
    run_experiment as run_exp138,
)
from exp138_response_handoff_preregistration_rules import (
    write_outputs as write_exp138_outputs,
)
from exp139_response_handoff_manifest_read_exact_sanity import (
    run_experiment as run_exp139,
)
from exp139_response_handoff_manifest_read_exact_sanity import (
    write_outputs as write_exp139_outputs,
)
from exp140_manifest_representability_exact_sanity import (
    run_experiment as run_exp140,
)
from exp140_manifest_representability_exact_sanity import (
    write_outputs as write_exp140_outputs,
)
from exp141_frozen_manifest_ordinal_representation import (
    ExperimentConfig as Exp141Config,
)
from exp141_frozen_manifest_ordinal_representation import (
    run_experiment as run_exp141,
)
from exp141_frozen_manifest_ordinal_representation import (
    save_figures as save_exp141_figures,
)
from exp141_frozen_manifest_ordinal_representation import (
    write_outputs as write_exp141_outputs,
)
from exp142_frozen_manifest_representation_nulls import (
    ExperimentConfig as Exp142Config,
)
from exp142_frozen_manifest_representation_nulls import (
    run_experiment as run_exp142,
)
from exp142_frozen_manifest_representation_nulls import (
    save_figures as save_exp142_figures,
)
from exp142_frozen_manifest_representation_nulls import (
    write_outputs as write_exp142_outputs,
)
from exp143_frozen_manifest_fit_stability import (
    ExperimentConfig as Exp143Config,
)
from exp143_frozen_manifest_fit_stability import (
    run_experiment as run_exp143,
)
from exp143_frozen_manifest_fit_stability import (
    save_figures as save_exp143_figures,
)
from exp143_frozen_manifest_fit_stability import (
    write_outputs as write_exp143_outputs,
)
from exp144_frozen_manifest_dimension_complexity_curve import (
    ExperimentConfig as Exp144Config,
)
from exp144_frozen_manifest_dimension_complexity_curve import (
    run_experiment as run_exp144,
)
from exp144_frozen_manifest_dimension_complexity_curve import (
    save_figures as save_exp144_figures,
)
from exp144_frozen_manifest_dimension_complexity_curve import (
    write_outputs as write_exp144_outputs,
)
from exp145_failed_manifest_no_fit_controls import (
    ExperimentConfig as Exp145Config,
)
from exp145_failed_manifest_no_fit_controls import (
    run_experiment as run_exp145,
)
from exp145_failed_manifest_no_fit_controls import (
    write_outputs as write_exp145_outputs,
)
from exp146_frozen_manifest_representation_summary import (
    run_experiment as run_exp146,
)
from exp146_frozen_manifest_representation_summary import (
    write_outputs as write_exp146_outputs,
)
from exp147_manifest_representation_no_metric_exact_sanity import (
    run_experiment as run_exp147,
)
from exp147_manifest_representation_no_metric_exact_sanity import (
    write_outputs as write_exp147_outputs,
)
from exp148_manifest_family_exact_sanity import (
    run_experiment as run_exp148,
)
from exp148_manifest_family_exact_sanity import (
    write_outputs as write_exp148_outputs,
)
from exp149_manifest_family_fit_comparison import (
    ExperimentConfig as Exp149Config,
)
from exp149_manifest_family_fit_comparison import (
    run_experiment as run_exp149,
)
from exp149_manifest_family_fit_comparison import (
    save_figures as save_exp149_figures,
)
from exp149_manifest_family_fit_comparison import (
    write_outputs as write_exp149_outputs,
)
from exp150_manifest_family_null_taxonomy import (
    ExperimentConfig as Exp150Config,
)
from exp150_manifest_family_null_taxonomy import (
    run_experiment as run_exp150,
)
from exp150_manifest_family_null_taxonomy import (
    save_figures as save_exp150_figures,
)
from exp150_manifest_family_null_taxonomy import (
    write_outputs as write_exp150_outputs,
)
from exp151_manifest_family_stricter_criteria import (
    ExperimentConfig as Exp151Config,
)
from exp151_manifest_family_stricter_criteria import (
    run_experiment as run_exp151,
)
from exp151_manifest_family_stricter_criteria import (
    save_figures as save_exp151_figures,
)
from exp151_manifest_family_stricter_criteria import (
    write_outputs as write_exp151_outputs,
)
from exp152_manifest_family_failed_manifest_accounting import (
    ExperimentConfig as Exp152Config,
)
from exp152_manifest_family_failed_manifest_accounting import (
    run_experiment as run_exp152,
)
from exp152_manifest_family_failed_manifest_accounting import (
    save_figures as save_exp152_figures,
)
from exp152_manifest_family_failed_manifest_accounting import (
    write_outputs as write_exp152_outputs,
)
from exp153_manifest_family_no_retuning_audit import (
    run_experiment as run_exp153,
)
from exp153_manifest_family_no_retuning_audit import (
    write_outputs as write_exp153_outputs,
)
from exp154_manifest_family_report_card import (
    run_experiment as run_exp154,
)
from exp154_manifest_family_report_card import (
    write_outputs as write_exp154_outputs,
)
from exp155_manifest_family_comparison_exact_sanity import (
    run_experiment as run_exp155,
)
from exp155_manifest_family_comparison_exact_sanity import (
    write_outputs as write_exp155_outputs,
)
from exp156_cross_family_robustness_exact_sanity import (
    run_experiment as run_exp156,
)
from exp156_cross_family_robustness_exact_sanity import (
    write_outputs as write_exp156_outputs,
)
from exp157_cross_family_robustness_criteria_table import (
    run_experiment as run_exp157,
)
from exp157_cross_family_robustness_criteria_table import (
    write_outputs as write_exp157_outputs,
)
from exp158_cross_family_robustness_decision import (
    ExperimentConfig as Exp158Config,
)
from exp158_cross_family_robustness_decision import (
    run_experiment as run_exp158,
)
from exp158_cross_family_robustness_decision import (
    save_figures as save_exp158_figures,
)
from exp158_cross_family_robustness_decision import (
    write_outputs as write_exp158_outputs,
)
from exp159_cross_family_robustness_threshold_sensitivity import (
    ExperimentConfig as Exp159Config,
)
from exp159_cross_family_robustness_threshold_sensitivity import (
    run_experiment as run_exp159,
)
from exp159_cross_family_robustness_threshold_sensitivity import (
    save_figures as save_exp159_figures,
)
from exp159_cross_family_robustness_threshold_sensitivity import (
    write_outputs as write_exp159_outputs,
)
from exp160_carry_forward_registry_export import (
    ExperimentConfig as Exp160Config,
)
from exp160_carry_forward_registry_export import (
    run_experiment as run_exp160,
)
from exp160_carry_forward_registry_export import (
    write_outputs as write_exp160_outputs,
)
from exp161_cross_family_failed_provisional_accounting import (
    ExperimentConfig as Exp161Config,
)
from exp161_cross_family_failed_provisional_accounting import (
    run_experiment as run_exp161,
)
from exp161_cross_family_failed_provisional_accounting import (
    save_figures as save_exp161_figures,
)
from exp161_cross_family_failed_provisional_accounting import (
    write_outputs as write_exp161_outputs,
)
from exp162_stress_test_handoff_plan import (
    ExperimentConfig as Exp162Config,
)
from exp162_stress_test_handoff_plan import (
    run_experiment as run_exp162,
)
from exp162_stress_test_handoff_plan import (
    write_outputs as write_exp162_outputs,
)
from exp163_cross_family_robustness_report_card import (
    ExperimentConfig as Exp163Config,
)
from exp163_cross_family_robustness_report_card import (
    run_experiment as run_exp163,
)
from exp163_cross_family_robustness_report_card import (
    write_outputs as write_exp163_outputs,
)
from exp164_cross_family_robustness_final_sanity import (
    run_experiment as run_exp164,
)
from exp164_cross_family_robustness_final_sanity import (
    write_outputs as write_exp164_outputs,
)
from exp165_carry_forward_failure_decomposition_exact_sanity import (
    run_experiment as run_exp165,
)
from exp165_carry_forward_failure_decomposition_exact_sanity import (
    write_outputs as write_exp165_outputs,
)
from exp166_carry_forward_failure_decomposition import (
    ExperimentConfig as Exp166Config,
)
from exp166_carry_forward_failure_decomposition import (
    run_experiment as run_exp166,
)
from exp166_carry_forward_failure_decomposition import (
    save_figures as save_exp166_figures,
)
from exp166_carry_forward_failure_decomposition import (
    write_outputs as write_exp166_outputs,
)
from exp167_cross_family_diagnostic_completeness_audit import (
    ExperimentConfig as Exp167Config,
)
from exp167_cross_family_diagnostic_completeness_audit import (
    run_experiment as run_exp167,
)
from exp167_cross_family_diagnostic_completeness_audit import (
    save_figures as save_exp167_figures,
)
from exp167_cross_family_diagnostic_completeness_audit import (
    write_outputs as write_exp167_outputs,
)
from exp168_stress_test_stop_condition_audit import (
    ExperimentConfig as Exp168Config,
)
from exp168_stress_test_stop_condition_audit import (
    run_experiment as run_exp168,
)
from exp168_stress_test_stop_condition_audit import (
    write_outputs as write_exp168_outputs,
)
from exp169_upstream_remediation_design_table import (
    ExperimentConfig as Exp169Config,
)
from exp169_upstream_remediation_design_table import (
    run_experiment as run_exp169,
)
from exp169_upstream_remediation_design_table import (
    write_outputs as write_exp169_outputs,
)
from exp170_missing_metric_impact_report import (
    ExperimentConfig as Exp170Config,
)
from exp170_missing_metric_impact_report import (
    run_experiment as run_exp170,
)
from exp170_missing_metric_impact_report import (
    write_outputs as write_exp170_outputs,
)
from exp171_failure_decomposition_no_retuning_audit import (
    run_experiment as run_exp171,
)
from exp171_failure_decomposition_no_retuning_audit import (
    write_outputs as write_exp171_outputs,
)
from exp172_carry_forward_failure_report_card import (
    ExperimentConfig as Exp172Config,
)
from exp172_carry_forward_failure_report_card import (
    run_experiment as run_exp172,
)
from exp172_carry_forward_failure_report_card import (
    write_outputs as write_exp172_outputs,
)
from exp173_carry_forward_failure_final_sanity import (
    run_experiment as run_exp173,
)
from exp173_carry_forward_failure_final_sanity import (
    write_outputs as write_exp173_outputs,
)
from exp174_remediation_plan_exact_sanity import (
    run_experiment as run_exp174,
)
from exp174_remediation_plan_exact_sanity import (
    write_outputs as write_exp174_outputs,
)
from exp175_diagnostic_complete_schema_export import (
    run_experiment as run_exp175,
)
from exp175_diagnostic_complete_schema_export import (
    write_outputs as write_exp175_outputs,
)
from exp176_failure_to_remediation_mapping import (
    ExperimentConfig as Exp176Config,
)
from exp176_failure_to_remediation_mapping import (
    run_experiment as run_exp176,
)
from exp176_failure_to_remediation_mapping import (
    write_outputs as write_exp176_outputs,
)
from exp177_new_manifest_family_design_v2 import (
    run_experiment as run_exp177,
)
from exp177_new_manifest_family_design_v2 import (
    write_outputs as write_exp177_outputs,
)
from exp178_preregistered_remediation_plan_export import (
    ExperimentConfig as Exp178Config,
)
from exp178_preregistered_remediation_plan_export import (
    run_experiment as run_exp178,
)
from exp178_preregistered_remediation_plan_export import (
    write_outputs as write_exp178_outputs,
)
from exp179_future_manifest_run_spec import (
    ExperimentConfig as Exp179Config,
)
from exp179_future_manifest_run_spec import (
    run_experiment as run_exp179,
)
from exp179_future_manifest_run_spec import (
    write_outputs as write_exp179_outputs,
)
from exp180_remediation_no_execution_audit import (
    ExperimentConfig as Exp180Config,
)
from exp180_remediation_no_execution_audit import (
    run_experiment as run_exp180,
)
from exp180_remediation_no_execution_audit import (
    write_outputs as write_exp180_outputs,
)
from exp181_remediation_plan_report_card import (
    ExperimentConfig as Exp181Config,
)
from exp181_remediation_plan_report_card import (
    run_experiment as run_exp181,
)
from exp181_remediation_plan_report_card import (
    write_outputs as write_exp181_outputs,
)
from exp182_remediation_plan_final_sanity import (
    run_experiment as run_exp182,
)
from exp182_remediation_plan_final_sanity import (
    write_outputs as write_exp182_outputs,
)
from exp183_v2_manifest_spec_exact_sanity import (
    run_experiment as run_exp183,
)
from exp183_v2_manifest_spec_exact_sanity import (
    write_outputs as write_exp183_outputs,
)
from exp184_v2_manifest_generation import (
    ExperimentConfig as Exp184Config,
)
from exp184_v2_manifest_generation import (
    run_experiment as run_exp184,
)
from exp184_v2_manifest_generation import (
    write_outputs as write_exp184_outputs,
)
from exp185_v2_family_fit_diagnostics import (
    ExperimentConfig as Exp185Config,
)
from exp185_v2_family_fit_diagnostics import (
    run_experiment as run_exp185,
)
from exp185_v2_family_fit_diagnostics import (
    save_figures as save_exp185_figures,
)
from exp185_v2_family_fit_diagnostics import (
    write_outputs as write_exp185_outputs,
)
from exp186_v2_null_taxonomy_diagnostics import (
    ExperimentConfig as Exp186Config,
)
from exp186_v2_null_taxonomy_diagnostics import (
    run_experiment as run_exp186,
)
from exp186_v2_null_taxonomy_diagnostics import (
    save_figures as save_exp186_figures,
)
from exp186_v2_null_taxonomy_diagnostics import (
    write_outputs as write_exp186_outputs,
)
from exp187_v2_stricter_criteria_diagnostics import (
    ExperimentConfig as Exp187Config,
)
from exp187_v2_stricter_criteria_diagnostics import (
    run_experiment as run_exp187,
)
from exp187_v2_stricter_criteria_diagnostics import (
    save_figures as save_exp187_figures,
)
from exp187_v2_stricter_criteria_diagnostics import (
    write_outputs as write_exp187_outputs,
)
from exp188_v2_failed_accounting import (
    ExperimentConfig as Exp188Config,
)
from exp188_v2_failed_accounting import (
    run_experiment as run_exp188,
)
from exp188_v2_failed_accounting import (
    save_figures as save_exp188_figures,
)
from exp188_v2_failed_accounting import (
    write_outputs as write_exp188_outputs,
)
from exp189_v2_coverage_metrics import (
    ExperimentConfig as Exp189Config,
)
from exp189_v2_coverage_metrics import (
    run_experiment as run_exp189,
)
from exp189_v2_coverage_metrics import (
    write_outputs as write_exp189_outputs,
)
from exp190_v2_restart_latent_order_stability import (
    ExperimentConfig as Exp190Config,
)
from exp190_v2_restart_latent_order_stability import (
    run_experiment as run_exp190,
)
from exp190_v2_restart_latent_order_stability import (
    save_figures as save_exp190_figures,
)
from exp190_v2_restart_latent_order_stability import (
    write_outputs as write_exp190_outputs,
)
from exp191_v2_no_retuning_audit import (
    ExperimentConfig as Exp191Config,
)
from exp191_v2_no_retuning_audit import (
    run_experiment as run_exp191,
)
from exp191_v2_no_retuning_audit import (
    write_outputs as write_exp191_outputs,
)
from exp192_v2_required_metric_aggregation import (
    ExperimentConfig as Exp192Config,
)
from exp192_v2_required_metric_aggregation import (
    run_experiment as run_exp192,
)
from exp192_v2_required_metric_aggregation import (
    save_figures as save_exp192_figures,
)
from exp192_v2_required_metric_aggregation import (
    write_outputs as write_exp192_outputs,
)
from exp193_v2_diagnostic_complete_bundle_report import (
    ExperimentConfig as Exp193Config,
)
from exp193_v2_diagnostic_complete_bundle_report import (
    run_experiment as run_exp193,
)
from exp193_v2_diagnostic_complete_bundle_report import (
    write_outputs as write_exp193_outputs,
)
from exp194_v2_manifest_generation_final_sanity import (
    run_experiment as run_exp194,
)
from exp194_v2_manifest_generation_final_sanity import (
    write_outputs as write_exp194_outputs,
)
from exp195_v2_carry_forward_exact_sanity import (
    run_experiment as run_exp195,
)
from exp195_v2_carry_forward_exact_sanity import (
    write_outputs as write_exp195_outputs,
)
from exp196_v2_bundle_input_audit import (
    ExperimentConfig as Exp196Config,
)
from exp196_v2_bundle_input_audit import (
    run_experiment as run_exp196,
)
from exp196_v2_bundle_input_audit import (
    write_outputs as write_exp196_outputs,
)
from exp197_v2_carry_forward_decision import (
    ExperimentConfig as Exp197Config,
)
from exp197_v2_carry_forward_decision import (
    run_experiment as run_exp197,
)
from exp197_v2_carry_forward_decision import (
    save_figures as save_exp197_figures,
)
from exp197_v2_carry_forward_decision import (
    write_outputs as write_exp197_outputs,
)
from exp198_v2_carry_forward_threshold_sensitivity import (
    ExperimentConfig as Exp198Config,
)
from exp198_v2_carry_forward_threshold_sensitivity import (
    run_experiment as run_exp198,
)
from exp198_v2_carry_forward_threshold_sensitivity import (
    save_figures as save_exp198_figures,
)
from exp198_v2_carry_forward_threshold_sensitivity import (
    write_outputs as write_exp198_outputs,
)
from exp199_v2_carry_forward_registry_export import (
    ExperimentConfig as Exp199Config,
)
from exp199_v2_carry_forward_registry_export import (
    run_experiment as run_exp199,
)
from exp199_v2_carry_forward_registry_export import (
    write_outputs as write_exp199_outputs,
)
from exp200_v2_stress_test_handoff_plan import (
    ExperimentConfig as Exp200Config,
)
from exp200_v2_stress_test_handoff_plan import (
    run_experiment as run_exp200,
)
from exp200_v2_stress_test_handoff_plan import (
    write_outputs as write_exp200_outputs,
)
from exp201_v2_failed_provisional_accounting import (
    ExperimentConfig as Exp201Config,
)
from exp201_v2_failed_provisional_accounting import (
    run_experiment as run_exp201,
)
from exp201_v2_failed_provisional_accounting import (
    save_figures as save_exp201_figures,
)
from exp201_v2_failed_provisional_accounting import (
    write_outputs as write_exp201_outputs,
)
from exp202_v2_carry_forward_failure_decomposition import (
    ExperimentConfig as Exp202Config,
)
from exp202_v2_carry_forward_failure_decomposition import (
    run_experiment as run_exp202,
)
from exp202_v2_carry_forward_failure_decomposition import (
    save_figures as save_exp202_figures,
)
from exp202_v2_carry_forward_failure_decomposition import (
    write_outputs as write_exp202_outputs,
)
from exp203_v2_carry_forward_no_retuning_audit import (
    ExperimentConfig as Exp203Config,
)
from exp203_v2_carry_forward_no_retuning_audit import (
    run_experiment as run_exp203,
)
from exp203_v2_carry_forward_no_retuning_audit import (
    write_outputs as write_exp203_outputs,
)
from exp204_v2_carry_forward_report_card import (
    ExperimentConfig as Exp204Config,
)
from exp204_v2_carry_forward_report_card import (
    run_experiment as run_exp204,
)
from exp204_v2_carry_forward_report_card import (
    write_outputs as write_exp204_outputs,
)
from exp205_v2_carry_forward_final_sanity import (
    run_experiment as run_exp205,
)
from exp205_v2_carry_forward_final_sanity import (
    write_outputs as write_exp205_outputs,
)
from exp206_v2_blocking_analysis_exact_sanity import (
    main as run_exp206_main,
)
from exp207_v2_blocked_root_cause_audit import (
    main as run_exp207_main,
)
from exp208_v2_criterion_margin_report import (
    main as run_exp208_main,
)
from exp209_v2_structural_vs_measured_blocking import (
    main as run_exp209_main,
)
from exp210_v3_manifest_family_design import (
    main as run_exp210_main,
)
from exp211_v3_preregistration_export import (
    main as run_exp211_main,
)
from exp212_v3_no_execution_audit import (
    main as run_exp212_main,
)
from exp213_v2_blocked_decision_report_card import (
    main as run_exp213_main,
)
from exp214_v2_blocked_v3_preregistration_final_sanity import (
    main as run_exp214_main,
)
from exp215_protocol_metadata_exact_sanity import (
    main as run_exp215_main,
)
from exp216_response_profile_invariance_exact_sanity import (
    main as run_exp216_main,
)
from exp217_current_profile_protocol_invariance_audit import (
    main as run_exp217_main,
)
from exp218_v3_protocol_invariant_patch_design import (
    main as run_exp218_main,
)
from exp219_v3_protocol_patched_preregistration_export import (
    main as run_exp219_main,
)
from exp220_v3_protocol_patch_audit import (
    main as run_exp220_main,
)
from exp221_protocol_invariance_language_audit import (
    main as run_exp221_main,
)
from exp222_protocol_patch_no_execution_audit import (
    main as run_exp222_main,
)
from exp223_protocol_invariance_final_sanity import (
    main as run_exp223_main,
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


def run_echo_response_signature_exact_sanity() -> None:
    rows = run_exp100()
    output_path = write_exp100_outputs(rows)
    print(f"exp100 wrote {output_path}")


def run_layered_echo_response_order_recovery() -> None:
    config = Exp101Config(
        reference_length=48,
        emission_position=6,
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer_values=(5,),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp101(config)
    data_path = write_exp101_outputs(rows, config.output_dir)
    figure_paths = save_exp101_figures(rows, config.output_dir)
    print(f"exp101 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_response_signature_coarse_protocol_stability() -> None:
    config = Exp102Config(
        reference_length=48,
        emission_position=6,
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        event_keep_probabilities=(1.0, 0.5),
        reference_strides=(1, 2),
        edge_removal_probabilities=(0.0, 0.1),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp102(config)
    data_path = write_exp102_outputs(rows, config.output_dir)
    figure_paths = save_exp102_figures(rows, config.output_dir)
    print(f"exp102 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_response_shortcut_robust_core() -> None:
    config = Exp103Config(
        reference_length=48,
        emission_position=6,
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        shortcut_probabilities=(0.0, 0.3),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp103(config)
    data_path = write_exp103_outputs(rows, config.output_dir)
    figure_paths = save_exp103_figures(rows, config.output_dir)
    print(f"exp103 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_response_reference_protocol_dependence() -> None:
    config = Exp104Config(
        num_systems=5,
        max_events=150,
        trigger_probability=0.20,
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        repetitions=2,
        seed=0,
        top_k=3,
        output_dir=Path("outputs"),
    )
    rows = run_exp104(config)
    data_path = write_exp104_outputs(rows, config.output_dir)
    figure_paths = save_exp104_figures(rows, config.output_dir)
    print(f"exp104 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_response_order_precondition_diagnostics() -> None:
    config = Exp105Config(
        reference_length=48,
        emission_position=6,
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        shortcut_probabilities=(0.0, 0.3),
        reference_strides=(1, 2),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp105(config)
    data_path = write_exp105_outputs(rows, config.output_dir)
    figure_paths = save_exp105_figures(rows, config.output_dir)
    print(f"exp105 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_response_signature_stability_exact_sanity() -> None:
    rows = run_exp106()
    output_path = write_exp106_outputs(rows)
    print(f"exp106 wrote {output_path}")


def run_echo_order_semantics_exact_sanity() -> None:
    rows = run_exp107()
    output_path = write_exp107_outputs(rows)
    print(f"exp107 wrote {output_path}")


def run_gated_echo_protocol_comparison() -> None:
    config = Exp108Config(
        reference_length=32,
        motif_count=20,
        delay_ranks=(5, 8, 13),
        shortcut_probability_values=(0.0, 0.5),
        gate_delay_ranks=(1, 5),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp108(config)
    data_path = write_exp108_outputs(rows, config.output_dir)
    figure_paths = save_exp108_figures(rows, config.output_dir)
    print(f"exp108 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_delay_spacing_tie_resolution() -> None:
    config = Exp109Config(
        reference_length_values=(64,),
        emission_position=8,
        layer_sets=("compact", "medium"),
        targets_per_layer=5,
        reference_strides=(1, 4),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp109(config)
    data_path = write_exp109_outputs(rows, config.output_dir)
    figure_paths = save_exp109_figures(rows, config.output_dir)
    print(f"exp109 wrote {data_path} and {len(figure_paths)} figures")


def run_response_order_scalar_representability() -> None:
    config = Exp110Config(
        reference_length=48,
        emission_position=6,
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        shortcut_probabilities=(0.0, 0.5),
        reference_strides=(1, 2),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp110(config)
    data_path = write_exp110_outputs(rows, config.output_dir)
    figure_paths = save_exp110_figures(rows, config.output_dir)
    print(f"exp110 wrote {data_path} and {len(figure_paths)} figures")


def run_echo_terminology_audit() -> None:
    rows = run_exp111()
    output_path = write_exp111_outputs(rows)
    print(f"exp111 wrote {output_path}")


def run_response_representability_exact_sanity() -> None:
    rows = run_exp112()
    output_path = write_exp112_outputs(rows)
    print(f"exp112 wrote {output_path}")


def run_response_order_underdetermination_exact_sanity() -> None:
    rows = run_exp113()
    output_path = write_exp113_outputs(rows)
    print(f"exp113 wrote {output_path}")


def run_single_reference_response_order_underdetermination() -> None:
    config = Exp114Config(
        target_counts=(8, 16),
        unique_rank_counts=(3, 5),
        layout_count=20,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp114(config)
    data_path = write_exp114_outputs(rows, config.output_dir)
    figure_paths = save_exp114_figures(rows, config.output_dir)
    print(f"exp114 wrote {data_path} and {len(figure_paths)} figures")


def run_multi_reference_response_profile_diagnostics() -> None:
    config = Exp115Config(
        reference_length=48,
        emission_positions=(6, 12),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp115(config)
    data_path = write_exp115_outputs(rows, config.output_dir)
    figure_paths = save_exp115_figures(rows, config.output_dir)
    print(f"exp115 wrote {data_path} and {len(figure_paths)} figures")


def run_response_representability_requirement_table() -> None:
    rows = run_exp116()
    output_path = write_exp116_outputs(rows)
    print(f"exp116 wrote {output_path}")


def run_response_profile_stability_under_protocol_variants() -> None:
    config = Exp117Config(
        reference_length=48,
        emission_positions=(6, 12),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        shortcut_probability_values=(0.0, 0.3),
        reference_strides=(1, 2),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp117(config)
    data_path = write_exp117_outputs(rows, config.output_dir)
    figure_paths = save_exp117_figures(rows, config.output_dir)
    print(f"exp117 wrote {data_path} and {len(figure_paths)} figures")


def run_response_profile_requirement_exact_sanity() -> None:
    rows = run_exp118()
    output_path = write_exp118_outputs(rows)
    print(f"exp118 wrote {output_path}")


def run_pairwise_response_profile_exact_sanity() -> None:
    rows = run_exp119()
    output_path = write_exp119_outputs(rows)
    print(f"exp119 wrote {output_path}")


def run_pairwise_response_protocol_comparison() -> None:
    config = Exp120Config(
        reference_length=48,
        emission_positions=(6, 12),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp120(config)
    data_path = write_exp120_outputs(rows, config.output_dir)
    figure_paths = save_exp120_figures(rows, config.output_dir)
    print(f"exp120 wrote {data_path} and {len(figure_paths)} figures")


def run_pairwise_response_null_baselines() -> None:
    config = Exp121Config(
        reference_length=48,
        emission_positions=(6, 12),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        repetitions=2,
        null_repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp121(config)
    data_path = write_exp121_outputs(rows, config.output_dir)
    figure_paths = save_exp121_figures(rows, config.output_dir)
    print(f"exp121 wrote {data_path} and {len(figure_paths)} figures")


def run_pairwise_response_protocol_variant_stability() -> None:
    config = Exp122Config(
        reference_length=48,
        emission_positions=(6, 12),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        shortcut_probability_values=(0.0, 0.3),
        reference_strides=(1, 2),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp122(config)
    data_path = write_exp122_outputs(rows, config.output_dir)
    figure_paths = save_exp122_figures(rows, config.output_dir)
    print(f"exp122 wrote {data_path} and {len(figure_paths)} figures")


def run_pairwise_response_missing_data_sensitivity() -> None:
    config = Exp123Config(
        target_count=30,
        protocol_count=4,
        reachable_probabilities=(0.5, 1.0),
        unique_rank_count=5,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp123(config)
    data_path = write_exp123_outputs(rows, config.output_dir)
    figure_paths = save_exp123_figures(rows, config.output_dir)
    print(f"exp123 wrote {data_path} and {len(figure_paths)} figures")


def run_pairwise_response_protocol_choice_dependence() -> None:
    config = Exp124Config(
        target_count=30,
        protocol_count=4,
        reachable_probability=0.8,
        unique_rank_count=5,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp124(config)
    data_path = write_exp124_outputs(rows, config.output_dir)
    figure_paths = save_exp124_figures(rows, config.output_dir)
    print(f"exp124 wrote {data_path} and {len(figure_paths)} figures")


def run_pairwise_response_null_admissibility_exact_sanity() -> None:
    rows = run_exp125()
    output_path = write_exp125_outputs(rows)
    print(f"exp125 wrote {output_path}")


def run_response_constraint_pool_exact_sanity() -> None:
    rows = run_exp126()
    output_path = write_exp126_outputs(rows)
    print(f"exp126 wrote {output_path}")


def run_response_constraint_heldout_protocol_validation() -> None:
    config = Exp127Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        train_fractions=(0.5,),
        max_constraints=1000,
        min_margins=(0.0, 0.05),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp127(config)
    data_path = write_exp127_outputs(rows, config.output_dir)
    figure_paths = save_exp127_figures(rows, config.output_dir)
    print(f"exp127 wrote {data_path} and {len(figure_paths)} figures")


def run_response_constraint_bootstrap_stability() -> None:
    config = Exp128Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        max_constraints=1000,
        min_margin=0.05,
        bootstrap_count=10,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp128(config)
    data_path = write_exp128_outputs(rows, config.output_dir)
    figure_paths = save_exp128_figures(rows, config.output_dir)
    print(f"exp128 wrote {data_path} and {len(figure_paths)} figures")


def run_response_constraint_null_separation() -> None:
    config = Exp129Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        max_constraints=1000,
        min_margin=0.05,
        null_repetitions=3,
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    summary_rows, by_type_rows = run_exp129(config)
    data_paths = write_exp129_outputs(summary_rows, by_type_rows, config.output_dir)
    figure_paths = save_exp129_figures(summary_rows, by_type_rows, config.output_dir)
    print(
        f"exp129 wrote {data_paths[0]}, {data_paths[1]}, "
        f"and {len(figure_paths)} figures"
    )


def run_response_constraint_pool_coverage() -> None:
    config = Exp130Config(
        target_counts=(20,),
        protocol_counts=(3, 5),
        unique_rank_count=5,
        reachable_probability=0.8,
        max_constraints_values=(500, 1000),
        min_margins=(0.0, 0.05),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp130(config)
    data_path = write_exp130_outputs(rows, config.output_dir)
    figure_paths = save_exp130_figures(rows, config.output_dir)
    print(f"exp130 wrote {data_path} and {len(figure_paths)} figures")


def run_response_constraint_validation_gate_summary() -> None:
    config = Exp131Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        max_constraints=1000,
        min_margin_values=(0.0, 0.05),
        repetitions=2,
        bootstrap_count=10,
        null_repetitions=3,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp131(config)
    data_path = write_exp131_outputs(rows, config.output_dir)
    figure_paths = save_exp131_figures(rows, config.output_dir)
    print(f"exp131 wrote {data_path} and {len(figure_paths)} figures")


def run_response_constraint_validation_exact_sanity() -> None:
    rows = run_exp132()
    output_path = write_exp132_outputs(rows)
    print(f"exp132 wrote {output_path}")


def run_response_handoff_exact_sanity() -> None:
    rows = run_exp133()
    output_path = write_exp133_outputs(rows)
    print(f"exp133 wrote {output_path}")


def run_response_handoff_threshold_sensitivity() -> None:
    config = Exp134Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        min_agreement_values=(0.7, 0.8),
        min_null_z_values=(0.0, 1.0),
        min_bootstrap_values=(0.6, 0.7),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp134(config)
    data_path = write_exp134_outputs(rows, config.output_dir)
    figure_paths = save_exp134_figures(rows, config.output_dir)
    print(f"exp134 wrote {data_path} and {len(figure_paths)} figures")


def run_response_handoff_protocol_selection() -> None:
    config = Exp135Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        min_margins=(0.0, 0.05),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp135(config)
    data_path = write_exp135_outputs(rows, config.output_dir)
    figure_paths = save_exp135_figures(rows, config.output_dir)
    print(f"exp135 wrote {data_path} and {len(figure_paths)} figures")


def run_response_handoff_manifest_export() -> None:
    config = Exp136Config(
        reference_length=48,
        emission_positions=(6, 12, 18),
        layer_delay_ranks=(3, 5, 8),
        targets_per_layer=5,
        max_manifests=3,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows, paths = run_exp136(config)
    data_path = write_exp136_outputs(rows, config.output_dir)
    print(f"exp136 wrote {data_path} and {len(paths)} manifests")


def run_response_handoff_failure_catalog() -> None:
    config = Exp137Config(
        target_count=30,
        protocol_count_values=(2, 4),
        reachable_probabilities=(0.3, 0.8),
        unique_rank_counts=(3, 5),
        repetitions=2,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp137(config)
    data_path = write_exp137_outputs(rows, config.output_dir)
    figure_paths = save_exp137_figures(rows, config.output_dir)
    print(f"exp137 wrote {data_path} and {len(figure_paths)} figures")


def run_response_handoff_preregistration_rules() -> None:
    rows = run_exp138()
    output_path = write_exp138_outputs(rows)
    print(f"exp138 wrote {output_path}")


def run_response_handoff_manifest_read_exact_sanity() -> None:
    rows = run_exp139()
    output_path = write_exp139_outputs(rows)
    print(f"exp139 wrote {output_path}")


def run_manifest_representability_exact_sanity() -> None:
    rows = run_exp140()
    output_path = write_exp140_outputs(rows)
    print(f"exp140 wrote {output_path}")


def run_frozen_manifest_ordinal_representation() -> None:
    config = Exp141Config(
        manifest_dir=Path("outputs/manifests"),
        dims=(1, 2),
        steps=300,
        restarts=1,
        learning_rate=0.05,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp141(config)
    data_path = write_exp141_outputs(rows, config.output_dir)
    figure_paths = save_exp141_figures(rows, config.output_dir)
    print(f"exp141 wrote {data_path} and {len(figure_paths)} figures")


def run_frozen_manifest_representation_nulls() -> None:
    config = Exp142Config(
        manifest_dir=Path("outputs/manifests"),
        embedding_dim=2,
        null_repetitions=2,
        steps=300,
        restarts=1,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp142(config)
    data_path = write_exp142_outputs(rows, config.output_dir)
    figure_paths = save_exp142_figures(rows, config.output_dir)
    print(f"exp142 wrote {data_path} and {len(figure_paths)} figures")


def run_frozen_manifest_fit_stability() -> None:
    config = Exp143Config(
        manifest_dir=Path("outputs/manifests"),
        embedding_dim=2,
        restart_count=3,
        steps=300,
        learning_rate=0.05,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp143(config)
    data_path = write_exp143_outputs(rows, config.output_dir)
    figure_paths = save_exp143_figures(rows, config.output_dir)
    print(f"exp143 wrote {data_path} and {len(figure_paths)} figures")


def run_frozen_manifest_dimension_complexity_curve() -> None:
    config = Exp144Config(
        manifest_dir=Path("outputs/manifests"),
        candidate_dims=(1, 2, 3),
        steps=300,
        restarts=1,
        learning_rate=0.05,
        complexity_lambda=0.01,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp144(config)
    data_path = write_exp144_outputs(rows, config.output_dir)
    figure_paths = save_exp144_figures(rows, config.output_dir)
    print(f"exp144 wrote {data_path} and {len(figure_paths)} figures")


def run_failed_manifest_no_fit_controls() -> None:
    config = Exp145Config(
        manifest_dir=Path("outputs/manifests"),
        generate_failed_controls=True,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp145(config)
    data_path = write_exp145_outputs(rows, config.output_dir)
    print(f"exp145 wrote {data_path}")


def run_frozen_manifest_representation_summary() -> None:
    rows = run_exp146()
    output_path = write_exp146_outputs(rows)
    print(f"exp146 wrote {output_path}")


def run_manifest_representation_no_metric_exact_sanity() -> None:
    rows = run_exp147()
    output_path = write_exp147_outputs(rows)
    print(f"exp147 wrote {output_path}")


def run_manifest_family_exact_sanity() -> None:
    rows = run_exp148()
    output_path = write_exp148_outputs(rows)
    print(f"exp148 wrote {output_path}")


def run_manifest_family_fit_comparison() -> None:
    config = Exp149Config(
        manifest_dir=Path("outputs/manifests"),
        dims=(1, 2),
        steps=300,
        restarts=1,
        seed=0,
        output_dir=Path("outputs"),
    )
    fit_rows, summary_rows, failed_rows = run_exp149(config)
    paths = write_exp149_outputs(
        fit_rows,
        summary_rows,
        failed_rows,
        config.output_dir,
    )
    figure_paths = save_exp149_figures(summary_rows, config.output_dir)
    print(f"exp149 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp149 wrote {len(figure_paths)} figures")


def run_manifest_family_null_taxonomy() -> None:
    config = Exp150Config(
        manifest_dir=Path("outputs/manifests"),
        embedding_dim=2,
        null_repetitions=2,
        steps=300,
        restarts=1,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows, taxonomy_rows = run_exp150(config)
    data_path = write_exp150_outputs(rows, config.output_dir)
    figure_paths = save_exp150_figures(taxonomy_rows, config.output_dir)
    print(f"exp150 wrote {data_path} and {len(figure_paths)} figures")


def run_manifest_family_stricter_criteria() -> None:
    config = Exp151Config(
        manifest_dir=Path("outputs/manifests"),
        dims=(1, 2),
        steps=300,
        restarts=1,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp151(config)
    data_path = write_exp151_outputs(rows, config.output_dir)
    figure_paths = save_exp151_figures(rows, config.output_dir)
    print(f"exp151 wrote {data_path} and {len(figure_paths)} figures")


def run_manifest_family_failed_manifest_accounting() -> None:
    config = Exp152Config(
        manifest_dir=Path("outputs/manifests"),
        generate_failed_controls=True,
        seed=0,
        output_dir=Path("outputs"),
    )
    rows = run_exp152(config)
    data_path = write_exp152_outputs(rows, config.output_dir)
    figure_paths = save_exp152_figures(rows, config.output_dir)
    print(f"exp152 wrote {data_path} and {len(figure_paths)} figures")


def run_manifest_family_no_retuning_audit() -> None:
    rows = run_exp153()
    output_path = write_exp153_outputs(rows)
    print(f"exp153 wrote {output_path}")


def run_manifest_family_report_card() -> None:
    rows = run_exp154()
    output_path = write_exp154_outputs(rows)
    print(f"exp154 wrote {output_path}")


def run_manifest_family_comparison_exact_sanity() -> None:
    rows = run_exp155()
    output_path = write_exp155_outputs(rows)
    print(f"exp155 wrote {output_path}")


def run_cross_family_robustness_exact_sanity() -> None:
    rows = run_exp156()
    output_path = write_exp156_outputs(rows)
    print(f"exp156 wrote {output_path}")


def run_cross_family_robustness_criteria_table() -> None:
    rows = run_exp157()
    output_path = write_exp157_outputs(rows)
    print(f"exp157 wrote {output_path}")


def run_cross_family_robustness_decision() -> None:
    config = Exp158Config(output_dir=Path("outputs"))
    metrics, decisions = run_exp158(config)
    paths = write_exp158_outputs(metrics, decisions, config.output_dir)
    figure_paths = save_exp158_figures(decisions, config.output_dir)
    print(f"exp158 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp158 wrote {len(figure_paths)} figures")


def run_cross_family_robustness_threshold_sensitivity() -> None:
    config = Exp159Config(
        output_dir=Path("outputs"),
        heldout_thresholds=(0.20, 0.25),
        null_gap_thresholds=(0.05, 0.10),
        stricter_pass_thresholds=(0.25, 0.50),
    )
    rows = run_exp159(config)
    output_path = write_exp159_outputs(rows, config.output_dir)
    figure_paths = save_exp159_figures(rows, config.output_dir)
    print(f"exp159 wrote {output_path} and {len(figure_paths)} figures")


def run_carry_forward_registry_export() -> None:
    config = Exp160Config(output_dir=Path("outputs"))
    registry_path, rows = run_exp160(config)
    output_path = write_exp160_outputs(rows, config.output_dir)
    print(f"exp160 wrote {output_path} and {registry_path}")


def run_cross_family_failed_provisional_accounting() -> None:
    config = Exp161Config(output_dir=Path("outputs"))
    rows = run_exp161(config)
    output_path = write_exp161_outputs(rows, config.output_dir)
    figure_paths = save_exp161_figures(rows, config.output_dir)
    print(f"exp161 wrote {output_path} and {len(figure_paths)} figures")


def run_stress_test_handoff_plan() -> None:
    config = Exp162Config(output_dir=Path("outputs"))
    rows = run_exp162(config)
    output_path = write_exp162_outputs(rows, config.output_dir)
    print(f"exp162 wrote {output_path}")


def run_cross_family_robustness_report_card() -> None:
    config = Exp163Config(output_dir=Path("outputs"))
    rows = run_exp163(config)
    output_path = write_exp163_outputs(rows, config.output_dir)
    print(f"exp163 wrote {output_path}")


def run_cross_family_robustness_final_sanity() -> None:
    rows = run_exp164()
    output_path = write_exp164_outputs(rows)
    print(f"exp164 wrote {output_path}")


def run_carry_forward_failure_decomposition_exact_sanity() -> None:
    rows = run_exp165()
    output_path = write_exp165_outputs(rows)
    print(f"exp165 wrote {output_path}")


def run_carry_forward_failure_decomposition() -> None:
    config = Exp166Config(output_dir=Path("outputs"))
    records, summary_rows = run_exp166(config)
    paths = write_exp166_outputs(records, summary_rows, config.output_dir)
    figure_paths = save_exp166_figures(summary_rows, config.output_dir)
    print(f"exp166 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp166 wrote {len(figure_paths)} figures")


def run_cross_family_diagnostic_completeness_audit() -> None:
    config = Exp167Config(output_dir=Path("outputs"))
    rows = run_exp167(config)
    output_path = write_exp167_outputs(rows, config.output_dir)
    figure_paths = save_exp167_figures(rows, config.output_dir)
    print(f"exp167 wrote {output_path} and {len(figure_paths)} figures")


def run_stress_test_stop_condition_audit() -> None:
    config = Exp168Config(output_dir=Path("outputs"))
    rows = run_exp168(config)
    output_path = write_exp168_outputs(rows, config.output_dir)
    print(f"exp168 wrote {output_path}")


def run_upstream_remediation_design_table() -> None:
    config = Exp169Config(output_dir=Path("outputs"))
    rows = run_exp169(config)
    output_path = write_exp169_outputs(rows, config.output_dir)
    print(f"exp169 wrote {output_path}")


def run_missing_metric_impact_report() -> None:
    config = Exp170Config(output_dir=Path("outputs"))
    rows = run_exp170(config)
    output_path = write_exp170_outputs(rows, config.output_dir)
    print(f"exp170 wrote {output_path}")


def run_failure_decomposition_no_retuning_audit() -> None:
    rows = run_exp171()
    output_path = write_exp171_outputs(rows)
    print(f"exp171 wrote {output_path}")


def run_carry_forward_failure_report_card() -> None:
    config = Exp172Config(output_dir=Path("outputs"))
    rows = run_exp172(config)
    output_path = write_exp172_outputs(rows, config.output_dir)
    print(f"exp172 wrote {output_path}")


def run_carry_forward_failure_final_sanity() -> None:
    rows = run_exp173()
    output_path = write_exp173_outputs(rows)
    print(f"exp173 wrote {output_path}")


def run_remediation_plan_exact_sanity() -> None:
    rows = run_exp174()
    output_path = write_exp174_outputs(rows)
    print(f"exp174 wrote {output_path}")


def run_diagnostic_complete_schema_export() -> None:
    rows, summary = run_exp175()
    paths = write_exp175_outputs(rows, summary)
    print(f"exp175 wrote {', '.join(str(path) for path in paths)}")


def run_failure_to_remediation_mapping() -> None:
    config = Exp176Config(output_dir=Path("outputs"))
    rows = run_exp176(config)
    output_path = write_exp176_outputs(rows, config.output_dir)
    print(f"exp176 wrote {output_path}")


def run_new_manifest_family_design_v2() -> None:
    rows = run_exp177()
    output_path = write_exp177_outputs(rows)
    print(f"exp177 wrote {output_path}")


def run_preregistered_remediation_plan_export() -> None:
    config = Exp178Config(output_dir=Path("outputs"))
    plan_path, rows = run_exp178(config)
    output_path = write_exp178_outputs(rows, config.output_dir)
    print(f"exp178 wrote {output_path} and {plan_path}")


def run_future_manifest_run_spec() -> None:
    config = Exp179Config(output_dir=Path("outputs"))
    spec_path, rows = run_exp179(config)
    output_path = write_exp179_outputs(rows, config.output_dir)
    print(f"exp179 wrote {output_path} and {spec_path}")


def run_remediation_no_execution_audit() -> None:
    config = Exp180Config(output_dir=Path("outputs"))
    rows = run_exp180(config)
    output_path = write_exp180_outputs(rows, config.output_dir)
    print(f"exp180 wrote {output_path}")


def run_remediation_plan_report_card() -> None:
    config = Exp181Config(output_dir=Path("outputs"))
    rows = run_exp181(config)
    output_path = write_exp181_outputs(rows, config.output_dir)
    print(f"exp181 wrote {output_path}")


def run_remediation_plan_final_sanity() -> None:
    rows = run_exp182()
    output_path = write_exp182_outputs(rows)
    print(f"exp182 wrote {output_path}")


def run_v2_manifest_spec_exact_sanity() -> None:
    rows = run_exp183()
    output_path = write_exp183_outputs(rows)
    print(f"exp183 wrote {output_path}")


def run_v2_manifest_generation() -> None:
    config = Exp184Config(
        output_dir=Path("outputs"),
        seed=0,
        max_constraints=1000,
        min_margin=0.05,
        bootstrap_count=5,
        null_repetitions=3,
    )
    paths, rows = run_exp184(config)
    output_path = write_exp184_outputs(rows, config.output_dir)
    print(f"exp184 wrote {output_path} and {len(paths)} manifests")


def run_v2_family_fit_diagnostics() -> None:
    config = Exp185Config(
        manifest_dir=Path("outputs/manifests_v2"),
        dims=(1, 2),
        steps=300,
        restarts=1,
        output_dir=Path("outputs"),
    )
    rows, summary_rows = run_exp185(config)
    paths = write_exp185_outputs(rows, summary_rows, config.output_dir)
    figure_paths = save_exp185_figures(summary_rows, config.output_dir)
    print(f"exp185 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp185 wrote {len(figure_paths)} figures")


def run_v2_null_taxonomy_diagnostics() -> None:
    config = Exp186Config(
        manifest_dir=Path("outputs/manifests_v2"),
        embedding_dim=2,
        null_repetitions=2,
        steps=300,
        restarts=1,
        output_dir=Path("outputs"),
    )
    rows = run_exp186(config)
    output_path = write_exp186_outputs(rows, config.output_dir)
    figure_paths = save_exp186_figures(rows, config.output_dir)
    print(f"exp186 wrote {output_path} and {len(figure_paths)} figures")


def run_v2_stricter_criteria_diagnostics() -> None:
    config = Exp187Config(
        manifest_dir=Path("outputs/manifests_v2"),
        dims=(1, 2),
        steps=300,
        restarts=1,
        output_dir=Path("outputs"),
    )
    rows = run_exp187(config)
    output_path = write_exp187_outputs(rows, config.output_dir)
    figure_paths = save_exp187_figures(rows, config.output_dir)
    print(f"exp187 wrote {output_path} and {len(figure_paths)} figures")


def run_v2_failed_accounting() -> None:
    config = Exp188Config(
        manifest_dir=Path("outputs/manifests_v2"),
        output_dir=Path("outputs"),
    )
    rows = run_exp188(config)
    output_path = write_exp188_outputs(rows, config.output_dir)
    figure_paths = save_exp188_figures(rows, config.output_dir)
    print(f"exp188 wrote {output_path} and {len(figure_paths)} figures")


def run_v2_coverage_metrics() -> None:
    config = Exp189Config(
        manifest_dir=Path("outputs/manifests_v2"),
        output_dir=Path("outputs"),
    )
    rows = run_exp189(config)
    output_path = write_exp189_outputs(rows, config.output_dir)
    print(f"exp189 wrote {output_path}")


def run_v2_restart_latent_order_stability() -> None:
    config = Exp190Config(
        manifest_dir=Path("outputs/manifests_v2"),
        embedding_dim=2,
        restart_count=3,
        steps=300,
        output_dir=Path("outputs"),
    )
    restart_rows, latent_rows = run_exp190(config)
    paths = write_exp190_outputs(restart_rows, latent_rows, config.output_dir)
    figure_paths = save_exp190_figures(restart_rows, latent_rows, config.output_dir)
    print(f"exp190 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp190 wrote {len(figure_paths)} figures")


def run_v2_no_retuning_audit() -> None:
    config = Exp191Config(output_dir=Path("outputs"))
    rows = run_exp191(config)
    output_path = write_exp191_outputs(rows, config.output_dir)
    print(f"exp191 wrote {output_path}")


def run_v2_required_metric_aggregation() -> None:
    config = Exp192Config(output_dir=Path("outputs"))
    metric_rows, completeness_rows = run_exp192(config)
    paths = write_exp192_outputs(metric_rows, completeness_rows, config.output_dir)
    figure_paths = save_exp192_figures(completeness_rows, config.output_dir)
    print(f"exp192 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp192 wrote {len(figure_paths)} figures")


def run_v2_diagnostic_complete_bundle_report() -> None:
    config = Exp193Config(output_dir=Path("outputs"))
    rows = run_exp193(config)
    output_path = write_exp193_outputs(rows, config.output_dir)
    print(f"exp193 wrote {output_path}")


def run_v2_manifest_generation_final_sanity() -> None:
    rows = run_exp194()
    output_path = write_exp194_outputs(rows)
    print(f"exp194 wrote {output_path}")


def run_v2_carry_forward_exact_sanity() -> None:
    rows = run_exp195()
    output_path = write_exp195_outputs(rows)
    print(f"exp195 wrote {output_path}")


def run_v2_bundle_input_audit() -> None:
    config = Exp196Config(output_dir=Path("outputs"))
    rows = run_exp196(config)
    output_path = write_exp196_outputs(rows, config.output_dir)
    print(f"exp196 wrote {output_path}")


def run_v2_carry_forward_decision() -> None:
    config = Exp197Config(output_dir=Path("outputs"))
    metrics, decisions = run_exp197(config)
    paths = write_exp197_outputs(metrics, decisions, config.output_dir)
    figure_paths = save_exp197_figures(decisions, config.output_dir)
    print(f"exp197 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp197 wrote {len(figure_paths)} figures")


def run_v2_carry_forward_threshold_sensitivity() -> None:
    config = Exp198Config(
        output_dir=Path("outputs"),
        heldout_thresholds=(0.20, 0.25),
        null_gap_thresholds=(0.05, 0.10),
        stricter_pass_thresholds=(0.25, 0.50),
    )
    rows = run_exp198(config)
    output_path = write_exp198_outputs(rows, config.output_dir)
    figure_paths = save_exp198_figures(rows, config.output_dir)
    print(f"exp198 wrote {output_path} and {len(figure_paths)} figures")


def run_v2_carry_forward_registry_export() -> None:
    config = Exp199Config(
        output_dir=Path("outputs"),
        manifest_dir=Path("outputs/manifests_v2"),
    )
    registry_path, rows = run_exp199(config)
    output_path = write_exp199_outputs(rows, config.output_dir)
    print(f"exp199 wrote {output_path}")
    if registry_path is not None:
        print(f"exp199 wrote {registry_path}")


def run_v2_stress_test_handoff_plan() -> None:
    config = Exp200Config(output_dir=Path("outputs"))
    rows = run_exp200(config)
    output_path = write_exp200_outputs(rows, config.output_dir)
    print(f"exp200 wrote {output_path}")


def run_v2_failed_provisional_accounting() -> None:
    config = Exp201Config(output_dir=Path("outputs"))
    rows = run_exp201(config)
    output_path = write_exp201_outputs(rows, config.output_dir)
    figure_paths = save_exp201_figures(rows, config.output_dir)
    print(f"exp201 wrote {output_path} and {len(figure_paths)} figures")


def run_v2_carry_forward_failure_decomposition() -> None:
    config = Exp202Config(output_dir=Path("outputs"))
    rows, summary = run_exp202(config)
    paths = write_exp202_outputs(rows, summary, config.output_dir)
    figure_paths = save_exp202_figures(summary, config.output_dir)
    print(f"exp202 wrote {', '.join(str(path) for path in paths)}")
    print(f"exp202 wrote {len(figure_paths)} figures")


def run_v2_carry_forward_no_retuning_audit() -> None:
    config = Exp203Config(output_dir=Path("outputs"))
    rows = run_exp203(config)
    output_path = write_exp203_outputs(rows, config.output_dir)
    print(f"exp203 wrote {output_path}")


def run_v2_carry_forward_report_card() -> None:
    config = Exp204Config(output_dir=Path("outputs"))
    rows = run_exp204(config)
    output_path = write_exp204_outputs(rows, config.output_dir)
    print(f"exp204 wrote {output_path}")


def run_v2_carry_forward_final_sanity() -> None:
    rows = run_exp205()
    output_path = write_exp205_outputs(rows)
    print(f"exp205 wrote {output_path}")


def run_v2_blocking_analysis_exact_sanity() -> None:
    run_exp206_main()


def run_v2_blocked_root_cause_audit() -> None:
    run_exp207_main()


def run_v2_criterion_margin_report() -> None:
    run_exp208_main()


def run_v2_structural_vs_measured_blocking() -> None:
    run_exp209_main()


def run_v3_manifest_family_design() -> None:
    run_exp210_main()


def run_v3_preregistration_export() -> None:
    run_exp211_main()


def run_v3_no_execution_audit() -> None:
    run_exp212_main()


def run_v2_blocked_decision_report_card() -> None:
    run_exp213_main()


def run_v2_blocked_v3_preregistration_final_sanity() -> None:
    run_exp214_main()


def run_protocol_metadata_exact_sanity() -> None:
    run_exp215_main()


def run_response_profile_invariance_exact_sanity() -> None:
    run_exp216_main()


def run_current_profile_protocol_invariance_audit() -> None:
    run_exp217_main()


def run_v3_protocol_invariant_patch_design() -> None:
    run_exp218_main()


def run_v3_protocol_patched_preregistration_export() -> None:
    run_exp219_main()


def run_v3_protocol_patch_audit() -> None:
    run_exp220_main()


def run_protocol_invariance_language_audit() -> None:
    run_exp221_main()


def run_protocol_patch_no_execution_audit() -> None:
    run_exp222_main()


def run_protocol_invariance_final_sanity() -> None:
    run_exp223_main()


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
    run_echo_response_signature_exact_sanity()
    run_layered_echo_response_order_recovery()
    run_echo_response_signature_coarse_protocol_stability()
    run_echo_response_shortcut_robust_core()
    run_echo_response_reference_protocol_dependence()
    run_echo_response_order_precondition_diagnostics()
    run_echo_response_signature_stability_exact_sanity()
    run_echo_order_semantics_exact_sanity()
    run_gated_echo_protocol_comparison()
    run_echo_delay_spacing_tie_resolution()
    run_response_order_scalar_representability()
    run_echo_terminology_audit()
    run_response_representability_exact_sanity()
    run_response_order_underdetermination_exact_sanity()
    run_single_reference_response_order_underdetermination()
    run_multi_reference_response_profile_diagnostics()
    run_response_representability_requirement_table()
    run_response_profile_stability_under_protocol_variants()
    run_response_profile_requirement_exact_sanity()
    run_pairwise_response_profile_exact_sanity()
    run_pairwise_response_protocol_comparison()
    run_pairwise_response_null_baselines()
    run_pairwise_response_protocol_variant_stability()
    run_pairwise_response_missing_data_sensitivity()
    run_pairwise_response_protocol_choice_dependence()
    run_pairwise_response_null_admissibility_exact_sanity()
    run_response_constraint_pool_exact_sanity()
    run_response_constraint_heldout_protocol_validation()
    run_response_constraint_bootstrap_stability()
    run_response_constraint_null_separation()
    run_response_constraint_pool_coverage()
    run_response_constraint_validation_gate_summary()
    run_response_constraint_validation_exact_sanity()
    run_response_handoff_exact_sanity()
    run_response_handoff_threshold_sensitivity()
    run_response_handoff_protocol_selection()
    run_response_handoff_manifest_export()
    run_response_handoff_failure_catalog()
    run_response_handoff_preregistration_rules()
    run_response_handoff_manifest_read_exact_sanity()
    run_manifest_representability_exact_sanity()
    run_frozen_manifest_ordinal_representation()
    run_frozen_manifest_representation_nulls()
    run_frozen_manifest_fit_stability()
    run_frozen_manifest_dimension_complexity_curve()
    run_failed_manifest_no_fit_controls()
    run_frozen_manifest_representation_summary()
    run_manifest_representation_no_metric_exact_sanity()
    run_manifest_family_exact_sanity()
    run_manifest_family_fit_comparison()
    run_manifest_family_null_taxonomy()
    run_manifest_family_stricter_criteria()
    run_manifest_family_failed_manifest_accounting()
    run_manifest_family_no_retuning_audit()
    run_manifest_family_report_card()
    run_manifest_family_comparison_exact_sanity()
    run_cross_family_robustness_exact_sanity()
    run_cross_family_robustness_criteria_table()
    run_cross_family_robustness_decision()
    run_cross_family_robustness_threshold_sensitivity()
    run_carry_forward_registry_export()
    run_cross_family_failed_provisional_accounting()
    run_stress_test_handoff_plan()
    run_cross_family_robustness_report_card()
    run_cross_family_robustness_final_sanity()
    run_carry_forward_failure_decomposition_exact_sanity()
    run_carry_forward_failure_decomposition()
    run_cross_family_diagnostic_completeness_audit()
    run_stress_test_stop_condition_audit()
    run_upstream_remediation_design_table()
    run_missing_metric_impact_report()
    run_failure_decomposition_no_retuning_audit()
    run_carry_forward_failure_report_card()
    run_carry_forward_failure_final_sanity()
    run_remediation_plan_exact_sanity()
    run_diagnostic_complete_schema_export()
    run_failure_to_remediation_mapping()
    run_new_manifest_family_design_v2()
    run_preregistered_remediation_plan_export()
    run_future_manifest_run_spec()
    run_remediation_no_execution_audit()
    run_remediation_plan_report_card()
    run_remediation_plan_final_sanity()
    run_v2_manifest_spec_exact_sanity()
    run_v2_manifest_generation()
    run_v2_family_fit_diagnostics()
    run_v2_null_taxonomy_diagnostics()
    run_v2_stricter_criteria_diagnostics()
    run_v2_failed_accounting()
    run_v2_coverage_metrics()
    run_v2_restart_latent_order_stability()
    run_v2_no_retuning_audit()
    run_v2_required_metric_aggregation()
    run_v2_diagnostic_complete_bundle_report()
    run_v2_manifest_generation_final_sanity()
    run_v2_carry_forward_exact_sanity()
    run_v2_bundle_input_audit()
    run_v2_carry_forward_decision()
    run_v2_carry_forward_threshold_sensitivity()
    run_v2_carry_forward_registry_export()
    run_v2_stress_test_handoff_plan()
    run_v2_failed_provisional_accounting()
    run_v2_carry_forward_failure_decomposition()
    run_v2_carry_forward_no_retuning_audit()
    run_v2_carry_forward_report_card()
    run_v2_carry_forward_final_sanity()
    run_v2_blocking_analysis_exact_sanity()
    run_v2_blocked_root_cause_audit()
    run_v2_criterion_margin_report()
    run_v2_structural_vs_measured_blocking()
    run_v3_manifest_family_design()
    run_v3_preregistration_export()
    run_v3_no_execution_audit()
    run_v2_blocked_decision_report_card()
    run_v2_blocked_v3_preregistration_final_sanity()
    run_protocol_metadata_exact_sanity()
    run_response_profile_invariance_exact_sanity()
    run_current_profile_protocol_invariance_audit()
    run_v3_protocol_invariant_patch_design()
    run_v3_protocol_patched_preregistration_export()
    run_v3_protocol_patch_audit()
    run_protocol_invariance_language_audit()
    run_protocol_patch_no_execution_audit()
    run_protocol_invariance_final_sanity()


if __name__ == "__main__":
    main()
