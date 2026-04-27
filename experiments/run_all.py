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


if __name__ == "__main__":
    main()
