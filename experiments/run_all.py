"""Run the lightweight validation and demonstration experiment suite."""

from __future__ import annotations

from pathlib import Path

import exp02_lorentz_length_contraction as exp02
import exp03_causalset_timelike_reconstruction as exp03
import exp05_finite_speed_lattice_counterexample as exp05
import exp06_spacelike_distance_reconstruction as exp06
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


def main() -> None:
    run_lorentz_length_contraction()
    run_legacy_timelike_reconstruction()
    run_finite_speed_lattice_counterexample()
    run_spacelike_proxy_experiment()
    run_timelike_pair_convergence()
    run_probe_pair_statistical_calibration()
    run_dimension_reconstruction()
    run_discrete_observer_radar_reconstruction()


if __name__ == "__main__":
    main()
