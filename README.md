# causal-spacetime-lab

`causal-spacetime-lab` is a small Python research simulation project for
testing operational reconstructions of spacetime quantities from causal order
and related information-accessibility structure.

The project is scientifically conservative: these simulations are sanity checks
for reconstruction procedures and known relativistic behavior. They do not prove
a new theory of spacetime.

## Current Milestone

The current milestone upgrades the project into a lightweight validation suite
for 1+1D special-relativistic and causal-set reconstruction experiments:

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
  -> operational time, distance, dimension, and spacetime-like geometry
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
