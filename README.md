# causal-spacetime-lab

`causal-spacetime-lab` is a small Python research simulation project for
testing operational reconstructions of spacetime quantities from causal order
and related information-accessibility structure.

The project is scientifically conservative: these simulations are sanity checks
for reconstruction procedures and known relativistic behavior. They do not prove
a new theory of spacetime.

## Current Milestone

The first milestone implements a minimal 1+1D Minkowski causal-set experiment:

- decompose events into radar coordinates for stationary and inertial observers,
- test Lorentz length contraction as lab-simultaneous endpoint event selection,
- sample events uniformly inside a causal diamond,
- build the causal order matrix,
- compute longest causal chains,
- count Alexandrov interval elements,
- estimate timelike proper time from interval cardinality when event density is
  supplied as extra metric scale information.

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

## Run The First Experiment

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

## What The Result Means

The experiment checks whether timelike separation in a 1+1D Minkowski causal
diamond can be reconstructed using causal interval cardinality once an event
density is specified. The longest chain is also reported as a causal-order
observable with a simple 1+1D asymptotic normalization.

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
