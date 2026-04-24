# causal-spacetime-lab

`causal-spacetime-lab` is a small Python research simulation project for
testing operational reconstructions of spacetime quantities from causal order
and related information-accessibility structure.

The project is scientifically conservative: these simulations are sanity checks
for reconstruction procedures and known relativistic behavior. They do not prove
a new theory of spacetime.

## Current Milestone

The first milestone implements a minimal 1+1D Minkowski causal-set experiment:

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

## What The Result Means

The experiment checks whether timelike separation in a 1+1D Minkowski causal
diamond can be reconstructed using causal interval cardinality once an event
density is specified. The longest chain is also reported as a causal-order
observable with a simple 1+1D asymptotic normalization.

This result does not show that spacetime is made of information, and it does not
derive metric scale from causal order alone. The interval-count estimate uses
event density as additional structure, which is exactly the point being tested.

