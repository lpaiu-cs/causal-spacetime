# AGENTS.md

## Project identity

This repository is a research-oriented simulation lab for studying whether spacetime quantities can be operationally reconstructed from causal information structure.

The working research hypothesis is:

> Time and space are not derived from velocity as displacement divided by time. Instead, operational time and distance can be reconstructed from primitive causal or null information-accessibility relations, with metric scale requiring additional structure such as event density, clocks, or free-fall information.

This project should remain scientifically conservative. Do not overclaim that the simulations prove a new theory of physics. The simulations should test reconstruction procedures, failure cases, and known relativistic behavior.

## Main goals

Build a Python simulation codebase for:

1. Radar coordinate decomposition in special relativity.
2. Lorentz kinematics and length contraction as observer-dependent event selection.
3. Causal set sprinkling in Minkowski spacetime.
4. Timelike proper-time reconstruction from longest chains.
5. Interval-volume reconstruction from Alexandrov interval cardinality.
6. Spacelike-distance reconstruction experiments, initially as exploratory code.
7. Finite-speed lattice counterexamples showing that finite signal speed alone does not imply Lorentzian spacetime.

## Scientific constraints

Use natural units `c = 1` by default.

Distinguish carefully between:

- signal speed,
- causal order,
- null relation,
- observer protocol,
- proper time,
- coordinate time,
- spatial distance,
- metric scale.

Avoid language implying that speed is primitive in the usual `v = dx/dt` sense.

Prefer this language:

- "null information-accessibility structure"
- "causal order"
- "operational reconstruction"
- "radar decomposition"
- "event density"
- "causal interval cardinality"

Avoid this language unless explicitly discussed as a rejected formulation:

- "time is information-transfer speed"
- "space is information-transfer speed"
- "this proves spacetime is information"

## Implementation standards

Use Python 3.11 or later.

Use:

- `numpy`
- `scipy` if needed
- `matplotlib`
- `networkx` only if useful and not performance-critical
- `pytest`
- `ruff`

Do not use heavyweight dependencies unless justified.

Keep the package modular.

Suggested structure:

```text
causal-spacetime-lab/
  README.md
  pyproject.toml
  AGENTS.md
  src/
    causal_spacetime_lab/
      __init__.py
      constants.py
      lorentz.py
      radar.py
      sprinkling.py
      causal.py
      chains.py
      intervals.py
      lattice.py
      metrics.py
      plotting.py
  experiments/
    exp01_radar_decomposition.py
    exp02_lorentz_length_contraction.py
    exp03_causalset_timelike_reconstruction.py
    exp04_interval_volume_reconstruction.py
    exp05_finite_speed_lattice_counterexample.py
  tests/
    test_lorentz.py
    test_radar.py
    test_sprinkling.py
    test_causal.py
    test_chains.py
  docs/
    research_notes.md
    references.md