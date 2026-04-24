from __future__ import annotations

import numpy as np

from causal_spacetime_lab.intervals import alexandrov_interval_indices


def test_alexandrov_interval_indices_returns_events_between_endpoints() -> None:
    C = np.zeros((4, 4), dtype=bool)
    C[0, 1] = True
    C[0, 2] = True
    C[0, 3] = True
    C[1, 3] = True
    C[2, 3] = True

    indices = alexandrov_interval_indices(C, 0, 3)

    assert set(indices.tolist()) == {1, 2}

