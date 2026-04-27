from __future__ import annotations

import numpy as np

from causal_spacetime_lab.chains import longest_chain_indices, longest_chain_length


def test_longest_chain_on_small_manual_dag() -> None:
    C = np.zeros((5, 5), dtype=bool)
    C[0, 1] = True
    C[1, 3] = True
    C[3, 4] = True
    C[0, 2] = True
    C[2, 4] = True
    C[0, 3] = True
    C[1, 4] = True
    C[0, 4] = True

    assert longest_chain_length(C) == 4
    assert longest_chain_length(C, start=0, end=4) == 4


def test_longest_chain_returns_zero_when_no_start_end_chain_exists() -> None:
    C = np.zeros((3, 3), dtype=bool)
    C[0, 1] = True

    assert longest_chain_length(C, start=1, end=2) == 0


def test_longest_chain_indices_returns_one_maximal_path() -> None:
    C = np.zeros((5, 5), dtype=bool)
    C[0, 1] = True
    C[1, 3] = True
    C[3, 4] = True
    C[0, 2] = True
    C[2, 4] = True
    C[0, 3] = True
    C[1, 4] = True
    C[0, 4] = True

    assert longest_chain_indices(C, start=0, end=4).tolist() == [0, 1, 3, 4]
