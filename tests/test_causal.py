from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1, causal_matrix_minkowski


def test_causal_matrix_has_no_self_causality() -> None:
    events = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
        ]
    )
    C = causal_matrix_1p1(events)

    assert not np.any(np.diag(C))


def test_causal_relation_matches_hand_coded_event_pairs() -> None:
    events = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.5, 1.0],
            [-1.0, 0.0],
        ]
    )

    C = causal_matrix_1p1(events)

    assert C[0, 1]
    assert C[0, 2]
    assert not C[0, 3]
    assert not C[1, 0]
    assert C[4, 0]


def test_generalized_causal_matrix_matches_hand_coded_2p1_events() -> None:
    events = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.5, 0.5],
            [1.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
        ]
    )

    C = causal_matrix_minkowski(events)

    assert C[0, 1]
    assert not C[0, 2]
    assert C[1, 3]
    assert not C[2, 3]

