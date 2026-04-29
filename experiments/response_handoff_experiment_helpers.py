"""Shared helpers for response handoff experiments."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def deterministic_handoff_profile() -> EchoResponseProfile:
    """Return a small deterministic response profile for exact checks."""

    return EchoResponseProfile(
        target_event_ids=np.arange(6, dtype=int),
        protocol_labels=["a", "b", "c", "d"],
        delay_rank_matrix=np.asarray(
            [
                [1, 1, 2, 2],
                [2, 2, 3, 3],
                [3, 3, 4, 4],
                [5, 5, 6, 6],
                [8, 8, 9, 9],
                [13, 13, 14, 14],
            ],
            dtype=int,
        ),
        reachable_matrix=np.ones((6, 4), dtype=bool),
    )


def failure_reason_rows(
    base: dict[str, float | str],
    failed_reasons: list[str],
) -> list[dict[str, float | str]]:
    """Expand failed reasons for simple CSV summaries."""

    if not failed_reasons:
        return [{**base, "failure_reason": "none"}]
    return [{**base, "failure_reason": reason} for reason in failed_reasons]
