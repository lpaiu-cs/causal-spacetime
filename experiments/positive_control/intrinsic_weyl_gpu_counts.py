"""Exact membership second moments on GPU; no generator or truth imports."""

import cupy as cp
import numpy as np


def clock_second_moment(after, before, a, b):
    """Same integer-weighted sum as the CPU helper, with float64 GPU GEMM.

    Every input, product and positive partial sum is an integer below 2**53.
    Moving this sum to GPU therefore adds no floating approximation. The CPU
    grouping and the resulting NumPy covariance/optimizer remain unchanged.
    """
    landmarks, bytes_per_landmark = after.shape
    n = 8 * bytes_per_landmark
    if n >= 2**53:
        raise ValueError("Float64 cannot exactly accumulate this event count")
    # ponytail: retain the CPU grouping; external grouping is needed only when
    # the packed archive and signatures no longer fit the measured RAM budget.
    signatures = np.empty((n, (2 * landmarks + 7) // 8), dtype=np.uint8)
    for start in range(0, bytes_per_landmark, 4096):
        end = min(start + 4096, bytes_per_landmark)
        bits = np.unpackbits(
            np.vstack((after[:, start:end], before[:, start:end])), axis=1
        )
        signatures[8 * start : 8 * end] = np.packbits(bits, axis=0).T
    width = signatures.shape[1]
    unique, weights = np.unique(
        signatures.view(np.dtype(f"V{width}")).ravel(), return_counts=True
    )
    del signatures
    unique = unique.view(np.uint8).reshape(-1, width)
    second = cp.zeros((len(a), len(a)), dtype=cp.float64)
    for start in range(0, len(unique), 4096):
        bits = np.unpackbits(unique[start : start + 4096], axis=1)
        membership = cp.asarray(
            (bits[:, a] & bits[:, landmarks + b]).T, dtype=cp.float64
        )
        weight = cp.asarray(weights[start : start + 4096], dtype=cp.float64)
        second += (membership * weight) @ membership.T
    return cp.asnumpy(second)


def selfcheck():
    from intrinsic_weyl_joint import _clock_second_moment

    rng = np.random.default_rng(719)
    for landmarks, n, distinct in ((13, 8192, 7), (9, 512, 512), (11, 256, 1)):
        patterns = rng.integers(0, 2, (2 * landmarks, distinct), dtype=np.uint8)
        bits = patterns[:, np.arange(n) % distinct]
        a, b = np.arange(landmarks - 1), np.arange(1, landmarks)
        after = np.packbits(bits[:landmarks], axis=1)
        before = np.packbits(bits[landmarks:], axis=1)
        membership = (bits[a] & bits[landmarks + b]).astype(float)
        expected = membership @ membership.T
        assert np.array_equal(clock_second_moment(after, before, a, b), expected)
        assert np.array_equal(_clock_second_moment(after, before, a, b), expected)
    return {"direct_and_grouped_second_moments_bit_exact": True, "cases": 3}


if __name__ == "__main__":
    import json

    print(json.dumps(selfcheck()))
