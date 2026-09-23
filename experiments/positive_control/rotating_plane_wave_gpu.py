"""Optional float64 GPU port of the four guarded relation-query methods.

Sampling remains in the original CPU generator. This module is never imported
by the estimator. The array formulas below are copied verbatim from the frozen
CPU reference; only the array and polynomial backends differ.
"""

import cupy as np
from cupy.polynomial.polynomial import polyval


class RotatingWaveGPU:
    def __init__(self, amplitude, rotation_rate, coefficients):
        if not (0 <= amplitude <= 0.75 and abs(rotation_rate) <= 0.5):
            raise ValueError("Outside the bounded rotating-wave generator family")
        self.amplitude = float(amplitude)
        self.rotation_rate = float(rotation_rate)
        self.coefficients = np.asarray(coefficients, dtype=np.float64)
        if self.coefficients.shape != (13, 6):
            raise ValueError("Expected the frozen 13-term propagator coefficients")

    def packed_relations(self, landmarks, events):
        if len(events) % 8:
            raise ValueError("Measurement batches must be divisible by eight")
        future, past = self.relations(
            np.asarray(landmarks, dtype=np.float64)[:, None],
            np.asarray(events, dtype=np.float64)[None],
        )
        return (
            np.asnumpy(np.packbits(future, axis=1)),
            np.asnumpy(np.packbits(past, axis=1)),
        )

    def blocks(self, length):
        """Return A, B/length, C, D blocks without subtracting close cosines.

        Taylor tails are <1e-24 in this bounded domain. Floating roundoff is
        handled separately by the action guard; this is not an interval library.
        """
        length = np.asarray(length, dtype=float)
        if not np.isfinite(length).all() or np.any(abs(length) > 1):
            raise ValueError("Finite propagator lengths in [-1,1] required")
        cp, cm, sp, sm, hc, hs = polyval(length**2, self.coefficients)
        w, a = self.rotation_rate, self.amplitude
        aa = np.zeros(length.shape + (2, 2))
        bb, cc, dd = (np.zeros_like(aa) for _ in range(3))
        aa[..., 0, 0] = dd[..., 0, 0] = cp
        aa[..., 1, 1] = dd[..., 1, 1] = cm
        aa[..., 0, 1], aa[..., 1, 0] = w * length * sm, -w * length * sp
        dd[..., 0, 1], dd[..., 1, 0] = w * length * sp, -w * length * sm
        bb[..., 0, 0] = sp - 2 * w * w * length**2 * hs
        bb[..., 1, 1] = sm - 2 * w * w * length**2 * hs
        bb[..., 0, 1], bb[..., 1, 0] = 2 * w * length * hc, -2 * w * length * hc
        cc[..., 0, 0], cc[..., 1, 1] = a * length * sp, -a * length * sm
        return aa, bb, cc, dd

    @staticmethod
    def inverse(matrix):
        det = (
            matrix[..., 0, 0] * matrix[..., 1, 1]
            - matrix[..., 0, 1] * matrix[..., 1, 0]
        )
        if np.any(det <= 0) or not np.isfinite(det).all():
            raise ArithmeticError("Nonpositive scaled Jacobi determinant")
        result = np.empty_like(matrix)
        result[..., 0, 0], result[..., 1, 1] = matrix[..., 1, 1], matrix[..., 0, 0]
        result[..., 0, 1], result[..., 1, 0] = -matrix[..., 0, 1], -matrix[..., 1, 0]
        return result / det[..., None, None]

    def action(self, length, left, right):
        """Exact-action formula and a conservative floating arithmetic guard."""
        length = np.asarray(length)
        if np.any(length <= 0):
            raise ValueError("Positive action lengths required")
        aa, bb, _, dd = self.blocks(length)
        inverse = self.inverse(bb)
        first, last = inverse @ aa, dd @ inverse
        terms = np.stack(
            [
                np.einsum("...i,...ij,...j->...", left, first, left),
                -2 * np.einsum("...i,...ij,...j->...", left, inverse, right),
                np.einsum("...i,...ij,...j->...", right, last, right),
            ]
        )
        value = terms.sum(axis=0) / (2 * length)
        norm = lambda m: np.max(np.sum(abs(m), axis=-1), axis=-1)  # noqa: E731
        g, b, an, dn = norm(inverse), norm(bb), norm(aa), norm(dd)
        eps, block_error = np.finfo(float).eps, 1e-12
        residual = 2 * block_error * g + 64 * eps * (1 + b * g)
        if np.any(residual >= 0.5):
            raise ArithmeticError("Jacobi inversion error cannot be controlled")
        inverse_error = g * residual / (1 - residual)
        efirst = (
            inverse_error * an
            + 2 * block_error * (g + inverse_error)
            + 64 * eps * g * an
        )
        elast = (
            inverse_error * dn
            + 2 * block_error * (g + inverse_error)
            + 64 * eps * g * dn
        )
        lp, rp = np.max(abs(left), axis=-1), np.max(abs(right), axis=-1)
        error = (
            lp * lp * efirst + 2 * lp * rp * inverse_error + rp * rp * elast
        ) / length
        error += 64 * eps * (1 + abs(terms).sum(axis=0)) / length
        return value, np.nextafter(4 * error, np.inf)

    def relations(self, left, right):
        """Broadcastable pair queries; unresolved margins stop the generator."""
        du, dv = right[..., 0] - left[..., 0], right[..., 1] - left[..., 1]
        future = du > 0
        first = np.where(future[..., None], left[..., 2:], right[..., 2:])
        last = np.where(future[..., None], right[..., 2:], left[..., 2:])
        cost, error = self.action(np.where(du == 0, 1.0, abs(du)), first, last)
        margin = -np.sign(du) * dv - cost
        guard = error + 256 * np.finfo(float).eps * (1 + abs(dv))
        near = (abs(margin) <= guard) & (du != 0)
        if near.any():
            raise ArithmeticError(f"{int(near.sum())} unresolved rotating-wave margins")
        related = (margin > guard) & (du != 0)
        return related & future, related & ~future
