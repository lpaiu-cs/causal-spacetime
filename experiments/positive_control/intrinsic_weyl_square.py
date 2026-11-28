"""Curvature-squared world-function and vacuum volume-clock correction.

No new fitted parameters. This is NOT a full sixth-degree model: covariant
second derivatives of curvature remain omitted. See square_correction.md.
"""

import numpy as np
from intrinsic_weyl import WEYL

SIGN = np.array([-1.0, 1.0, 1.0, 1.0])
RAISING = np.einsum("a,b,c,d->abcd", SIGN, SIGN, SIGN, SIGN)


def correction(points, coefficients, pairs, *, volume_clock=True):
    """Return value and analytic left/right/coefficient derivatives."""
    m = len(pairs)
    if not np.any(coefficients[:10]):
        return (
            np.zeros(m),
            np.zeros((m, 4)),
            np.zeros((m, 4)),
            np.zeros((m, len(coefficients))),
        )
    x, y = points[pairs[:, 0]], points[pairs[:, 1]]
    tensor = np.einsum("k,kabcd->abcd", coefficients[:10], WEYL)
    u = np.einsum("abcd,pa,pb,pc->pd", tensor, x, y, x, optimize=True)
    v = np.einsum("abcd,pa,pb,pc->pd", tensor, x, y, y, optimize=True)
    value = np.sum((4 * u * v - u * u - v * v) * SIGN, axis=1) / 45
    gu, gv = (4 * v - 2 * u) * SIGN / 45, (4 * u - 2 * v) * SIGN / 45
    ux = np.einsum("jbcd,pb,pc->pjd", tensor, y, x, optimize=True) + np.einsum(
        "abjd,pa,pb->pjd", tensor, x, y, optimize=True
    )
    uy = np.einsum("ajcd,pa,pc->pjd", tensor, x, x, optimize=True)
    vx = np.einsum("jbcd,pb,pc->pjd", tensor, y, y, optimize=True)
    vy = np.einsum("ajcd,pa,pc->pjd", tensor, x, y, optimize=True) + np.einsum(
        "abjd,pa,pb->pjd", tensor, x, y, optimize=True
    )
    dl = np.einsum("pd,pjd->pj", gu, ux) + np.einsum("pd,pjd->pj", gv, vx)
    dr = np.einsum("pd,pjd->pj", gu, uy) + np.einsum("pd,pjd->pj", gv, vy)
    uc = np.einsum("kabcd,pa,pb,pc->pkd", WEYL, x, y, x, optimize=True)
    vc = np.einsum("kabcd,pa,pb,pc->pkd", WEYL, x, y, y, optimize=True)
    dc = np.einsum("pd,pkd->pk", gu, uc) + np.einsum("pd,pkd->pk", gv, vc)
    if volume_clock:
        # The final row is the fixed outer-tip displacement in normalized units.
        delta = np.vstack((x - y, [1, 0, 0, 0]))
        squared = np.sum(delta * delta * SIGN, axis=1)
        ds = 2 * delta * SIGN
        electric = np.einsum("acbd,pc,pd->pab", tensor, delta, delta, optimize=True)
        raised = electric * SIGN[None, :, None] * SIGN[None, None, :]
        invariant = np.sum(tensor * tensor * RAISING) / 8
        invariant_c = np.einsum("kabcd,abcd,abcd->k", WEYL, tensor, RAISING) / 4
        e_delta = np.einsum("ajbd,pd->pjab", tensor, delta) + np.einsum(
            "acbj,pc->pjab", tensor, delta
        )
        e_c = np.einsum("kacbd,pc,pd->pkab", WEYL, delta, delta, optimize=True)
        q = (
            20 * np.sum(electric * raised, axis=(1, 2)) - 3 * squared**2 * invariant
        ) / 10080
        q_delta = (
            40 * np.einsum("pab,pjab->pj", raised, e_delta)
            - 6 * squared[:, None] * ds * invariant
        ) / 10080
        q_c = (
            40 * np.einsum("pab,pkab->pk", raised, e_c)
            - 3 * squared[:, None] ** 2 * invariant_c
        ) / 10080
        value += 0.5 * squared[:-1] * (q[:-1] - q[-1])
        derivative = 0.5 * (
            ds[:-1] * (q[:-1] - q[-1])[:, None] + squared[:-1, None] * q_delta[:-1]
        )
        dl += derivative
        dr -= derivative
        dc += 0.5 * squared[:-1, None] * (q_c[:-1] - q_c[-1])
    if len(coefficients) > 10:
        dc = np.column_stack((dc, np.zeros((m, len(coefficients) - 10))))
    return value, dl, dr, dc


def selfcheck():
    # Complex-step verification checks every input derivative on a generic Weyl
    # tensor with electric AND magnetic components, not a plane-wave type prior.
    points = np.random.default_rng(1701).normal(size=(8, 4)) * 0.15
    c = np.linspace(-0.7, 0.8, 34)
    pairs = np.column_stack(np.triu_indices(len(points), 1))
    maximum = 0.0
    for volume_clock in (False, True):
        value, dl, dr, dc = correction(points, c, pairs, volume_clock=volume_clock)
        for k in range(34):
            changed = c.astype(complex)
            changed[k] += 1e-25j
            numeric = (
                correction(points, changed, pairs, volume_clock=volume_clock)[0].imag
                / 1e-25
            )
            maximum = max(maximum, float(np.max(abs(numeric - dc[:, k]))))
        for i in range(len(points)):
            for a in range(4):
                changed = points.astype(complex)
                changed[i, a] += 1e-25j
                numeric = (
                    correction(changed, c, pairs, volume_clock=volume_clock)[0].imag
                    / 1e-25
                )
                expected = (pairs[:, 0] == i) * dl[:, a] + (pairs[:, 1] == i) * dr[:, a]
                maximum = max(maximum, float(np.max(abs(numeric - expected))))
        assert np.isfinite(value).all()
        assert (
            np.max(
                abs(
                    correction(points, c, pairs[:, ::-1], volume_clock=volume_clock)[0]
                    - value
                )
            )
            < 1e-14
        )
    assert maximum < 1e-12
    assert all(np.all(x == 0) for x in correction(points, np.zeros(34), pairs))
    axial = np.array([[-0.5, 0, 0, 0], [-0.25, 0, 0, 0]])
    for electric_only, magnetic_only in ((True, False), (False, True), (True, True)):
        coefficients = np.r_[c[:5] * electric_only, c[5:10] * magnetic_only]
        electric = np.einsum("k,kij->ij", coefficients[:5], WEYL[:5, 0, 1:, 0, 1:])
        magnetic = np.einsum("k,kij->ij", coefficients[5:], WEYL[:5, 0, 1:, 0, 1:])
        expected = (17 * np.sum(electric**2) + 3 * np.sum(magnetic**2)) / 10080
        got = correction(axial, coefficients, np.array([[0, 1]]))[0][0]
        assert abs(got / (-0.5 * 0.25**2 * (0.25**4 - 1)) - expected) < 1e-15
    return {
        "maximum_complex_step_jacobian_error": maximum,
        "pair_symmetry_and_flat_zero": "passed",
        "electric_and_magnetic_volume_weights": "passed",
    }


if __name__ == "__main__":
    import json

    print(json.dumps(selfcheck()))
