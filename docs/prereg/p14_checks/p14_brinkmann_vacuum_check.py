"""P14 design check 1: the Brinkmann quadratic profile is exact vacuum.

For  ds^2 = 2 du dv + A(u) (x^2 - y^2) du^2 + dx^2 + dy^2  with A(u) an
UNDETERMINED function, three exact statements carry the design
(p14_weyl_curvature.md Section 4), and this script is what [VERIFIED]
refers to -- an output pasted into a document is not a check, per the
repository's own reproducibility rule:

  1. det g = -1 identically, so sqrt(-g) = 1: the volume form is
     Minkowski's for every profile, and a uniform-coordinate Poisson
     sprinkling is the covariant Poisson process exactly.
  2. Ricci = 0 identically: vacuum, no matter term anywhere.
  3. Riemann != 0: eight non-vanishing components, all proportional to
     A(u). With Ricci zero the curvature is pure Weyl.

What this script does NOT license, recorded because the first draft
claimed it: sqrt(-g) = 1 says nothing about the volume of a causal
diamond, whose boundary moves with the light cones. See design check 2
(p14_interval_volume_constant_a.py), which measures that shift.

Aborts on failure rather than recording it (review C12/C21 lineage).
Runs in a few seconds; requires sympy (not a package dependency -- this
is a design check, run by hand, and its frozen copy moves to
docs/prereg/frozen/p14/ with a stamp if and when P14 freezes).
"""

import sympy as sp

u, v, x, y = sp.symbols("u v x y", real=True)
A = sp.Function("A")(u)
coords = [u, v, x, y]
n = 4

H = A * (x ** 2 - y ** 2)
g = sp.Matrix([
    [H, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
])
ginv = g.inv()

det = sp.simplify(g.det())
assert det == -1, f"det g = {det}, expected -1"

Gamma = [[[0] * n for _ in range(n)] for _ in range(n)]
for a in range(n):
    for b in range(n):
        for c in range(n):
            s = sum(
                ginv[a, d] * (sp.diff(g[d, b], coords[c])
                              + sp.diff(g[d, c], coords[b])
                              - sp.diff(g[b, c], coords[d]))
                for d in range(n))
            Gamma[a][b][c] = sp.simplify(s / 2)

Riem = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
for a in range(n):
    for b in range(n):
        for c in range(n):
            for d in range(n):
                t = sp.diff(Gamma[a][b][d], coords[c]) \
                    - sp.diff(Gamma[a][b][c], coords[d])
                for e in range(n):
                    t += Gamma[a][c][e] * Gamma[e][b][d] \
                         - Gamma[a][d][e] * Gamma[e][b][c]
                Riem[a][b][c][d] = sp.simplify(t)

Ric = sp.zeros(n)
for b in range(n):
    for d in range(n):
        Ric[b, d] = sp.simplify(sum(Riem[a][b][a][d] for a in range(n)))

assert sp.simplify(Ric) == sp.zeros(n), f"Ricci != 0:\n{Ric}"

nonzero = [(a, b, c, d)
           for a in range(n) for b in range(n)
           for c in range(n) for d in range(n)
           if Riem[a][b][c][d] != 0]
assert len(nonzero) == 8, f"expected 8 Riemann components, got {len(nonzero)}"
for a, b, c, d in nonzero:
    ratio = sp.simplify(Riem[a][b][c][d] / A)
    assert ratio in (1, -1), (
        f"component {(a, b, c, d)} is {Riem[a][b][c][d]}, not +/-A(u)")

print("det g              = -1                     [exact]")
print("sqrt(-g)           =  1                     [exact]")
print("Ricci              =  0 (all 16 components) [exact]")
print(f"Riemann != 0       : {len(nonzero)} components, all +/-A(u)")
print("=> vacuum, curvature pure Weyl, volume form exactly flat: PASS")
