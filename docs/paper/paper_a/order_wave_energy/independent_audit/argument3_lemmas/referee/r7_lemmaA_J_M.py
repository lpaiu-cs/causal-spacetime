# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r7: Lemma A (index count, moments, Hoeffding representation), Lemma J, Lemma M count.

Exact enumeration only (Fractions), small n.
 (a) #{tuples in [n]^(2r) with every value occurring >= 2 times} <= B_(2r) n^r.
 (b) Lone-index count for degenerate kernels: #{(j1<k1, ..., j_r<k_r) pairs-tuples with every
     index >= 2 occurrences} <= B_(2r... ) n^(r)  (r factors, 2r slots -> <= r distinct).
 (c) E(U - EU)^4 for a two-reference U-statistic on a 3-point space vs the Lemma A route
     2^3 (E(lin)^4 + E(deg)^4), and the Hoeffding-representation identity U = mean_pi W_pi.
 (d) Lemma J: fourth norm of a centred mean of n iid D in {-1,0,1} with P(D!=0)=pi against
     (4 pi)^(1/4) n^(-3/4) + 3^(1/4) pi^(1/2) n^(-1/2); exact.
 (e) Lemma M: number of distinct pilot raw moments needed to certify every labelling of a
     triple (p,f,a per anchor + b for each ORDERED pair).
"""
import itertools
import math
from fractions import Fraction as Fr


def bell(n):
    B = [1]
    row = [1]
    for _ in range(n):
        new = [row[-1]]
        for x in row:
            new.append(new[-1] + x)
        row = new
        B.append(row[0])
    return B[n]


out = []
# (a)
ok = True
for n in range(1, 6):
    for r in range(1, 4):
        cnt = 0
        for t in itertools.product(range(n), repeat=2 * r):
            c = [t.count(v) for v in set(t)]
            if all(x >= 2 for x in c):
                cnt += 1
        ok &= cnt <= bell(2 * r) * n**r
out.append(f"(a) tuple count <= B_(2r) n^r for n<=5, r<=3: {ok}")

# (b) degenerate pairs: r factors h2(j,k) with j<k; survivors have each index >= 2 times
ok = True
for n in range(2, 7):
    pairs = [(j, kk) for j in range(n) for kk in range(j + 1, n)]
    for r in (2, 3):
        cnt = 0
        for tup in itertools.product(pairs, repeat=r):
            flat = [x for p in tup for x in p]
            if all(flat.count(v) >= 2 for v in set(flat)):
                cnt += 1
        ok &= cnt <= bell(2 * r) * n**r
out.append(f"(b) degenerate-pair survivor count <= B_(2r) n^r (r factors, <= r distinct indices): {ok}")

# (c) U-statistic on a 3-point space with weights; exact moments for n = 2..5
vals = [0, 1, 2]
pr = [Fr(1, 5), Fr(1, 2), Fr(3, 10)]
h = {(a, b): Fr(int((a + b) % 3 == 0) + int(a == b), 2) for a in vals for b in vals}  # symmetric, in [0,1]
Eh = sum(pr[a] * pr[b] * h[(a, b)] for a in vals for b in vals)
h1 = {a: sum(pr[b] * h[(a, b)] for b in vals) - Eh for a in vals}
for n in (2, 3, 4, 5):
    m4 = Fr(0)
    lin4 = Fr(0)
    deg4 = Fr(0)
    rep_ok = True
    npairs = n * (n - 1) // 2
    for x in itertools.product(vals, repeat=n):
        p = Fr(1)
        for v in x:
            p *= pr[v]
        U = sum(h[(x[i], x[j])] for i in range(n) for j in range(i + 1, n)) / npairs
        lin = Fr(2, n) * sum(h1[v] for v in x)
        deg = U - Eh - lin
        m4 += p * (U - Eh) ** 4
        lin4 += p * lin**4
        deg4 += p * deg**4
        if n <= 4:
            mp_ = n // 2
            Ws = []
            for perm in itertools.permutations(range(n)):
                Ws.append(sum(h[(x[perm[2 * l]], x[perm[2 * l + 1]])] for l in range(mp_)) / mp_)
            rep_ok &= sum(Ws) / len(Ws) == U
    out.append(f"(c) n={n}: E(U-EU)^4={float(m4):.3e}  route 8(E lin^4 + E deg^4)={float(8*(lin4+deg4)):.3e}"
               f"  holds={m4 <= 8 * (lin4 + deg4)}  E deg^4 * n^4 = {float(deg4 * n**4):.3f}"
               + (f"  Hoeffding representation exact={rep_ok}" if n <= 4 else ""))

# (d) Lemma J
worst = 0.0
for n in (1, 2, 5, 20, 100, 1000):
    for pi in (Fr(1, 1000), Fr(1, 100), Fr(1, 10), Fr(1, 2)):
        # D = +1 w.p. pi/2, -1 w.p. pi/4, 0 else  (generic, mean not zero)
        dist = {1: pi / 2, -1: pi / 4, 0: 1 - 3 * pi / 4}
        mu = sum(v * p for v, p in dist.items())
        c2 = sum((v - mu) ** 2 * p for v, p in dist.items())
        c4 = sum((v - mu) ** 4 * p for v, p in dist.items())
        m4 = n * c4 + 3 * n * (n - 1) * c2**2
        lhs = (float(m4) / n**4) ** 0.25
        rhs = (4 * float(pi)) ** 0.25 * n**-0.75 + 3**0.25 * float(pi) ** 0.5 * n**-0.5
        worst = max(worst, lhs / rhs)
out.append(f"(d) Lemma J fourth-norm bound: max lhs/rhs = {worst:.4f} (<= 1)")

# (e) Lemma M
n_moments_one_labelling = 9 + 3
n_moments_all = 9 + 6
out.append(f"(e) one labelling uses {n_moments_one_labelling} raw moments (union-bound factor 2*12 = 24);"
           f" certifying every labelling of an unordered triple needs {n_moments_all} (factor 30)")
print("\n".join(out))
