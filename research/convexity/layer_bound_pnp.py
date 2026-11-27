"""Theorem P1e for opt PnP (a3 hidden): every valid template has |G| >= 32 (generic data).

Face w = (1, 1, 1): the initial forms are the cubic parts g1, g2, g3, all supported on the full
degree-3 triangle T3 (coordinates (a1, a2), a3 = D - a1 - a2).  The top layer L of any valid
template lies in G and is complete for (g1, g2, g3) (Theorem P1a, Lean: TopLayer).

  (a)-(d)  data-independent necessary conditions (layer_search.py);
  (e)      shifted into N^3 the rows of L are polynomial multiples of the g_i, so every
           monomial of L lies in I = (g1, g2, g3).  If no degree-6 monomial lies in I (checked
           below: an open condition, so one instance certifies it for generic data), then the
           triangle bounding L has side >= 7.

The size-ordered search finds the first set passing (a)-(e) at |L| = 32, and that set is
complete (exact rank 32 at an instance => full rank for generic data).  So lambda = 32.
"""
import itertools, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symgen.problem import load_matlab_problem
from symgen import ff
from face_test import instance_polys, initial_system_w
from layer_search import layer_search_by_size

add = lambda a, b: tuple(x + y for x, y in zip(a, b))
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))
pb = load_matlab_problem('problems/problem_opt_pnp_hesch.m')


def face_system(seed):
    init = initial_system_w(instance_polys(pb, seed=seed, p=ff.P), (1, 1, 1))
    return [{(m[0], m[1]): c for m, c in f.items()} for f in init]


def rows_on(init, L):
    idx = {q: k for k, q in enumerate(sorted(L))}
    rows = []
    for f in init:
        S = list(f)
        for t in {sub(c, S[0]) for c in L}:
            if all(add(t, b) in idx for b in S):
                r = [0] * len(L)
                for m, c in f.items():
                    r[idx[add(t, m)]] = c % ff.P
                rows.append(r)
    return rows, sorted(L)


T3 = sorted(p for p in itertools.product(range(4), repeat=2) if sum(p) <= 3)
assert all(sorted(f) == T3 for f in face_system(0)), 'cubic parts are not dense'

# (e): no degree-6 monomial in the ideal (rank 27 of 28, and e_b outside the row space)
L6 = {p for p in itertools.product(range(7), repeat=2) if sum(p) <= 6}
for seed in (0, 1, 2):
    rows, cols = rows_on(face_system(seed), L6)
    r0 = ff.rank(rows)
    bad = [b for b in cols if ff.rank(rows + [[int(q == b) for q in cols]]) == r0]
    assert r0 == 27 and not bad
print('(e) certified: the degree-6 Macaulay matrix has rank 27/28 and contains no monomial (3 instances)')

side = lambda L: max(a + b for a, b in L) - min(a for a, b in L) - min(b for a, b in L)
size, L = layer_search_by_size([T3, T3, T3], 0, 36, extra_ok=lambda L: side(L) >= 7, verbose=False)
rows, cols = rows_on(face_system(0), set(L))
print(f'smallest layer passing (a)-(e): |L| = {size}; exact rank of its {len(rows)} rows: {ff.rank(rows)}')
print(f'=> lambda_(1,1,1) = {size}, so |G| >= {size} for every valid opt PnP template (symgen: 87)')
