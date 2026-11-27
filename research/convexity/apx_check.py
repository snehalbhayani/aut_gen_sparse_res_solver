"""Sanity check of Theorem C4 (APX-hardness): the C2 gadget system of a cubic graph.

Elements X = edges of a 3-regular graph, one set per vertex (its 3 incident edges), so every
element lies in exactly two sets and min set cover = min vertex cover tau.  Theorem C2 gives
min |G| = 1 + |X| + tau; C4 only adds the arithmetic |X| = 3n/2 <= 3 tau.

Checked here, exactly over GF(p):
  * K4 and the triangular prism: brute force over ALL flat variable sets (h, x's and y's free);
  * the cube Q3 and the Petersen graph: brute force over the set variables (X forced by f0).
"""
import itertools, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from symgen import ff

GRAPHS = {
    'K4': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    'prism': [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    'cube Q3': [(a, a ^ b) for a in range(8) for b in (1, 2, 4) if a < a ^ b],
    'Petersen': [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)]
                + [(5 + i, 5 + (i + 2) % 5) for i in range(5)],
}


def system(edges, rng):
    n = 1 + max(max(e) for e in edges)
    X = list(range(len(edges)))
    sets = [[k for k, e in enumerate(edges) if v in e] for v in range(n)]
    assert all(len(S) == 3 for S in sets), 'not cubic'
    eqs = [{'h': 1, **{('x', x): -1 for x in X}}]
    for v, S in enumerate(sets):
        for x in S:
            eqs.append({('y', v): 1, ('x', x): rng.randint(1, 9)})
        eqs.append({('y', v): 1, ('x', S[0]): rng.randint(11, 19)})   # unbalanced double edge
    return X, sets, eqs


def valid(U, eqs):
    cols = sorted(U, key=str)
    rows = [[e.get(c, 0) % ff.P for c in cols] for e in eqs if set(e) <= U]
    return bool(rows) and ff.rank(rows) == len(U)


def min_flat_all(X, sets, eqs):
    others = [('x', x) for x in X] + [('y', v) for v in range(len(sets))]
    for k in range(len(others) + 1):
        for extra in itertools.combinations(others, k):
            if valid({'h', *extra}, eqs):
                return 1 + k


def min_flat_y(X, sets, eqs):
    base = {'h', *[('x', x) for x in X]}
    for k in range(len(sets) + 1):
        for Y in itertools.combinations(range(len(sets)), k):
            if valid(base | {('y', v) for v in Y}, eqs):
                return len(base) + k


def vertex_cover(edges, n):
    for k in range(n + 1):
        for C in itertools.combinations(range(n), k):
            if all(a in C or b in C for a, b in edges):
                return k


rng = random.Random(0)
for name, edges in GRAPHS.items():
    n = 1 + max(max(e) for e in edges)
    X, sets, eqs = system(edges, rng)
    tau = vertex_cover(edges, n)
    full = len(X) + len(sets) <= 15
    m = (min_flat_all if full else min_flat_y)(X, sets, eqs)
    assert m == 1 + len(X) + tau, (name, m, tau)
    print(f'{name:9s} n={n:2d} |X|={len(X):2d} tau={tau}: min |G| = {m} = 1 + |X| + tau '
          f'({"all variable sets" if full else "set variables"}); |X| <= 3 tau: {len(X) <= 3 * tau}')
