"""Verify the SET COVER -> MIN-ELIM (flat templates, linear system) reduction by brute force.

Variables: h, x in X, y_S.  Equations (linear parts):
  f0       = h - sum_x x
  f_{S,x}  = y_S + r_{S,x} x          (x in S)
  f'_S     = y_S + r'_S x_1(S)        (r' != r: unbalanced 2-cycle)
Flat template: U = G (variables, h in U); rows = equations supported in U; valid iff rank = |U|.
Claim: min |U| = 1 + |X| + (minimum set cover).
"""
import itertools, random
from fractions import Fraction
import sympy


def rank_Q(rows, cols):
    if not rows:
        return 0
    return sympy.Matrix([[r.get(c, 0) for c in cols] for r in rows]).rank()


def instance(nX, nS, rng, coeffs):
    X = list(range(nX))
    while True:
        sets = [tuple(sorted(rng.sample(X, rng.randint(1, min(3, nX))))) for _ in range(nS)]
        if set().union(*sets) == set(X):
            break
    eqs = [{'h': 1, **{('x', x): -1 for x in X}}]
    for s, S in enumerate(sets):
        for x in S:
            eqs.append({('y', s): 1, ('x', x): coeffs(rng)})
        eqs.append({('y', s): 1, ('x', S[0]): coeffs(rng) + 7})   # different ratio on the double edge
    return X, sets, eqs


def min_flat(X, sets, eqs):
    allvars = ['h'] + [('x', x) for x in X] + [('y', s) for s in range(len(sets))]
    best = None
    others = allvars[1:]
    for size in range(0, len(others) + 1):
        for extra in itertools.combinations(others, size):
            U = set(['h', *extra])
            rows = [e for e in eqs if set(e) <= U]
            if rank_Q(rows, sorted(U, key=str)) == len(U):
                return len(U)
    return None


def min_cover(X, sets):
    for k in range(1, len(sets) + 1):
        for C in itertools.combinations(range(len(sets)), k):
            if set().union(*[sets[i] for i in C]) == set(X):
                return k


if __name__ == '__main__':
    rng = random.Random(0)
    ok = 0
    for trial in range(40):
        nX, nS = rng.randint(2, 4), rng.randint(2, 5)
        regime = trial % 2
        coeffs = (lambda r: Fraction(r.randint(1, 3))) if regime == 0 else (lambda r: Fraction(r.randint(1, 10**6), r.randint(1, 10**6)))
        X, sets, eqs = instance(nX, nS, rng, coeffs)
        m = min_flat(X, sets, eqs)
        c = min_cover(X, sets)
        assert m == 1 + nX + c, (X, sets, m, c)
        ok += 1
    print(f'{ok}/40 random instances: min flat |G| = 1 + |X| + min set cover  (small integer and random rational coefficients)')
