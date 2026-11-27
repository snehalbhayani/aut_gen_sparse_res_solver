"""Exhaustive lifting test: every column set C of degree-<=2 monomials with |C| <= CMAX.

For a column set C the optimal eigen basis is maximal, E = {e in C : e + e_h in C},
because enlarging G only makes validity harder.  So C admits a valid template with
|G| = |C \\ E| iff the boundary columns are independent in the matrix of ALL rows
fitting in C.  We search for any C with |G| < flat optimum.
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from symgen import ff

# variables: 0=h, 1=x0, 2=x1, 3=y0, 4=y1 ; X={0,1}, sets {0},{1} (double edges y_s - x_s)
n = 5
def u(i): return tuple(1 if j == i else 0 for j in range(n))
Z = tuple([0] * n)
pt = [3, 5, 7, 2, 4]                                  # the unique solution
lin = [{u(0): 1, u(1): -1, u(2): -1},                 # f0 = h - x0 - x1
       {u(3): 1, u(1): 1}, {u(3): 1, u(1): 2},        # set {0}: y0 + x0, y0 + 2 x0
       {u(4): 1, u(2): 1}, {u(4): 1, u(2): 3}]        # set {1}: y1 + x1, y1 + 3 x1
eqs = []
for f in lin:
    g = dict(f); g[Z] = -sum(c * pt[[i for i in range(n) if m[i]][0]] for m, c in f.items())
    eqs.append({m: c % ff.P for m, c in g.items()})
mons = [m for m in itertools.product(range(3), repeat=n) if sum(m) <= 2]
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))
eh = u(0)
FLAT = int(os.environ.get("FLAT", 5))
CMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9
t0 = time.time(); checked = 0; found = []
for size in range(2, CMAX + 1):
    for C in itertools.combinations(mons, size):
        Cs = set(C)
        bd = [c for c in C if add(c, eh) not in Cs]           # G for the maximal E
        if not bd or len(bd) >= FLAT or len(bd) == size:      # need |G| < flat and E nonempty
            continue
        rows = []
        for f in eqs:
            base = next(iter(f))
            for c in C:
                t = sub(c, base)
                if min(t) >= 0 and all(add(t, m) in Cs for m in f):
                    rows.append({add(t, m): v for m, v in f.items()})
        rows = [dict(r) for r in {tuple(sorted(r.items())) for r in rows}]
        checked += 1
        if len(rows) < len(bd):
            continue
        Mx = np.array([[r.get(b, 0) for b in bd] for r in rows], dtype=float)
        if ff.rank(Mx) == len(bd):
            found.append((len(bd), C))
            print('VALID with |G| =', len(bd), 'C =', C, flush=True)
            break
    print(f'|C| = {size}: done, {checked} column sets checked so far, {time.time() - t0:.0f}s', flush=True)
    if found:
        break
print('result:', 'counterexample to lifting found' if found else f'no template with |G| < {FLAT} among degree-<=2 column sets of size <= {CMAX}')
