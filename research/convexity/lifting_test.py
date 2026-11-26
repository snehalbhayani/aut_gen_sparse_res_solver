"""Lifting test: can non-flat templates (degree <= 2 columns) beat the flat optimum on
set-cover gadget systems?  Exact GF(p) search with symgen (eigen size free)."""
import sys, os, itertools, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import sympy
from symgen.problem import Problem
from symgen.template import TemplateBuilder
from symgen import ff
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from setcover_reduction import min_cover


def build(X, sets, rng):
    h = sympy.Symbol('h'); xs = {x: sympy.Symbol(f'x{x}') for x in X}; ys = {s: sympy.Symbol(f'y{s}') for s in range(len(sets))}
    vars_ = [h] + [xs[x] for x in X] + [ys[s] for s in range(len(sets))]
    pt = {v: rng.randint(1, 9) for v in vars_}
    lin = [h - sum(xs.values())]
    for s, S in enumerate(sets):
        for x in S:
            lin.append(ys[s] + rng.randint(1, 3) * xs[x])
        lin.append(ys[s] + (rng.randint(1, 3) + 7) * xs[S[0]])
    eqs = [sympy.expand(l - l.subs(pt)) for l in lin]
    return Problem('sc', vars_, [], eqs), vars_


rng = random.Random(5)
for trial in range(6):
    X = list(range(rng.randint(2, 3)))
    while True:
        sets = [tuple(sorted(rng.sample(X, rng.randint(1, len(X))))) for _ in range(rng.randint(2, 4))]
        if set().union(*sets) == set(X):
            break
    pb, vars_ = build(X, sets, rng)
    n = len(vars_)
    flat = 1 + len(X) + min_cover(X, sets)
    best = None
    for deg in (1, 2):
        B = [m for m in itertools.product(range(deg + 1), repeat=n) if sum(m) <= deg]
        for W in [tuple([1] * n)] + [tuple(2 if j == i else 1 for j in range(n)) for i in range(n)]:
            tb = TemplateBuilder(pb, 0, weights=W, eigen_target=10**6)
            rep = tb.repair(*tb.from_basis(B))
            if rep is None:
                continue
            top, E = rep
            for _ in range(3):
                top = tb.reduce_columns(top, E); top = tb.squarify(top, E); top = tb.reduce_rows(top, E); top = tb.squarify(top, E)
            G = sorted(tb.columns(top, E) - set(E))
            assert ff.rank(tb.matrix(top, G)) == len(G)
            if best is None or len(G) < best[0]:
                best = (len(G), deg, len(E))
    print(f'X={X} sets={sets}: flat optimum |G| = {flat}; best found (deg<=2 search) |G| = {best[0]} '
          f'(column degree {best[1]}, eigen size {best[2]})', flush=True)
