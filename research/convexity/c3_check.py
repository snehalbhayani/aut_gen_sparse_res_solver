"""Sanity check of Theorem C3 (N > 1): f0' = prod_i (h - a_i) - (sum x - sum x*), gadgets as in C2.
Template: E = {1, h, ..., h^(N-1)}, G = {h^N} + X + Y(cover).  Check exact validity, and that
dropping a cover set (or using a non-cover) makes it invalid; also solve numerically."""
import itertools, random
import sympy
import numpy as np

rng = random.Random(3)
for trial in range(5):
    N = rng.randint(2, 4)
    X = list(range(rng.randint(2, 3)))
    while True:
        sets = [tuple(sorted(rng.sample(X, rng.randint(1, len(X))))) for _ in range(rng.randint(2, 4))]
        if set().union(*sets) == set(X):
            break
    h = sympy.Symbol('h'); xs = sympy.symbols(f'x0:{len(X)}'); ys = sympy.symbols(f'y0:{len(sets)}')
    xstar = {x: rng.randint(1, 5) for x in xs}; ystar = {y: rng.randint(1, 5) for y in ys}
    roots = rng.sample(range(1, 20), N)
    f0 = sympy.expand(sympy.prod([h - a for a in roots]) - (sum(xs) - sum(xstar.values())))
    gad = []
    for s, S in enumerate(sets):
        for i, x in enumerate(S):
            r = rng.randint(1, 3); gad.append((s, ys[s] + r * xs[x] - (ystar[ys[s]] + r * xstar[xs[x]])))
        r2 = rng.randint(4, 9); x1 = xs[S[0]]
        gad.append((s, ys[s] + r2 * x1 - (ystar[ys[s]] + r2 * xstar[x1])))

    def valid(chosen):
        rows = [f0] + [g for s, g in gad if s in chosen]
        G = [h ** N] + list(xs) + [ys[s] for s in chosen]
        M = sympy.Matrix([[sympy.Poly(rw, h, *xs, *ys).coeff_monomial(m) for m in G] for rw in rows])
        return M.rank() == len(G)

    covers = [c for k in range(1, len(sets) + 1) for c in itertools.combinations(range(len(sets)), k)
              if set().union(*[sets[i] for i in c]) == set(X)]
    noncovers = [c for k in range(0, len(sets)) for c in itertools.combinations(range(len(sets)), k)
                 if set().union(*[sets[i] for i in c]) != set(X)]
    assert all(valid(set(c)) for c in covers), 'a cover template is invalid'
    assert not any(valid(set(c)) for c in noncovers), 'a non-cover template is valid'
    best = min(len(c) for c in covers)
    # numeric solve with the minimal template: action matrix on E = {1,h,...,h^(N-1)}
    sol = sympy.solve([f0] + [g for _, g in gad], [h, *xs, *ys], dict=True)
    hs = sorted(float(sv[h]) for sv in sol)
    print(f'N={N} X={X} sets={sets}: covers valid, non-covers invalid; min |G| = {1 + len(X) + best}; '
          f'{len(sol)} solutions, h = {hs} (roots {sorted(roots)})')
