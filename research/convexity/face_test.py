"""Bernstein face test: predict spurious eigenvalues / infeasible hidden variables.

For a weight w (here the coordinate directions +-e_j), the initial form in_w(f)
keeps the terms of f on the face of NP(f) maximising w.alpha.  By Bernstein's
second theorem, F has (generically) exactly MV solutions in the torus iff no
initial system in_w(F), w != 0, has a torus solution.  Solutions of in_w(F)
are solutions "at infinity" of the toric compactification.  For the coordinate
direction w = +e_h they are points with x_h = infinity, which can never be
eigenvalues of the finite action matrix X (hidden x_h)  =>  A12 singular.

For w = +-e_j, in_w(f) = x_j^d * g(x without x_j); we count torus solutions of
the g's over GF(p) at the same seeded instance used by symgen (Groebner basis
with a Rabinowitsch variable for the torus).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import sympy
from symgen import ff

P_SMALL = 32003   # sympy Groebner over GF(p) is fast for small p


def instance_polys(problem, seed=0, p=P_SMALL):
    vals = problem.instance_mod_p(seed=seed, p=p)
    polys = []
    for k, ids in enumerate(problem.coeff_ids):
        polys.append({m: vals[c] for m, c in ids.items() if vals[c] % p})
    return polys


def initial_system(polys, j, sign):
    out = []
    for f in polys:
        if not f:
            continue
        key = max(sign * m[j] for m in f)
        out.append({m: c for m, c in f.items() if sign * m[j] == key})
    return out


def torus_solution_count(polys, n, skip, p=P_SMALL):
    xs = sympy.symbols(f'y0:{n}')
    t = sympy.Symbol('t_rab')
    keep = [i for i in range(n) if i != skip]
    gens = [xs[i] for i in keep]
    eqs = []
    for f in polys:
        e = 0
        for m, c in f.items():
            mon = 1
            for i in keep:
                mon *= xs[i] ** m[i]
            e += c * mon
        e = sympy.expand(e)
        if e != 0:
            # divide out the monomial content (irrelevant in the torus)
            P = sympy.Poly(e, *gens, modulus=p)
            mins = [min(mm[i] for mm in P.monoms()) for i in range(len(gens))]
            e = sympy.expand(e / sympy.prod([g ** k for g, k in zip(gens, mins)]))
            if e.free_symbols:
                eqs.append(e)
            elif e % p:
                return 0  # nonzero constant: inconsistent
    if not eqs:
        return 'positive-dimensional'
    eqs.append(1 - t * sympy.prod(gens))
    G = sympy.groebner(eqs, *gens, t, modulus=p, order='grevlex')
    if list(G.exprs) == [1]:
        return 0
    if not G.is_zero_dimensional:
        return 'positive-dimensional'
    # count standard monomials (= #solutions with multiplicity)
    lead = [sympy.Poly(g, *gens, t, modulus=p).monoms(order='grevlex')[0] for g in G.exprs]
    nv = len(gens) + 1
    count, stack, seen = 0, [tuple([0] * nv)], set()
    while stack:
        m = stack.pop()
        if m in seen:
            continue
        seen.add(m)
        if any(all(m[i] >= l[i] for i in range(nv)) for l in lead):
            continue
        count += 1
        for i in range(nv):
            q = list(m); q[i] += 1; stack.append(tuple(q))
    return count


def face_report(problem, seed=0, p=None):
    p = p or P_SMALL
    polys = instance_polys(problem, seed, p)
    rep = {}
    for j in range(problem.n):
        for sign in (1, -1):
            init = initial_system(polys, j, sign)
            rep[('+' if sign > 0 else '-') + f'a{j + 1}'] = torus_solution_count(init, problem.n, j, p)
    return rep


if __name__ == '__main__':
    from symgen.problem import load_matlab_problem
    import glob
    for f in sys.argv[1:] or sorted(glob.glob('problems/problem_[o9r]*.m')):
        pb = load_matlab_problem(f)
        print(pb.name, face_report(pb), flush=True)


# ---------------------------------------------------------------- general weights
def unimodular_with_first_row(w):
    """Integer matrix U, det +-1, whose first row is the primitive vector w."""
    import sympy
    n = len(w)
    V = sympy.eye(n)
    v = sympy.Matrix(w)
    # integer row operations reducing v to (+-1, 0, ..., 0); apply the same to V
    while True:
        nz = [i for i in range(n) if v[i] != 0]
        if len(nz) == 1 and nz[0] == 0:
            break
        piv = min(nz, key=lambda i: abs(v[i]))
        if piv != 0:
            v.row_swap(0, piv); V.row_swap(0, piv)
        for i in range(1, n):
            if v[i] != 0:
                q = v[i] // v[0]
                v[i] -= q * v[0]; V[i, :] = V[i, :] - q * V[0, :]
    if v[0] < 0:
        V[0, :] = -V[0, :]
    U = (V.inv()).T
    assert list(U[0, :]) == list(w) and abs(U.det()) == 1
    return [[int(x) for x in U.row(i)] for i in range(n)]


def initial_system_w(polys, w):
    out = []
    for f in polys:
        if not f:
            continue
        key = max(sum(a * b for a, b in zip(w, m)) for m in f)
        out.append({m: c for m, c in f.items() if sum(a * b for a, b in zip(w, m)) == key})
    return out


def face_count_w(polys, w, n, p=None):
    p = p or P_SMALL
    U = unimodular_with_first_row(w)
    init = initial_system_w(polys, w)
    moved = []
    for f in init:
        g = {}
        for m, c in f.items():
            beta = tuple(sum(U[i][k] * m[k] for k in range(n)) for i in range(n))
            g[beta] = c
        lo = [min(b[i] for b in g) for i in range(n)]
        moved.append({tuple(b[i] - lo[i] for i in range(n)): c for b, c in g.items()})
    return torus_solution_count(moved, n, 0, p)


def facet_report(problem, Q, seed=0, p=None):
    p = p or P_SMALL
    polys = instance_polys(problem, seed, p)
    out = {}
    for a, b in Q.facets:
        cnt = face_count_w(polys, a, problem.n, p)
        if cnt:
            out[a] = cnt
    return out
