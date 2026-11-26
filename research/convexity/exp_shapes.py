"""Does convexity of B matter for feasibility / eigen size? (thesis Sec. 4.2.1)

For a fixed problem and hidden direction, sample monomial sets B of a fixed size s:
  * the Minkowski-sum set (thesis choice),
  * other convex lattice sets (B = conv(B) cap Z^n) of size s,
  * non-convex sets obtained by random boundary swaps from the convex one,
and record, with exact GF(p) rank tests:
  C1  rank Top = |B|-corank is a valid template at all (rank Top[:,G] = |G| after repair-free check)
  eig minimal eigen size (= |cols| - rank Top, matroid step) if feasible
  rows total number of rows (translates of all supports)
"""
import sys, os, random, json
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.spatial import ConvexHull
from symgen.template import TemplateBuilder
from symgen.polytope import Polytope
from symgen import ff


def convexity_defect(B):
    try:
        P = Polytope(B)
    except Exception:
        return None
    if P.facets is None:
        return None
    return len(P.lattice_points()) - len(B)


def evaluate(tb, B):
    top, E = tb.from_basis(sorted(B))
    if not top or not E:
        return {'feasible': False, 'rows': len(top)}
    cols = tb.columns(top, E)
    G = sorted(cols - set(E))
    ok = len(top) >= len(G) and ff.rank(tb.matrix(top, G)) == len(G) and tb.recoverable(E, cols)
    out = {'feasible': bool(ok), 'rows': len(top), 'cols': len(cols)}
    if ok:
        E1 = tb.reduce_eigen_basis(top, E)
        out['eig'] = len(E1)
        out['corank'] = len(cols) - ff.rank(tb.matrix(top, sorted(cols)))
    return out


def neighbours(p):
    n = len(p)
    for i in range(n):
        for d in (-1, 1):
            q = list(p); q[i] += d
            if min(q) >= 0:
                yield tuple(q)


def random_swap(B, rng, steps):
    """Non-convex perturbation: move boundary points to outside neighbours (keeps size)."""
    B = set(B)
    for _ in range(steps):
        boundary = [p for p in B if any(q not in B for q in neighbours(p))]
        outside = sorted({q for p in B for q in neighbours(p) if q not in B})
        p = rng.choice(sorted(boundary)); q = rng.choice(outside)
        B.remove(p); B.add(q)
    return B


def random_convex(n_target, dim, rng, box):
    """Random convex lattice set: lattice points of the hull of random points, trimmed to size."""
    for _ in range(2000):
        pts = [tuple(rng.randint(0, box) for _ in range(dim)) for _ in range(rng.randint(dim + 1, 2 * dim + 4))]
        try:
            P = Polytope(pts)
        except Exception:
            continue
        if P.facets is None:
            continue
        L = P.lattice_points()
        if len(L) == n_target:
            return set(L)
    return None


def run(problem, h, Bconv, trials=200, seed=0, label=''):
    rng = random.Random(seed)
    tb = TemplateBuilder(problem, h)
    res = {'minkowski': evaluate(tb, Bconv)}
    s = len(Bconv)
    stats = {}
    for kind in ('convex', 'swap2', 'swap5', 'swap10'):
        rows = []
        for t in range(trials):
            if kind == 'convex':
                B = random_convex(s, problem.n, rng, box=max(max(p) for p in Bconv) + 2)
                if B is None:
                    continue
            else:
                B = random_swap(Bconv, rng, int(kind[4:]))
            r = evaluate(tb, B)
            r['defect'] = convexity_defect(B)
            rows.append(r)
        feas = [r for r in rows if r['feasible']]
        stats[kind] = {
            'n': len(rows), 'feasible_rate': round(len(feas) / max(1, len(rows)), 3),
            'mean_rows': round(float(np.mean([r['rows'] for r in rows])), 1) if rows else None,
            'eig_hist': dict(sorted(Counter(r['eig'] for r in feas).items())),
            'mean_defect': round(float(np.mean([r['defect'] for r in rows if r['defect'] is not None])), 2) if rows else None,
        }
    res['samples'] = stats
    print(label, json.dumps(res, indent=1))
    return res


if __name__ == '__main__':
    import sympy
    from symgen.problem import Problem
    x1, x2 = sympy.symbols('a1 a2')
    f1 = (-45*x1**3*x2**3 - 9*x1**3*x2**2 + 375*x1**2*x2**3 + 765*x1**3*x2 + 228*x1**2*x2**2
          - 660*x1*x2**3 + 626*x1**3 - 307*x1**2*x2 - 366*x1*x2**2 - 180*x2**3 + 150*x1**2
          - 1052*x1*x2 + 464*x2**2 - 1096*x1 + 1156*x2 + 156)
    f2 = 3*x1*x2**2 + 10*x1*x2 - 9*x2**2 + 3*x1 - 4*x2 - 8
    pb = Problem('F_lambda_8pt_instance', [x1, x2], [], [f1, f2])
    Q = Polytope(pb.supports[0]) + Polytope(pb.supports[1])
    for h in (0, 1):
        run(pb, h, set(Q.lattice_points()), trials=300, label=f'F+lambda 8pt (thesis Ex. 4.2.1), hidden a{h+1}:')
