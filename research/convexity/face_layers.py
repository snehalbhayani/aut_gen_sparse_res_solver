"""P1a face layers for a benchmark: which initial systems in_w(F), w_h > 0, force big layers.

For each weight w with w_h > 0 the top w-layer L of C lies in G and must be *complete* for
in_w(F) (Theorem P1a).  A complete L has no kernel vector, so the in_w-equations used in L
have no common torus zero.  Hence L contains a translate of supp in_w(f) for some f in every
subset of in_w(F) that has no torus zero:

    lambda_w >= min over zero-free subsets U of max_{f in U} |supp in_w f|.

Combined with the global anchor A (a translate of supp f_a, Lemma L) the top layer either
avoids A or meets it exactly in A's top face, so

    |C| >= |A| + lambda_w - |face_w(supp f_a)|.

Torus zeros of each subset are counted with a Groebner basis over GF(32003) (face_test.py).

    python3 research/convexity/face_layers.py problem_9pt2radial 3 4
"""
import itertools, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from symgen.problem import load_matlab_problem
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from face_test import instance_polys, initial_system_w, unimodular_with_first_row, torus_solution_count


def moved_initial(polys, w, n):
    U = unimodular_with_first_row(w)
    out = []
    for f in initial_system_w(polys, w):
        g = {tuple(sum(U[i][k] * m[k] for k in range(n)) for i in range(n)): c for m, c in f.items()}
        lo = [min(b[i] for b in g) for i in range(n)]
        out.append({tuple(b[i] - lo[i] for i in range(n)): c for b, c in g.items()})
    return out


def primitive(w):
    from math import gcd
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    return g == 1


def main():
    name, h, anchor = sys.argv[1], int(sys.argv[2]) - 1, int(sys.argv[3]) - 1
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    pb = load_matlab_problem(f'problems/{name}.m')
    n, m = pb.n, len(pb.eqs)
    polys = instance_polys(pb)
    res = []
    for w in itertools.product(range(-R, R + 1), repeat=n):
        if w[h] <= 0 or not primitive(w):
            continue
        init = moved_initial(polys, w, n)
        sizes = [len(f) for f in init]
        free = []
        for k in range(1, m + 1):
            for U in itertools.combinations(range(m), k):
                if torus_solution_count([init[i] for i in U], n, 0) == 0:
                    free.append(U)
        if not free:
            print(f'w={w}: in_w(F) has torus zeros -> no template at all', flush=True)
            continue
        lam, U = min((max(sizes[i] for i in U), U) for U in free)
        gain = lam - sizes[anchor]
        res.append((gain, lam, w, U, sizes))
        print(f'w={w}: |supp in_w f| = {sizes}; cheapest zero-free subset {tuple(i+1 for i in U)} '
              f'-> lambda_w >= {lam}; |face_w(anchor)| = {sizes[anchor]}', flush=True)
    best = max(res)
    print(f'best single-layer gain over the anchor: {best[0]} at w = {best[2]} (lambda_w >= {best[1]})')


if __name__ == '__main__':
    main()
