"""Support bound from Lemma L (sound version).

Every valid template uses a set U of equations on whose common torus zeros x_h takes finitely
many values; C then contains a translate t_f + supp f for each f in U, so

    |C| >= min over translations |U_f (t_f + supp f)| >= max_{f,g in U} min_t |supp f u (t + supp g)|,

and |G| = |C| - |E| at eigen size |E|.  The pair term is |S_f| + |S_g| - (max overlap of S_f
with a translate of S_g), computed exactly from the difference multiset.

We CERTIFY that a subset U is infinite (so it cannot be the used set): at a true solution
v (from symgen's solver, random data) the Jacobian of U has full row rank, so V(U) is
locally smooth of dimension n - |U| > 0, and some tangent vector moves x_h.  Subsets that
are not certified infinite are conservatively treated as finite (this can only weaken the
bound).
"""
import sys, os, itertools
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import sympy
from symgen.problem import load_matlab_problem
from symgen.template import search
from symgen.numeric import NumericSolver

for name, eig, best_g in (('problem_relpose_6p_focal', 18, 12), ('problem_relpose_7p_fr_1s_partial_elim', 19, 22), ('problem_9pt2radial', 27, 90),
                          ('problem_opt_pnp_hesch', 27, 87)):
    pb = load_matlab_problem(f'problems/{name}.m')
    T = search(pb, keep=1, verbose=False, cfg={'try_orders': False})
    h = T.h
    S = NumericSolver(pb, T)
    rng = np.random.default_rng(1)
    d = rng.standard_normal(len(pb.params))
    sols = S.solve(d); res = S.residuals(sols, d)
    v = sols[np.argmin(res)]                         # an accurate solution
    J = sympy.Matrix(pb.eqs).jacobian(pb.vars)
    Jf = sympy.lambdify(pb.vars + pb.params, J, 'numpy')
    Jv = np.array(Jf(*v, *d), dtype=complex)
    m = len(pb.eqs); sizes = [len(s) for s in pb.supports]
    cert_inf = []
    for k in range(1, m + 1):
        for U in itertools.combinations(range(m), k):
            JU = Jv[list(U)]
            sv = np.linalg.svd(JU, compute_uv=False)
            if k < pb.n and sv[-1] > 1e-8 * sv[0]:        # full row rank -> smooth, dim n-k > 0
                _, _, Vh = np.linalg.svd(JU)
                ker = Vh[k:].conj().T                       # tangent directions
                if np.linalg.norm(ker[h]) > 1e-6:           # x_h moves along V(U)
                    cert_inf.append(U)
    finite = [U for k in range(1, m + 1) for U in itertools.combinations(range(m), k) if U not in cert_inf]
    S = [set(map(tuple, sp)) for sp in pb.supports]
    def pair(i, j):
        diffs = Counter(tuple(a - b for a, b in zip(p, q)) for p in S[i] for q in S[j])
        return len(S[i]) + len(S[j]) - max(diffs.values())
    def bound(U):
        return max([sizes[i] for i in U] + [pair(i, j) for i, j in itertools.combinations(U, 2)])
    lbC, U = min((bound(U), U) for U in finite)
    print(f'{name} (hidden a{h+1}): supports {sizes}; certified-infinite subsets: {len(cert_inf)}; '
          f'cheapest possibly-finite subset {tuple(i+1 for i in U)} => |C| >= {lbC}, '
          f'|G| >= {lbC - eig} at eigen {eig}  (symgen |G| = {best_g})', flush=True)
