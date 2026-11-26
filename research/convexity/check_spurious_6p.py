"""Numerical check of the face theorem (w_h = 0 case) on 6pt focal relative pose.

Face w = (0,0,1): solutions at a3 = infinity satisfy F33(a1,a2) = 0 and det F(a1,a2) = 0
(3 points).  Theorem 2 predicts that, with a2 hidden, the action matrix X has these
points' a2-values as (spurious) eigenvalues.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import sympy
from symgen.problem import load_matlab_problem
from symgen.template import search
from symgen.numeric import NumericSolver

pb = load_matlab_problem('problems/problem_relpose_6p_focal.m')
T = search(pb, keep=1, verbose=False, cfg={'try_orders': False})   # hidden a2, eigen 18
S = NumericSolver(pb, T)
a1, a2 = sympy.symbols('a1 a2')
c = sympy.symbols('c1:28')
M = lambda k: sympy.Matrix(3, 3, list(c[9 * k: 9 * k + 9]))
F = M(0) + a1 * M(1) + a2 * M(2)
a1_of_a2 = sympy.solve(F[2, 2], a1)[0]
cubic = sympy.lambdify(c, sympy.Poly(sympy.numer(sympy.together(F.det().subs(a1, a1_of_a2))), a2).all_coeffs())
rng = np.random.default_rng(3)
worst = 0
for trial in range(20):
    d = rng.standard_normal(27)
    sols = S.solve(d)
    res = S.residuals(sols, d)
    spurious = np.sort_complex(sols[np.argsort(res)[-3:], T.h])        # 3 worst-residual eigenpairs
    predicted = np.sort_complex(np.roots(np.array(cubic(*d), dtype=float)).astype(complex))
    err = max(min(abs(s - p) for p in predicted) / max(1, abs(s)) for s in spurious)
    worst = max(worst, err)
    if trial < 3:
        print('eigenvalues of X with large residual:', np.round(spurious, 6))
        print('a2 of {F33=0, det F=0}             :', np.round(predicted, 6))
print('max relative mismatch over 20 instances: %.2e' % worst)
