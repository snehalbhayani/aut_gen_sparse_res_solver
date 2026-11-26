"""Exact-count theorem check: corank(kQ) and feasibility for dilations kQ, Q = NP(f1)+NP(f2).

Thesis Ex. 4.2.1 instance: Z = 8 torus solutions + 1 point on the facet w=(0,1) orbit
(x2 = infinity, x1 = 3).  Prediction: corank(kQ) = 9 for k >> 0; hidden x1 feasible
(x1 regular on Z) with eigen size 9; hidden x2 never feasible (x2 has a pole on Z).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import sympy
from symgen.problem import Problem
from symgen.polytope import Polytope
from symgen.template import TemplateBuilder
from symgen import ff

x1, x2 = sympy.symbols('a1 a2')
f1 = (-45*x1**3*x2**3 - 9*x1**3*x2**2 + 375*x1**2*x2**3 + 765*x1**3*x2 + 228*x1**2*x2**2
      - 660*x1*x2**3 + 626*x1**3 - 307*x1**2*x2 - 366*x1*x2**2 - 180*x2**3 + 150*x1**2
      - 1052*x1*x2 + 464*x2**2 - 1096*x1 + 1156*x2 + 156)
f2 = 3*x1*x2**2 + 10*x1*x2 - 9*x2**2 + 3*x1 - 4*x2 - 8
pb = Problem('F_lambda_8pt_instance', [x1, x2], [], [f1, f2])
Q = Polytope(pb.supports[0]) + Polytope(pb.supports[1])
for k in range(1, 7):
    Bk = sorted({tuple(k * v for v in p) for p in Q.vertices})
    B = Polytope(Bk).lattice_points()
    line = f'k={k}  |B|={len(B):4d}'
    for h in (0, 1):
        tb = TemplateBuilder(pb, h)
        top, E = tb.from_basis(B)
        cols = sorted(tb.columns(top, E))
        cork = len(cols) - ff.rank(tb.matrix(top, cols))
        G = sorted(set(cols) - set(E))
        feas = ff.rank(tb.matrix(top, G)) == len(G)
        eig = len(tb.reduce_eigen_basis(top, E)) if feas else None
        line += f' | hidden a{h+1}: corank {cork:2d}, feasible {str(feas):5s}, eigen {eig}'
    print(line, flush=True)
