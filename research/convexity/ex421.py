"""Reproduce and probe thesis Ex. 4.2.1 (F+lambda 8pt instance, eq. 32)."""
import sys, os, itertools, random
from math import comb
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
F1, F2 = pb.supports
Q = Polytope(F1) + Polytope(F2)
B = Q.lattice_points()
print('|B| =', len(B), ' (thesis: 30)')
Bs = set(B)
def translates(A, Bs):
    return {(b[0]-A[0][0], b[1]-A[0][1]) for b in Bs if all((b[0]-A[0][0]+a[0], b[1]-A[0][1]+a[1]) in Bs for a in A)}
print('translates of F1, F2 inside B:', len(translates(F1, Bs)), len(translates(F2, Bs)), ' (thesis: 6, 20)')
shiftable = [b for b in B if (b[0], b[1]+1) in Bs]
print('points whose x2-shift stays in B:', len(shiftable), ' C(N,8) =', comb(len(shiftable), 8), ' (thesis: 1081575)')
print('C(24,8) =', comb(24, 8), ' (thesis non-convex count 735471)')
# mixed volume / BKK
from fractions import Fraction
def area(P):
    import numpy as np
    from scipy.spatial import ConvexHull
    return ConvexHull(np.array(P.vertices, float)).volume
mv = area(Q) - area(Polytope(F1)) - area(Polytope(F2))
print('mixed volume (BKK bound) =', round(mv, 6), ' true #solutions (thesis) = 8')
