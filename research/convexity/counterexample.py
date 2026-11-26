"""Minkowski sums are not necessary: a valid template from a non-convex column set.

Thesis Ex. 4.2.1 instance (F + lambda 8pt), hidden x1.  Search random boundary swaps of
the Minkowski-sum set B (30 points) for a non-convex set B' whose template is valid
(exact GF(p) rank test: full column rank of the eliminated block certifies validity).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy
from symgen.problem import Problem
from symgen.polytope import Polytope
from symgen.template import TemplateBuilder
from symgen import ff
from exp_shapes import random_swap, convexity_defect

x1, x2 = sympy.symbols('a1 a2')
f1 = (-45*x1**3*x2**3 - 9*x1**3*x2**2 + 375*x1**2*x2**3 + 765*x1**3*x2 + 228*x1**2*x2**2
      - 660*x1*x2**3 + 626*x1**3 - 307*x1**2*x2 - 366*x1*x2**2 - 180*x2**3 + 150*x1**2
      - 1052*x1*x2 + 464*x2**2 - 1096*x1 + 1156*x2 + 156)
f2 = 3*x1*x2**2 + 10*x1*x2 - 9*x2**2 + 3*x1 - 4*x2 - 8
pb = Problem('F_lambda_8pt_instance', [x1, x2], [], [f1, f2])
Q = Polytope(pb.supports[0]) + Polytope(pb.supports[1])
B0 = set(Q.lattice_points())
tb = TemplateBuilder(pb, 0)
rng = random.Random(11)
for trial in range(2000):
    B = random_swap(B0, rng, 4)
    d = convexity_defect(B)
    if not d:
        continue
    top, E = tb.from_basis(sorted(B))
    cols = tb.columns(top, E); G = sorted(cols - set(E))
    if top and E and ff.rank(tb.matrix(top, G)) == len(G) and tb.recoverable(E, cols):
        E1 = tb.reduce_eigen_basis(top, E)
        print(f'found after {trial+1} tries: |B| = {len(B)}, convexity defect = {d} '
              f'(B misses {d} lattice points of its own convex hull)')
        print(f'rows {len(top)}, eliminated columns {len(G)} with full column rank -> VALID, eigen size {len(E1)}')
        for y in range(max(p[1] for p in B), -1, -1):
            print('   ' + ''.join('#' if (x, y) in B else ('o' if (x, y) in B0 else '.') for x in range(0, 7)))
        print("   ('#' = monomial in B', 'o' = Minkowski point left out, '.' = outside)")
        break
