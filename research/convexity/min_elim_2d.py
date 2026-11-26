"""Part 1 probe: smallest elimination block |G| for the thesis Ex. 4.2.1 instance (x1 hidden).

Face-certificate lower bound: the right-edge leading forms have widths 3 and 2, so every
valid template has |G| >= 5.  Search for templates with small |G| (eigen size free):
exact greedy row/column removal from several starting sets, deterministic orders.
"""
import sys, os, itertools
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
best = None
for W in [(1, 1), (2, 1), (1, 2), (1, 3), (3, 1)]:
    tb = TemplateBuilder(pb, 0, weights=W, eigen_target=10**6)   # keep the eigen basis large
    for kx, ky in itertools.product(range(3, 8), range(3, 9)):     # rectangles (incl. ones smaller than the Minkowski sum)
        B = [(a, b) for a in range(kx + 1) for b in range(ky + 1)]
        top, E = tb.from_basis(B)
        rep = tb.repair(top, E)
        if rep is None:
            continue
        top, E = rep
        top = tb.reduce_columns(top, E)
        top = tb.squarify(top, E)
        top = tb.reduce_rows(top, E)
        top = tb.squarify(top, E)
        cols = tb.columns(top, E)
        G = sorted(cols - set(E))
        assert ff.rank(tb.matrix(top, G)) == len(G)
        if best is None or len(G) < best[0]:
            best = (len(G), len(E), W, (kx, ky), G)
print(f'smallest |G| found: {best[0]} (eigen size {best[1]}), start box {best[3]}, order weights {best[2]}')
print('eliminated columns:', best[4])
print('face-certificate lower bound: 5')
