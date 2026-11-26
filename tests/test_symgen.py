"""Tests for symgen.  Run:  python3 -m pytest tests  (or python3 tests/test_symgen.py)"""
import os
import random
import shutil
import sys
import unittest

import numpy as np
import sympy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from symgen import ff  # noqa: E402
from symgen.codegen import matlab_solver  # noqa: E402
from symgen.numeric import evaluate  # noqa: E402
from symgen.polytope import Polytope  # noqa: E402
from symgen.problem import Problem, load_matlab_problem  # noqa: E402
from symgen.template import TemplateBuilder, search  # noqa: E402


class TestFiniteField(unittest.TestCase):
    def test_rank_matches_sympy(self):
        rng = random.Random(3)
        for _ in range(20):
            m, n, r = rng.randint(1, 8), rng.randint(1, 8), rng.randint(0, 5)
            L = sympy.Matrix(m, r, lambda i, j: rng.randint(-5, 5))
            R = sympy.Matrix(r, n, lambda i, j: rng.randint(-5, 5))
            A = L * R if r else sympy.zeros(m, n)
            Ap = np.array([[int(x) % ff.P for x in A.row(i)] for i in range(m)], dtype=float)
            self.assertEqual(ff.rank(Ap), A.rank())

    def test_incremental_span_is_exact_for_long_vectors(self):
        rng = np.random.default_rng(0)
        n = 3000  # longer than one chunk: exercises the overflow-safe dot product
        S = ff.IncrementalSpan(n)
        vs = rng.integers(0, ff.P, size=(5, n)).astype(float)
        for v in vs:
            self.assertTrue(S.add(v))
        combo = np.mod(vs[0] * 7 + vs[3] * (ff.P - 2), ff.P)
        self.assertTrue(S.contains(combo))


class TestPolytope(unittest.TestCase):
    def test_square_shifted_lattice_points(self):
        sq = Polytope([(0, 0), (2, 0), (0, 2), (2, 2)])
        self.assertEqual(len(sq.lattice_points()), 9)
        # p - eps*(1,1) in [0,2]^2  <=>  p in {1,2}^2
        self.assertEqual(sorted(sq.lattice_points((1, 1))), [(1, 1), (1, 2), (2, 1), (2, 2)])

    def test_minkowski_sum(self):
        tri = Polytope([(0, 0), (1, 0), (0, 1)])
        self.assertEqual(len((tri + tri).lattice_points()), 6)


def conics_problem():
    """Two generic conics in the plane: 4 solutions."""
    x, y = sympy.symbols('a1 a2')
    c = sympy.symbols('c1:13')
    mons = [x**2, x*y, y**2, x, y, 1]
    f1 = sum(ci * m for ci, m in zip(c[:6], mons))
    f2 = sum(ci * m for ci, m in zip(c[6:], mons))
    return Problem('conics', [x, y], list(c), [f1, f2], {'hidden_var': 1, 'size_of_combs': [2, 3]})


class TestTemplate(unittest.TestCase):
    def test_conics(self):
        pb = conics_problem()
        T = search(pb, verbose=False)
        # Bezout: 4 solutions, and the eigen basis reduction must reach it
        self.assertEqual(T.size[0], 4)
        stats = evaluate(pb, T, trials=20)
        self.assertGreater(stats['good_solutions_per_instance'], 3.9)
        self.assertLess(stats['median_log10_residual'], -10)

    def test_deterministic(self):
        pb = conics_problem()
        a, b = search(pb, verbose=False), search(pb, verbose=False)
        self.assertEqual((a.top, a.E, a.G), (b.top, b.E, b.G))

    def test_eigen_reduction_reaches_corank(self):
        pb = conics_problem()
        tb = TemplateBuilder(pb, 0)
        Q = Polytope(pb.supports[0]) + Polytope(pb.supports[1]) + Polytope([(0, 0), (1, 0), (0, 1)])
        top, E = tb.repair(*tb.from_basis(Q.lattice_points((1, 1))))
        cols = sorted(tb.columns(top, E))
        corank = len(cols) - ff.rank(tb.matrix(top, cols))
        self.assertEqual(len(tb.reduce_eigen_basis(top, E)), corank)

    def test_relpose_6p_focal_regression(self):
        pb = load_matlab_problem(os.path.join(ROOT, 'problems', 'problem_relpose_6p_focal.m'))
        T = search(pb, verbose=False, keep=1, cfg={'try_orders': False})
        self.assertEqual(T.size, (18, 12))  # same as the Maple generated solver (30 columns)
        stats = evaluate(pb, T, trials=10)
        self.assertGreater(stats['good_solutions_per_instance'], 14)


class TestSaturatedFocal(unittest.TestCase):
    def test_parsed_file_matches_derivation(self):
        from symgen.derive_6p_focal_sat import derive
        vars_, params, orig, new = derive()
        pb = load_matlab_problem(os.path.join(ROOT, 'problems', 'problem_relpose_6p_focal_sat.m'))
        self.assertEqual(len(pb.eqs), len(orig) + len(new))
        for e, f in zip(pb.eqs, orig + new):
            self.assertEqual(sympy.expand(e - f), 0)
        self.assertTrue(all(sympy.degree(g, vars_[2]) == 1 for g in new))

    def test_spurious_solutions_removed(self):
        pb = load_matlab_problem(os.path.join(ROOT, 'problems', 'problem_relpose_6p_focal_sat.m'))
        T = search(pb, verbose=False, keep=1, cfg={'combs': [[10, 11, 19]], 'try_orders': False})
        self.assertEqual(T.size[0], 15)  # was 18 = 15 + 3 at a3 = infinity
        stats = evaluate(pb, T, trials=10)
        self.assertGreater(stats['good_solutions_per_instance'], 14)


class TestEigenTarget(unittest.TestCase):
    def test_tradeoff_keeps_validity(self):
        pb = conics_problem()
        small = search(pb, verbose=False)
        big = search(pb, verbose=False, cfg={'eigen_target': small.size[0] + 2})
        self.assertEqual(big.size[0], small.size[0] + 2)
        self.assertLessEqual(big.size[1], small.size[1])
        self.assertGreater(evaluate(pb, big, trials=10)['good_solutions_per_instance'], 3.9)


@unittest.skipUnless(shutil.which('octave'), 'octave not installed')
class TestMatlabCodegen(unittest.TestCase):
    def test_conics_in_octave(self):
        import tempfile
        from symgen.bench import octave_accuracy
        pb = conics_problem()
        T = search(pb, verbose=False)
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'solver.m'), 'w') as fh:
                fh.write(matlab_solver(pb, T))
            r = octave_accuracy(pb, d, T.h + 1, trials=10)
        self.assertGreater(r['good_solutions_per_instance'], 3.9)


if __name__ == '__main__':
    unittest.main()
