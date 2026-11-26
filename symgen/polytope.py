"""Exact integer polytope helpers (Newton polytopes, Minkowski sums, lattice points).

Qhull (via scipy) is only used to find *which* points span each facet; the
facet inequalities themselves are recomputed exactly in integers, so the
lattice-point enumeration contains no floating point decisions.
"""
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd

import numpy as np
import sympy
from scipy.spatial import ConvexHull


def _affine_dim(pts):
    pts = np.asarray(pts, dtype=float)
    if len(pts) <= 1:
        return 0
    return int(np.linalg.matrix_rank(pts[1:] - pts[0]))


def _integer_normal(face_pts):
    """Exact primitive integer normal of the hyperplane through n points."""
    v0 = face_pts[0]
    D = sympy.Matrix([[int(x) - int(y) for x, y in zip(v, v0)] for v in face_pts[1:]])
    ns = D.nullspace()  # face_pts may be more than n points; nullspace is still 1-dim
    if len(ns) != 1:
        return None
    vec = ns[0]
    den = reduce(lambda a, b: a * b // gcd(a, b), [Fraction(str(x)).denominator for x in vec], 1)
    a = [int(Fraction(str(x)) * den) for x in vec]
    g = reduce(gcd, [abs(x) for x in a])
    return tuple(x // g for x in a)


class Polytope:
    """Convex hull of a finite set of integer points in Z^n."""

    def __init__(self, points):
        pts = sorted({tuple(int(x) for x in p) for p in points})
        self.n = len(pts[0])
        self.dim = _affine_dim(pts)
        self.facets = None
        if self.dim == self.n and self.n >= 2:
            hull = ConvexHull(np.array(pts, dtype=float))
            self.vertices = [pts[i] for i in sorted(hull.vertices)]
            facets = set()
            P = np.array(pts, dtype=float)
            for eq in hull.equations:
                on = [pts[i] for i in np.flatnonzero(np.abs(P @ eq[:-1] + eq[-1]) < 1e-7)]
                a = _integer_normal(on)
                if a is None:
                    raise ValueError('could not recover an exact facet normal')
                b = sum(x * y for x, y in zip(a, on[0]))
                # orient as a.x <= b
                if any(sum(x * y for x, y in zip(a, p)) > b for p in self.vertices):
                    a, b = tuple(-x for x in a), -b
                facets.add((a, b))
            self.facets = sorted(facets)
        elif self.dim == self.n == 1:
            lo, hi = pts[0][0], pts[-1][0]
            self.vertices = [(lo,), (hi,)] if hi > lo else [(lo,)]
            self.facets = [((1,), hi), ((-1,), -lo)]
        else:
            self.vertices = pts  # lower dimensional: keep all points
        self.lo = [min(p[i] for p in self.vertices) for i in range(self.n)]
        self.hi = [max(p[i] for p in self.vertices) for i in range(self.n)]

    def __add__(self, other):
        return Polytope([tuple(a + b for a, b in zip(p, q)) for p in self.vertices for q in other.vertices])

    def lattice_points(self, shift=None):
        """Integer points p with p - eps*shift in the polytope, for infinitesimal eps > 0.

        This is the Canny-Emiris basis Z^n cap (Q + delta) with delta = -eps*shift
        evaluated exactly: a point on a facet a.x = b is kept iff a.shift > 0
        (or a.shift == 0, i.e. the shift is parallel to that facet).
        Lower dimensional polytopes have no points for a generic shift.
        """
        if self.facets is None:
            return [] if shift is not None and any(shift) else list(self.vertices)
        out = []
        for p in product(*[range(l, h + 1) for l, h in zip(self.lo, self.hi)]):
            ok = True
            for a, b in self.facets:
                s = sum(x * y for x, y in zip(a, p)) - b
                if s > 0 or (s == 0 and shift is not None and sum(x * y for x, y in zip(a, shift)) < 0):
                    ok = False
                    break
            if ok:
                out.append(p)
        return out


def minkowski(polys):
    return reduce(lambda a, b: a + b, polys)
