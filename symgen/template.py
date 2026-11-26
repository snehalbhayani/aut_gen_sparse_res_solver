"""Sparse-resultant elimination templates with exact, deterministic reduction.

Setting (same as the Maple generator).  Unknowns x = (x_1..x_n), equations
f_1..f_m, hidden variable x_h, and the extra linear polynomial

        l(x, lambda) = x_h - lambda .

A template is
    * top rows   : monomial multiples  t * f_k             (lambda free)
    * l-rows     : monomial multiples  e * l   for e in E   (E = eigen basis)
and its columns are all monomials touched by the rows.  Write the columns as
E (the l-row multipliers) and G (everything else).  For a solution xi with
lambda = xi_h, the monomial vector v(xi) lies in the kernel of all rows:

        top:    A1 v_E + A2 v_G = 0
        l-rows: lambda v_e = v_{e + e_h}          (e in E)

so as soon as A2 = Top[:, G] has full column rank, v_G = -A2^+ A1 v_E and
v_E is an eigenvector of the |E| x |E| action matrix X with eigenvalue xi_h.

    VALID TEMPLATE  <=>  rank Top[:, G] = |G|   (plus variable recovery)

This single exact condition replaces the random search.  It is checked over
GF(p) (see ff.py), where a full rank is a proof for generic data.

Reduction
---------
1. Eigen-basis reduction.  Dropping the l-row of e moves column e from E to G;
   the template stays valid iff Top[:, e] is independent of the span of the
   current G columns.  That is a matroid, so one greedy pass gives the minimum
   eigenvalue-problem size |E| = |cols| - rank(Top) for this row set, without
   any search.  The order only decides *which* monomials stay in E (we keep
   low degree ones).
2. Template reduction.  Top rows are removed in a fixed order whenever the
   exact rank condition still holds; columns touched only by removed rows
   disappear with them.  Finally a square A2 is extracted.
Both steps are repeated until (|E|, |G|) stops decreasing.
"""
from dataclasses import dataclass, field
from itertools import combinations, product

import numpy as np

from . import ff
from .polytope import Polytope, minkowski


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _unit(n, i):
    return tuple(1 if j == i else 0 for j in range(n))


@dataclass
class Template:
    h: int                      # hidden variable index (0-based)
    top: list                   # [(k, t)]  rows t * f_k
    E: list                     # eigen basis monomials (l-rows)
    info: dict = field(default_factory=dict)

    @property
    def size(self):
        return (len(self.E), len(self.G))

    # filled by TemplateBuilder.finalize
    G: list = field(default_factory=list)
    recovery: dict = field(default_factory=dict)


class TemplateBuilder:
    def __init__(self, problem, h, seed=0, verbose=False, weights=None, eigen_target=1):
        self.pb = problem
        # stop the eigen basis reduction at this size (trade a bigger eigenvalue
        # problem for a smaller elimination A2); 1 = reduce to the minimum
        self.eigen_target = max(1, eigen_target or 1)
        # removal priority of monomials: heavier (w.m) monomials are removed first
        self.w = tuple(weights) if weights else (1,) * problem.n
        self.n = problem.n
        self.h = h
        self.eh = _unit(self.n, h)
        self.vals = problem.instance_mod_p(seed)
        self.verbose = verbose
        self._rowcache = {}

    def _prio(self, m):
        """Sort key: monomials that should be removed first come first."""
        return (-sum(w * x for w, x in zip(self.w, m)), -sum(m), tuple(-x for x in m))

    # ------------------------------------------------------------------ rows
    def row(self, r):
        """{column monomial: value mod p} of top row r = (k, t)."""
        if r not in self._rowcache:
            k, t = r
            self._rowcache[r] = {_add(t, a): self.vals[c] for a, c in self.pb.coeff_ids[k].items()}
        return self._rowcache[r]

    def columns(self, top, E):
        cols = set(E) | {_add(e, self.eh) for e in E}
        for r in top:
            cols.update(self.row(r))
        return cols

    def matrix(self, rows, cols):
        idx = {c: j for j, c in enumerate(cols)}
        M = np.zeros((len(rows), len(cols)))
        for i, r in enumerate(rows):
            for c, v in self.row(r).items():
                j = idx.get(c)
                if j is not None:
                    M[i, j] = v
        return M

    def top_column(self, rows, c):
        return np.array([self.row(r).get(c, 0) for r in rows], dtype=float)

    # ------------------------------------------------------------- checks
    def recovery_pairs(self, E, cols):
        """For each unhidden variable i: list of (e, e + e_i) with both monomials available."""
        Es = set(E)
        out = {}
        for i in range(self.n):
            if i == self.h:
                continue
            ei = _unit(self.n, i)
            out[i] = [(e, _add(e, ei)) for e in sorted(Es, key=lambda m: (sum(m), m)) if _add(e, ei) in cols]
        return out

    def recoverable(self, E, cols):
        return all(self.recovery_pairs(E, cols).values())

    def feasible(self, top, E):
        cols = self.columns(top, E)
        G = sorted(cols - set(E), key=lambda m: (sum(m), m))
        if len(top) < len(G) or not E:
            return False
        if not self.recoverable(E, cols):
            return False
        return ff.rank(self.matrix(top, G)) == len(G)

    # -------------------------------------------------------- construction
    def from_basis(self, B):
        Bs = set(B)
        top = []
        for k, A in enumerate(self.pb.supports):
            T = set()
            for b in B:
                t = _sub(b, A[0])
                if min(t) >= 0 and all(_add(t, a) in Bs for a in A):
                    T.add(t)
            top += [(k, t) for t in sorted(T)]
        E = sorted(t for t in B if _add(t, self.eh) in Bs)
        return top, E

    def repair(self, top, E, max_iter=None):
        """Make Top[:, G] full column rank by deleting the rows of non-eliminable columns.

        Columns are scanned from low to high degree, so the columns found to be
        dependent (and removed together with every row touching them) are the
        high degree ones.  Deterministic.
        """
        top, E = list(top), list(E)
        max_iter = max_iter or 10 * (len(top) + len(E))
        for _ in range(max_iter):
            if not top or not E:
                return None
            cols = self.columns(top, E)
            G = sorted(cols - set(E), key=lambda m: (sum(m), m))
            span = ff.IncrementalSpan(len(top))
            bad = None
            for g in G:
                if not span.add(self.top_column(top, g)):
                    bad = g
                    break
            if bad is None:
                if len(top) >= len(G) and self.recoverable(E, cols):
                    return top, E
                return None
            top = [r for r in top if bad not in self.row(r)]
            E = [e for e in E if e != bad and _add(e, self.eh) != bad]
        return None

    # ---------------------------------------------------------- reduction
    def reduce_eigen_basis(self, top, E):
        """Matroid-greedy removal of l-rows (minimal |E| for this row set)."""
        cols = self.columns(top, E)
        Ecur = set(E)
        G = sorted(cols - Ecur)
        span = ff.IncrementalSpan(len(top))
        for g in G:
            span.add(self.top_column(top, g))
        changed = True
        while changed:
            changed = False
            for e in sorted(Ecur, key=self._prio):
                if len(Ecur) <= self.eigen_target:
                    break
                v = self.top_column(top, e)
                touched_by_l = _sub(e, self.eh) in Ecur
                newE = Ecur - {e}
                if not np.any(v):
                    if touched_by_l:
                        continue
                    newcols = self.columns(top, newE)
                    if self.recoverable(newE, newcols):
                        Ecur = newE
                        changed = True
                    continue
                if span.contains(v):
                    continue
                newcols = self.columns(top, newE)
                if not self.recoverable(newE, newcols):
                    continue
                span.add(v)
                Ecur = newE
                changed = True
        return sorted(Ecur)

    def reduce_rows(self, top, E):
        """Greedy exact removal of top rows (high degree multipliers first)."""
        S = list(top)
        changed = True
        while changed:
            changed = False
            for r in sorted(S, key=lambda r: (self._prio(r[1]), -r[0])):
                S2 = [q for q in S if q != r]
                if self.feasible(S2, E):
                    S = S2
                    changed = True
        return S

    def reduce_columns(self, top, E):
        """Greedy exact removal of columns: drop a G column with *all* rows touching it.

        Works on the redundant (non-square) row set, where removing several rows
        at once can keep the rank condition although no single row removal does.
        Columns that l-rows need (e + e_h) cannot be dropped.
        """
        S = list(top)
        needed = {_add(e, self.eh) for e in E} | set(E)
        changed = True
        while changed:
            changed = False
            cols = self.columns(S, E)
            for c in sorted(cols - needed, key=self._prio):
                S2 = [r for r in S if c not in self.row(r)]
                if len(S2) < len(S) and self.feasible(S2, E):
                    S = S2
                    changed = True
        return S

    def squarify(self, top, E):
        cols = self.columns(top, E)
        G = sorted(cols - set(E), key=lambda m: (sum(m), m))
        rows = sorted(top, key=lambda r: (sum(r[1]), r[0], r[1]))
        keep = ff.independent_rows(self.matrix(rows, G))
        return [rows[i] for i in keep]

    def reduce(self, top, E, rows=True):
        best = None
        while True:
            E = self.reduce_eigen_basis(top, E)
            if rows:
                top = self.reduce_columns(top, E)
            top = self.squarify(top, E)
            if rows:
                top = self.reduce_rows(top, E)
                top = self.squarify(top, E)
            cols = self.columns(top, E)
            size = (len(E), len(cols) - len(E))
            if best is not None and size >= best[0]:
                return best[1], best[2]
            best = (size, top, E)

    def finalize(self, top, E, info=None):
        cols = self.columns(top, E)
        G = sorted(cols - set(E), key=lambda m: (sum(m), m))
        assert len(G) == len(top) and ff.rank(self.matrix(top, G)) == len(G), 'template is not valid'
        T = Template(self.h, list(top), list(E), dict(info or {}))
        T.G = G
        T.recovery = self.recovery_pairs(E, cols)
        return T


# ---------------------------------------------------------------------------
# Candidate generation and search
# ---------------------------------------------------------------------------

def _extra_polytopes(n, h):
    seg = Polytope([tuple([0] * n), _unit(n, h)])
    simplex = Polytope([tuple([0] * n)] + [_unit(n, i) for i in range(n)])
    return {'segment': seg, 'simplex': simplex}


def candidate_bases(problem, h, combs, shifts=None, extra='both'):
    """Deterministic list of (description, basis) from Minkowski sums of Newton polytopes."""
    n = problem.n
    polys = [Polytope(A) for A in problem.supports]
    extras = _extra_polytopes(n, h)
    if extra != 'both':
        extras = {extra: extras[extra]}
    shifts = shifts if shifts is not None else [s for s in product((-1, 0, 1), repeat=n)]
    seen = set()
    for comb in combs:
        for ename, ep in extras.items():
            parts = [polys[i - 1] if i <= len(polys) else ep for i in comb]
            if all(i <= len(polys) for i in comb) and ename != 'segment':
                continue  # extra polytope not used; avoid duplicates
            Q = minkowski(parts)
            for d in shifts:
                B = tuple(sorted(Q.lattice_points(d)))
                if not B or B in seen:
                    continue
                seen.add(B)
                yield {'comb': list(comb), 'extra': ename, 'shift': d}, list(B)


def default_combs(problem, cfg):
    m = len(problem.eqs)
    if cfg.get('combs'):
        return [list(c) for c in cfg['combs']]
    if cfg.get('poly_comb'):
        return [cfg['poly_comb']]
    sizes = cfg.get('size_of_combs') or [problem.n, problem.n + 1]
    combs = []
    for s in sizes:
        combs += [list(c) for c in combinations(range(1, m + 2), s)]
    return combs[: cfg.get('max_combs', 200)]


def search(problem, cfg=None, seed=0, keep=4, verbose=True, log=print):
    """Find a small valid template. Returns the best Template."""
    cfg = dict(problem.cfg, **(cfg or {}))
    hv = cfg.get('hidden_var', -1)
    hs = range(problem.n) if hv in (-1, None) else [hv - 1]
    pool = []
    for h in hs:
        tb = TemplateBuilder(problem, h, seed=seed, eigen_target=cfg.get('eigen_target'))
        ncand = 0
        for desc, B in candidate_bases(problem, h, default_combs(problem, cfg), extra=cfg.get('extra', 'both')):
            ncand += 1
            top, E = tb.from_basis(B)
            rep = tb.repair(top, E)
            if rep is None:
                continue
            top, E = rep
            E1 = tb.reduce_eigen_basis(top, E)
            top1 = tb.squarify(top, E1)
            ncols = len(tb.columns(top1, E1))
            score = (len(E1), ncols - len(E1))
            pool.append((score, h, desc, len(B), top, E))
        if verbose:
            log(f'hidden var a{h + 1}: {ncand} candidate bases, {sum(1 for p in pool if p[1] == h)} valid')
    if not pool:
        raise RuntimeError('no valid template found; try other poly_comb / size_of_combs')
    pool.sort(key=lambda p: (p[0], p[3]))
    n = problem.n
    weightings = [(1,) * n] + [tuple(2 if j == i else 1 for j in range(n)) for i in range(n)]
    if not cfg.get('try_orders', True):
        weightings = weightings[:1]
    best = None
    for score, h, desc, nB, top, E in pool[:keep]:
        for w in weightings:
            tb = TemplateBuilder(problem, h, seed=seed, weights=w, eigen_target=cfg.get('eigen_target'))
            top2, E2 = tb.reduce(top, E)
            T = tb.finalize(top2, E2, dict(desc, basis_size=nB, hidden_var=h + 1, weights=w,
                                            initial=(len(E), len(tb.columns(top, E)) - len(E)),
                                            after_eigen_reduction=score))
            if verbose:
                log(f'  {desc} |B|={nB} w={w}: eigen {score[0]}, template {score[1]} '
                    f'-> reduced eigen {T.size[0]}, template {T.size[1]}')
            if best is None or T.size < best.size:
                best = T
    return best
