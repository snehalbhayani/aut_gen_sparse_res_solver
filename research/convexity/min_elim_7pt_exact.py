"""Exact (1b) search for the 7pt focal+radial (one-sided) problem, a3 hidden.

Claim: every valid template with the minimal eigen size 19 has |G| >= 22, i.e. |C| >= 41.
symgen's 22 + 19 template is therefore optimal.

Proof structure (see research/proof/lower_bounds.md, Theorem P1d):

1. Lemma L + support_bound.py: the used equations include f3, f4 or f5, so C contains a
   translate of S3 (35 points), S4 or S5 (38 points).  Translate it to the origin.
2. Validity at |E| = 19 forces rank M_C >= |C| - 19 (the kernel injects into K^E), hence
   #rows(C) >= |C| - 19, where #rows(C) counts ALL translates t + S_i inside C.
3. Points of C hit by no fitting row can be deleted: the rows stay the same and |C| drops,
   so the inequality of step 2 only gets harder.  Hence it suffices to treat C = anchor + X
   where X is a union of exteriors (row minus anchor) of rows fitting in C.
4. |C| <= 40 bounds |X| (<= 5 for S3, <= 2 for S4/S5).  We enumerate every such union and
   count rows.  Every candidate has #rows(C) < |C| - 19.

Steps 2-4 are pure counting, so the bound holds for every data instance for which steps 1
and "19 is the minimal eigen size" hold (generic data).  The numeric part below only
re-checks the latter: 19 solutions with pairwise distinct a3 (Theorem 2 then gives |E| >= 19).
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from symgen.problem import load_matlab_problem
from symgen.template import search
from symgen.numeric import NumericSolver

NAME, H, EIG, MAXC = 'problem_relpose_7p_fr_1s_partial_elim', 2, 19, 40

pb = load_matlab_problem(f'problems/{NAME}.m')
S = [frozenset(map(tuple, s)) for s in pb.supports]
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))

def fitting_rows(C):
    """All rows x^t f_i with t + S_i inside C (t ranges over C - S_i)."""
    rows = []
    for i, Si in enumerate(S):
        b0 = next(iter(Si))
        for t in {sub(c, b0) for c in C}:
            if all(add(t, b) in C for b in Si):
                rows.append((i, t))
    return rows


# ---- numeric check: 19 solutions with distinct a3 (so the minimal eigen size is 19) ----
T = search(pb, keep=1, verbose=False, cfg={'try_orders': False})
assert T.h == H and len(T.E) == EIG, (T.h, len(T.E))
sol = NumericSolver(pb, T)
d = np.random.default_rng(1).standard_normal(len(pb.params))
X = sol.solve(d)
res = sol.residuals(X, d)
print(f'{NAME}: {len(X)} solutions at random data, max residual {max(res):.1e}, '
      f'min |a3| gap {min(abs(a - b) for a, b in itertools.combinations(X[:, H], 2)):.2e}')

# ---- exhaustive search (pure counting) -------------------------------------------------
total = 0
for anc in (2, 3, 4):
    A = S[anc]
    budget = MAXC - len(A)
    ext = set()
    for i, Si in enumerate(S):
        for t in {sub(a, b) for a in A for b in Si}:
            e = frozenset(add(t, b) for b in Si) - A
            if 0 < len(e) <= budget:
                ext.add(e)
    unions, frontier = {frozenset()}, {frozenset()}
    while frontier:
        frontier = {U | e for U in frontier for e in ext
                    if len(U | e) <= budget and (U | e) not in unions}
        unions |= frontier
    slack = []
    for Xs in unions:
        C = A | Xs
        slack.append(len(fitting_rows(C)) - (len(C) - EIG))
        total += 1
    print(f'  anchor S{anc+1} (|S| = {len(A)}, up to {budget} extra points): {len(ext)} row exteriors, '
          f'{len(unions)} candidate column sets; max of #rows - (|C| - {EIG}) = {max(slack)}'
          + ('  -> none feasible' if max(slack) < 0 else '  -> FEASIBLE CANDIDATE'))
print(f'checked {total} column sets; |C| >= {MAXC + 1}, so |G| >= {MAXC + 1 - EIG} at eigen size {EIG}')
