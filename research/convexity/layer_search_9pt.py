"""Syzygy-aware exact search for the e3 face layer of 9pt two-radial (a3 hidden).

The top a3-layer L of a valid template is complete for in_{e3}(F) (Theorem P1a): the rows
x^t in_{e3}(f_i) inside L have rank |L|.  Necessary conditions, all data-independent:

  (a) every point of L lies in such a row (a zero column breaks injectivity);
  (b) L contains a translate of S4' (no proper subset of in_{e3}(F) is zero-free);
  (c) rank <= #rows - #syz, where #syz counts Koszul syzygies with distinct leading rows.
      For i < j and t + S_i' + S_j' inside L, the identity
          in f_j * (x^t in f_i) - in f_i * (x^t in f_j) = 0
      is a linear relation among rows of L.  Ordering rows by (equation index, multiplier),
      its largest row is (i, t + max S_j'); relations with distinct largest rows are
      independent (triangular).

A minimal complete L is row-connected (its block-diagonal pieces are complete on their own),
so we grow L from S4' by adding rows that meet it, and test (c) on every set up to the size
limit.  If none passes, lambda_{e3} > limit.

    python3 research/convexity/layer_search_9pt.py 14
"""
import sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from layer_bound_9pt import SUPPORTS, S4

add = lambda a, b: (a[0] + b[0], a[1] + b[1], a[2] + b[2])
sub = lambda a, b: (a[0] - b[0], a[1] - b[1], a[2] - b[2])
SUPP = [sorted(S) for S in SUPPORTS]
PAIRS = [(i, j, sorted({add(a, b) for a in SUPP[i] for b in SUPP[j]}), max(SUPP[j]))
         for i in range(4) for j in range(i + 1, 4)]


def translates(L, S):
    b0 = S[0]
    return [t for t in {sub(c, b0) for c in L} if all(add(t, b) in L for b in S)]


def slack(L):
    rows = sum(len(translates(L, S)) for S in SUPP)
    lead = set()
    for i, j, Sij, mj in PAIRS:
        for t in translates(L, Sij):
            lead.add((i, add(t, mj)))
    return rows - len(lead) - len(L), rows, len(lead)


def main(extra):
    A = frozenset(S4)
    limit = len(A) + extra
    level, seen = {A}, {A}
    best = None
    t0 = time.time()
    while level:
        nxt = set()
        for L in level:
            sl = slack(L)
            if sl[0] >= 0 and (best is None or len(L) < best[0]):
                best = (len(L), sl, sorted(L))
            for S in SUPP:                       # rows meeting L
                for t in {sub(c, b) for c in L for b in S}:
                    r = frozenset(add(t, b) for b in S)
                    if r <= L:
                        continue
                    L2 = L | r
                    if len(L2) <= limit and L2 not in seen:
                        seen.add(L2)
                        nxt.add(L2)
        level = nxt
        print(f'  frontier {len(level):8d}, seen {len(seen):8d}, {time.time()-t0:6.0f}s', flush=True)
    if best:
        print(f'smallest set passing the syzygy-aware count: |L| = {best[0]} (slack, rows, syz) = {best[1]}')
    else:
        print(f'no row-connected set with |L| <= {limit} passes: lambda_e3 >= {limit + 1}')
    return best


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
