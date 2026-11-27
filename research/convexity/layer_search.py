"""Syzygy-aware lower bound on the size of a complete face layer (any dimension).

Given the supports S_1..S_m (lattice points in Z^d, coordinates of the face slice) of an
initial system in_w(F) with no torus zero, every complete set L (Theorem P1a: the top
w-layer of a valid template, w_h > 0) satisfies, independently of the coefficients:

  (a) every point of L lies in a row (translate t + S_i inside L);
  (b) L contains a translate of the anchor support (Lemma L for the face system);
  (c) |L| <= rank <= #rows - #syz, where #syz counts the Koszul syzygies
      in f_j * (x^t in f_i) - in f_i * (x^t in f_j), i < j, t + S_i + S_j inside L, with
      distinct largest rows (i, t + max S_j); such relations are independent (triangular);
  (d) P1a inside L: for every direction v, the v-top part L' of L is met by at least |L'|
      rows of L (a row meeting L' has its v-top part there; M[:, L'] must be injective).

A minimal complete L is row-connected, so we grow L from the anchor by adding rows that
meet it.  `layer_search_by_size` processes sets in increasing size (every set is generated
from a smaller one), so it stops at the first size where some set passes (a)-(d) and an
optional extra necessary condition `extra_ok`.  Returns (size, example) or None.

    from layer_search import layer_search
    layer_search([S1, S2, S3], anchor=0, extra=20)
"""
import time


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _translates(L, S):
    b0 = S[0]
    return [t for t in {_sub(c, b0) for c in L} if all(_add(t, b) in L for b in S)]


def _directions(d):
    import itertools
    from math import gcd
    out = []
    for v in itertools.product((-1, 0, 1), repeat=d):
        if any(v):
            out.append(v)
    return out


def make_slack(supports):
    supp = [sorted(map(tuple, S)) for S in supports]
    dirs = _directions(len(supp[0][0]))
    pairs = [(i, sorted({_add(a, b) for a in supp[i] for b in supp[j]}), max(supp[j]))
             for i in range(len(supp)) for j in range(i + 1, len(supp))]

    def slack(L):
        rows = sum(len(_translates(L, S)) for S in supp)
        lead = {(i, _add(t, mj)) for i, Sij, mj in pairs for t in _translates(L, Sij)}
        return rows - len(lead) - len(L), rows, len(lead)

    def sublayers_ok(L):
        R = [frozenset(_add(t, b) for b in S) for S in supp for t in _translates(L, S)]
        for v in dirs:
            top = max(sum(a * b for a, b in zip(v, p)) for p in L)
            Lp = {p for p in L if sum(a * b for a, b in zip(v, p)) == top}
            if sum(1 for r in R if r & Lp) < len(Lp):
                return False
        return True
    slack.sublayers_ok = sublayers_ok
    return supp, slack


def layer_search(supports, anchor, extra, verbose=True):
    supp, slack = make_slack(supports)
    A = frozenset(supp[anchor])
    limit = len(A) + extra
    level, seen, best, t0 = {A}, {A}, None, time.time()
    while level:
        nxt = set()
        for L in level:
            s = slack(L)
            if s[0] >= 0 and (best is None or len(L) < best[0]) and slack.sublayers_ok(L):
                best = (len(L), s, sorted(L))
            for S in supp:
                for t in {_sub(c, b) for c in L for b in S}:
                    r = frozenset(_add(t, b) for b in S)
                    if r <= L:
                        continue
                    L2 = L | r
                    if len(L2) <= limit and L2 not in seen:
                        seen.add(L2)
                        nxt.add(L2)
        level = nxt
        if verbose:
            print(f'  frontier {len(level):8d}, seen {len(seen):8d}, {time.time() - t0:6.0f}s', flush=True)
    return best, limit


def layer_search_by_size(supports, anchor, limit, extra_ok=None, verbose=True):
    """Smallest |L| <= limit passing (a)-(d) (and extra_ok), processing sets by size."""
    supp, slack = make_slack(supports)
    A = frozenset(supp[anchor])
    buckets = {len(A): {A}}
    seen = {A}
    t0 = time.time()
    for size in range(len(A), limit + 1):
        level = buckets.pop(size, set())
        hits = [L for L in level if slack(L)[0] >= 0 and slack.sublayers_ok(L)
                and (extra_ok is None or extra_ok(L))]
        if verbose:
            print(f'  |L| = {size}: {len(level)} sets, {len(hits)} pass, {time.time() - t0:.0f}s', flush=True)
        if hits:
            return size, sorted(hits[0])
        for L in level:
            for S in supp:
                for t in {_sub(c, b) for c in L for b in S}:
                    r = frozenset(_add(t, b) for b in S)
                    if r <= L:
                        continue
                    L2 = L | r
                    if len(L2) <= limit and L2 not in seen:
                        seen.add(L2)
                        buckets.setdefault(len(L2), set()).add(L2)
        seen -= level            # sets of this size can no longer be regenerated
    return None
