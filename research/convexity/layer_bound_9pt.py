"""Counting bound for the e3 face layer of 9pt two-radial (a3 hidden), via compression.

Theorem P1a: the top a3-layer L of any valid template lies in G and is complete for
in_{e3}(F).  Completeness needs #rows(L) >= |L|, where the rows are the translates of the
four initial supports (coordinates a1, a2, a4) inside L:

    S1' = S2' = {a1 + a2 <= 1, a1 + a2 + a4 <= 2}   (7 points)
    S3'       = {a1 + a2 + a4 <= 1}                  (4 points)
    S4'       = {a1 + a2 + a4 <= 3, a4 <= 2}          (19 points)

All four are needed (no proper subset of in_{e3}(F) is zero-free, face_layers.py), so L
contains a translate of S4'.  All supports are down-sets, so by the compression lemma
(lower_bounds.md) compressing L along each axis never decreases a translate count; the
minimum of |L| subject to the counting condition is therefore attained by a down-set
L of N^3 containing S4'.  We enumerate those by size.
"""
import itertools

S1 = {p for p in itertools.product(range(3), repeat=3) if p[0] + p[1] <= 1 and sum(p) <= 2}
S3 = {p for p in itertools.product(range(2), repeat=3) if sum(p) <= 1}
S4 = {p for p in itertools.product(range(4), repeat=3) if sum(p) <= 3 and p[2] <= 2}
SUPPORTS = [S1, S1, S3, S4]
assert [len(s) for s in SUPPORTS] == [7, 7, 4, 19]


def fits(L, S):
    return sum(1 for t in L if all((t[0] + a, t[1] + b, t[2] + c) in L for a, b, c in S))


def slack(L):
    return sum(fits(L, S) for S in SUPPORTS) - len(L)


def addable(L):
    """Points whose addition keeps L a down-set."""
    out = set()
    for p in L:
        for k in range(3):
            q = list(p); q[k] += 1; q = tuple(q)
            if q not in L and all(tuple(q[j] - (j == i) for j in range(3)) in L
                                  for i in range(3) if q[i] > 0):
                out.add(q)
    for k in range(3):
        q = tuple(int(j == k) for j in range(3))
        if q not in L and False:
            out.add(q)
    return out


def main(max_extra=30):
    start = frozenset(S4)
    level = {start}
    seen = {start}
    best = max(slack(start), -10**9)
    print(f'|L| = {len(start)}: 1 down-set, max slack {slack(start)}')
    for k in range(1, max_extra + 1):
        nxt = set()
        for L in level:
            for q in addable(L):
                L2 = L | {q}
                if L2 not in seen:
                    seen.add(L2)
                    nxt.add(L2)
        level = nxt
        sl = max(slack(L) for L in level)
        print(f'|L| = {len(start) + k}: {len(level)} down-sets, max #rows - |L| = {sl}', flush=True)
        if sl >= 0:
            ex = max(level, key=slack)
            print('first feasible size', len(start) + k, 'example', sorted(ex))
            return len(start) + k


if __name__ == '__main__':
    main()
