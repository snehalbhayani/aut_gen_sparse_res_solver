"""Close the (1b) gap on the thesis Ex. 4.2.1 instance (x1 hidden, eigen size 9) by exhaustion.

Rigorous reductions (see lower_bounds.md):
  * Lemma L: both f1 and f2 must be used, so C contains a translate of supp f1 = 4x4 box;
    translate it to [0,3]^2 (validity is translation invariant).
  * P1a: the rightmost column of C (x1 = D) lies in G and holds >= 5 points.
    If D > 3, that column misses the box, |C| >= 21, |G| >= 12.  So D = 3.
  * |G| = |C| - 9 <= 10 needs |C| <= 19: the box plus <= 3 further points with x1 <= 3.
Every such C is checked with EXACT rational ranks.  No valid one => min |G| = 11 (achieved).
"""
import itertools
from fractions import Fraction

f1 = {(3,3):-45,(3,2):-9,(2,3):375,(3,1):765,(2,2):228,(1,3):-660,(3,0):626,(2,1):-307,(1,2):-366,
      (0,3):-180,(2,0):150,(1,1):-1052,(0,2):464,(1,0):-1096,(0,1):1156,(0,0):156}
f2 = {(1,2):3,(1,1):10,(0,2):-9,(1,0):3,(0,1):-4,(0,0):-8}
EQS = [f1, f2]


def rank_exact(rows, cols):
    M = [[Fraction(r.get(c, 0)) for c in cols] for r in rows]
    rk, ncol = 0, len(cols)
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                fac = M[i][c] / M[rk][c]
                M[i] = [a - fac * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def rows_in(C):
    Cs = set(C); out = []
    for f in EQS:
        base = next(iter(f))
        for c in C:
            t = (c[0] - base[0], c[1] - base[1])
            if all((t[0] + a, t[1] + b) in Cs for (a, b) in f):
                out.append({(t[0] + a, t[1] + b): v for (a, b), v in f.items()})
    return out


def check(C):
    """valid template with eigen size 9 on column set C (x1 hidden)?"""
    Cs = set(C)
    rows = rows_in(C)
    cols = sorted(C)
    rk = rank_exact(rows, cols)
    if len(C) - rk > 9:
        return False                       # corank > 9: eigen size cannot be 9
    bd = [c for c in C if (c[0] + 1, c[1]) not in Cs]
    return rank_exact(rows, bd) == len(bd)   # boundary eliminable => min eigen = corank (= 9)


box = [(a, b) for a in range(4) for b in range(4)]
region = [(a, b) for a in range(-4, 4) for b in range(-4, 8) if (a, b) not in box]
checked = found = 0
for k in range(1, 4):
    for extra in itertools.combinations(region, k):
        C = box + list(extra)
        col3 = [p for p in C if p[0] == 3]
        if len(col3) < 5:
            continue                        # P1a: the rightmost column needs >= 5 points
        checked += 1
        if check(C):
            found += 1
            print('VALID with |G| =', len(C) - 9, ':', sorted(extra))
print(f'checked {checked} column sets with |C| <= 19 (exact rational ranks); valid ones: {found}')
print('=> minimal |G| at eigen size 9 is 11' if found == 0 else '=> gap not closed')
print('sanity: the 4x5 box [0,3]x[0,4] (|C| = 20):', check([(a, b) for a in range(4) for b in range(5)]))
