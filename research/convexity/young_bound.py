"""Rigorous P1b bound in 2-D via compression: enumerate Young diagrams.

Compression lemma: for down-set supports S_i (e.g. full boxes), column/row compressions of
any finite C do not decrease |C - S_i| (number of fitting translates).  Hence
  min{ |C| : |C| - sum_i |C - S_i| <= N }  is attained by a Young diagram,
and every valid template with eigen size N has |G| >= L(N) - N.
"""
from functools import lru_cache


def partitions(s, maxpart=None):
    maxpart = s if maxpart is None else maxpart
    if s == 0:
        yield ()
        return
    for k in range(min(s, maxpart), 0, -1):
        for rest in partitions(s - k, k):
            yield (k,) + rest


def fits(lam, box):
    """# translates of the full box [0,p]x[0,q] inside the Young diagram with row lengths lam."""
    p, q = box                      # widths in x (columns) and y (rows), as point counts - 1
    # down-set: t fits iff its top-right corner t + (p, q) is a cell
    return sum(max(0, lam[y + q] - p) for y in range(len(lam) - q)) if len(lam) > q else 0


def L(N, boxes, smax=60):
    for s in range(1, smax + 1):
        best = None
        for lam in partitions(s):
            d = s - sum(fits(lam, b) for b in boxes)
            if best is None or d < best[0]:
                best = (d, lam)
        if best[0] <= N:
            return s, best
    return None


if __name__ == '__main__':
    # thesis Ex. 4.2.1: f1 support = full 4x4 box (p=q=3), f2 support = full 2x3 box (p=1, q=2)
    boxes = [(3, 3), (1, 2)]
    for N in (9, 15):
        s, (d, lam) = L(N, boxes)
        print(f'N = {N}: smallest column set with deficiency <= N has {s} cells '
              f'(Young diagram rows {lam}, deficiency {d})  =>  every valid template with eigen size {N} '
              f'has |G| >= {s - N}')
