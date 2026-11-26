"""Exact linear algebra over a prime field GF(p).

All template decisions (is this column eliminable? can this row be dropped?)
are made with exact rank computations over GF(p), with the coefficient
polynomials evaluated at a seeded random point.  A *full* rank observed mod p
is a proof that the corresponding minor is a non-zero polynomial in the data
coefficients, i.e. the property holds for generic data.  (A rank *deficiency*
can, with probability ~deg/p, be a false alarm; that only ever makes a
template slightly less reduced, never invalid.)

Arithmetic uses float64 with p < 2^21: a single product is < 2^42 and a dot
product of up to 2048 terms stays below 2^53, so everything is exact.
"""
import numpy as np

# Largest prime below 2^21.
P = 2097143
_CHUNK = 1024


def inv(a, p=P):
    return pow(int(a) % p, p - 2, p)


def rank(A, p=P):
    """Rank of an integer matrix (entries already reduced mod p)."""
    A = np.array(A, dtype=np.float64, copy=True)
    if A.size == 0:
        return 0
    m, n = A.shape
    r = 0
    for c in range(n):
        if r == m:
            break
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        i = r + nz[0]
        if i != r:
            A[[r, i]] = A[[i, r]]
        A[r] = np.mod(A[r] * inv(A[r, c], p), p)
        below = r + 1 + np.flatnonzero(A[r + 1:, c])
        if below.size:
            A[below] = np.mod(A[below] - np.outer(A[below, c], A[r]), p)
        r += 1
    return r


def independent_rows(A, p=P):
    """Indices of a maximal set of linearly independent rows (greedy in row order)."""
    return IncrementalSpan(A.shape[1] if A.ndim == 2 else 0, p).extend(A)


class IncrementalSpan:
    """Fully reduced echelon basis of a growing subspace of GF(p)^n.

    `add(v)` returns True (and enlarges the span) iff v is not already in it.
    Adding vectors in a fixed order and keeping the independent ones is the
    matroid greedy algorithm, so the number accepted is order independent.
    """

    def __init__(self, n, p=P):
        self.n, self.p = n, p
        self.B = np.zeros((0, n))
        self.piv = []

    def reduce(self, v):
        v = np.mod(np.asarray(v, dtype=np.float64), self.p)
        for s in range(0, len(self.piv), _CHUNK):
            piv = self.piv[s:s + _CHUNK]
            v = np.mod(v - np.mod(v[piv] @ self.B[s:s + _CHUNK], self.p), self.p)
        return v

    def contains(self, v):
        return not np.any(self.reduce(v))

    def add(self, v):
        v = self.reduce(v)
        nz = np.flatnonzero(v)
        if nz.size == 0:
            return False
        q = nz[0]
        v = np.mod(v * inv(v[q], self.p), self.p)
        if self.piv:
            self.B = np.mod(self.B - np.outer(self.B[:, q], v), self.p)
        self.B = np.vstack([self.B, v])
        self.piv.append(q)
        return True

    def extend(self, rows):
        return [i for i, v in enumerate(rows) if self.add(v)]
