#!/usr/bin/env python3
"""
rms_base.py - Rebuilds CA(2q^3-1; 3, q^2+q+1, q) of Raaphorst-Moura-Stevens (2014)
==================================================================================

Construction (q a prime power, arithmetic in GF(q), see gf.py):
  - f: primitive polynomial of degree 3 over GF(q), alpha its root in GF(q^3).
  - Column j (j = 0..k-1, k = q^2+q+1) <-> point p_j = alpha^j (coordinates in GF(q)^3).
  - Block 1: rows x (non-zero vectors of GF(q)^3), entry <x, p_j>  (= windows of the m-sequence).
  - Block 2: same, with the column order reversed (p_{k-1-j}): reciprocal polynomial.
  - Plus the all-zero row.
Three columns that are collinear in block 1 are not collinear in block 2.
Total: 2(q^3-1) + 1 = 2q^3 - 1 rows.

For prime q the output is identical to the earlier prime-only version.
Usage: python rms_base.py [q]     ->  base_q{q}.txt
"""
import sys
import itertools
import numpy as np
from gf import GF


def mulmod(a, b, f, F):
    """Product of polynomials over GF(q) (coefficient lists, degree < 3) modulo
    the monic cubic x^3 + f2 x^2 + f1 x + f0."""
    A, M = F.add, F.mul
    r = [0] * 5
    for i in range(3):
        for j in range(3):
            r[i + j] = A[r[i + j], M[a[i], b[j]]]
    for d in (4, 3):                     # x^d = x^(d-3) * x^3, x^3 = -(f0 + f1 x + f2 x^2)
        c = r[d]
        if c:
            r[d] = 0
            for i in range(3):
                r[d - 3 + i] = A[r[d - 3 + i], F.neg[M[c, f[i]]]]
    return [int(v) for v in r[:3]]


def powers(f, F):
    """List of alpha^j, j = 0..q^3-2; None if f is not primitive."""
    q = F.q
    x, cur, out, seen = [0, 1, 0], [1, 0, 0], [], set()
    for _ in range(q**3 - 1):
        t = tuple(cur)
        if t in seen:
            return None
        seen.add(t); out.append(t)
        cur = mulmod(cur, x, f, F)
    return out if cur == [1, 0, 0] else None


def primitive_polynomial(F):
    for f in itertools.product(range(F.q), repeat=3):
        if f[0] == 0:
            continue
        p = powers(list(f), F)
        if p:
            return list(f), p
    raise ValueError("no primitive polynomial found")


def build(q):
    F = GF(q)
    f, pw = primitive_polynomial(F)
    k = q * q + q + 1
    P = np.array(pw[:k])                                   # points of PG(2,q) (Singer cycle)
    X = np.array([x for x in itertools.product(range(q), repeat=3) if any(x)])
    B1 = F.dot(X, P)
    B2 = F.dot(X, P[::-1])
    return np.vstack([np.zeros((1, k), dtype=int), B1, B2]), f


if __name__ == "__main__":
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    A, f = build(q)
    np.savetxt(f"base_q{q}.txt", A, fmt="%d")
    print(f"q={q}, f(x) = x^3 + {f[2]}x^2 + {f[1]}x + {f[0]}: {A.shape[0]} rows x {A.shape[1]} columns")
