#!/usr/bin/env python3
"""
rms_base.py - Rebuilds CA(2q^3-1; 3, q^2+q+1, q) of Raaphorst-Moura-Stevens (2014)
==================================================================================

Construction (q prime):
  - f: primitive polynomial of degree 3 over GF(q), alpha its root in GF(q^3).
  - Column j (j = 0..k-1, k = q^2+q+1) <-> point p_j = alpha^j (coordinates in GF(q)^3).
  - Block 1: rows x (non-zero vectors of GF(q)^3), entry <x, p_j>  (= windows of the m-sequence).
  - Block 2: same, with the column order reversed (p_{k-1-j}): reciprocal polynomial.
  - Plus the all-zero row.
Three columns that are collinear in block 1 are not collinear in block 2.
Total: 2(q^3-1) + 1 = 2q^3 - 1 rows.

Usage: python rms_base.py [q]     ->  base_q{q}.txt
"""
import sys
import itertools
import numpy as np


def mulmod(a, b, f, q):
    """Product of polynomials (coefficient lists, degree < 3) modulo monic f of degree 3."""
    r = [0] * 5
    for i in range(3):
        for j in range(3):
            r[i + j] = (r[i + j] + a[i] * b[j]) % q
    for d in (4, 3):                     # reduce x^d using x^3 = -(f0 + f1 x + f2 x^2)
        c = r[d]
        if c:
            r[d] = 0
            for i in range(3):
                r[d - 3 + i] = (r[d - 3 + i] - c * f[i]) % q
    return r[:3]


def powers(f, q):
    """List of alpha^j, j = 0..q^3-2; None if f is not primitive."""
    x, cur, out = [0, 1, 0], [1, 0, 0], []
    for _ in range(q**3 - 1):
        out.append(tuple(cur))
        cur = mulmod(cur, x, f, q)
    return out if cur == [1, 0, 0] and len(set(out)) == q**3 - 1 else None


def primitive_polynomial(q):
    for f in itertools.product(range(q), repeat=3):
        if f[0] == 0:
            continue
        p = powers(list(f), q)
        if p:
            return list(f), p
    raise ValueError("no primitive polynomial found")


def build(q):
    f, pw = primitive_polynomial(q)
    k = q * q + q + 1
    P = np.array(pw[:k])                                   # points of PG(2,q) (Singer cycle)
    X = np.array([x for x in itertools.product(range(q), repeat=3) if any(x)])
    B1 = (X @ P.T) % q
    B2 = (X @ P[::-1].T) % q
    return np.vstack([np.zeros((1, k), dtype=int), B1, B2]), f


if __name__ == "__main__":
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    A, f = build(q)
    np.savetxt(f"base_q{q}.txt", A, fmt="%d")
    print(f"q={q}, f(x) = x^3 + {f[2]}x^2 + {f[1]}x + {f[0]}: {A.shape[0]} rows x {A.shape[1]} columns")
