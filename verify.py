#!/usr/bin/env python3
"""
verify.py - Independent covering array checker
Usage: python verify.py file.txt t v
Checks that for EVERY choice of t columns, all v^t value combinations
appear in at least one row. Deliberately simple code (brute force).
"""
import sys
import itertools
import numpy as np


def check(A, t, v):
    N, k = A.shape
    assert A.min() >= 0 and A.max() < v, "values outside 0..v-1"
    weights = v ** np.arange(t - 1, -1, -1)
    missing = 0
    for cols in itertools.combinations(range(k), t):
        codes = A[:, cols] @ weights
        missing += v**t - len(np.unique(codes))
    return N, k, missing


if __name__ == "__main__":
    f, t, v = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    A = np.loadtxt(f, dtype=int, ndmin=2)
    N, k, m = check(A, t, v)
    print(f"{f}: N={N} rows, k={k} columns, strength t={t}, v={v} -> "
          + ("VALID (all combinations covered)" if m == 0 else f"INVALID ({m} combinations missing)"))
