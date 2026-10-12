#!/usr/bin/env python3
"""
construct_one_copy.py - Explicit construction (no search) of
CA(2q^3 - 1 + q^2 (q-2) ; 3, q^2+q+2, q)   for q a prime power (arithmetic in GF(q), see gf.py).   q=5 -> CA(324; 3, 32, 5)

1. RMS base: CA(2q^3-1; 3, q^2+q+1, q) (zero row, blocks B1, B2).
2. New column n: copy of column 0 in the zero row and in B1,
   copy of column 0 PLUS 1 in B2.
   -> every triple {i, j, n} with i, j != 0 is covered (like {i, j, 0}, up to a
      bijection in each block); for {0, w, n} the values with n - col0 in {0, 1}
      are already covered.
3. Added rows: for each a in GF(q), e in GF(q) minus {0, 1} (labelled 2..q-1), s in GF(q):
      col0 = a, n = a + e, all other columns = s.
   -> covers the missing (a, b, a+e) for every other column w.
   q * (q-2) * q rows.
Usage: python construct_one_copy.py [q]
"""
import sys
import numpy as np
from rms_base import build
from gf import GF

q = int(sys.argv[1]) if len(sys.argv) > 1 else 5
A, _ = build(q)
F = GF(q)
assert q >= 3, "the shifted-copy extension needs q >= 3"
n0 = 1 + (q**3 - 1)                       # rows: zero + B1, then B2
new = A[:, 0].copy()
new[n0:] = F.add[new[n0:], 1]
A = np.column_stack([A, new])
k = A.shape[1]
extra = []
for a in range(q):
    for e in range(2, q):
        for s in range(q):
            r = np.full(k, s); r[0] = a; r[-1] = F.add[a, e]
            extra.append(r)
A = np.vstack([A, np.array(extra)])
np.savetxt(f"CA_{A.shape[0]}_3_{k}_{q}.txt", A, fmt="%d")
print(f"CA({A.shape[0]}; 3, {k}, {q}) written")
